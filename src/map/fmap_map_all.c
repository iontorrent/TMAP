#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <config.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#include <unistd.h>
#endif
#include <unistd.h>
#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "../util/fmap_progress.h"
#include "../util/fmap_sam.h"
#include "../util/fmap_sort.h"
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt_gen.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_bwt_match.h"
#include "../index/fmap_sa.h"
#include "../io/fmap_seq_io.h"
#include "../server/fmap_shm.h"
#include "fmap_map1.h"
#include "fmap_map1_aux.h"
#include "fmap_map2.h"
#include "fmap_map2_aux.h"
#include "fmap_map3.h"
#include "fmap_map3_aux.h"
#include "fmap_map_util.h"
#include "fmap_map_all.h"

/* Notes:
   We could avoid some computation give we know which algorithms to 
   run. This includes stacks, memory pools, as well as not loading 
   in all reference data.
   */

// sort by min-seqid, min-position, max-score
#define __fmap_map_all_hit_sort_lt(a, b) ( ((a).seqid < (b).seqid \
                                            || ( (a).seqid == (b).seqid && (a).pos < (b).pos ) \
                                            || ( (a).seqid == (b).seqid && (a).pos == (b).pos && (a).score < (b).score )) \
                                          ? 1 : 0 )

FMAP_SORT_INIT(fmap_map_all_hit_t, fmap_map_all_hit_t, __fmap_map_all_hit_sort_lt)

#ifdef HAVE_LIBPTHREAD
     static pthread_mutex_t fmap_map_all_read_lock = PTHREAD_MUTEX_INITIALIZER;
     static int32_t fmap_map_all_read_lock_low = 0;
#define FMAP_MAP_ALL_THREAD_BLOCK_SIZE 512
#endif


     static char *algos[6] = {
         "dummy",
         "map1", // 0x1 
         "map2", // 0x2 
         "dummy",
         "map3", // 0x4
         "dummy"
     };

static inline fmap_map_all_aln_t *
fmap_map_all_aln_init()
{
  return fmap_calloc(1, sizeof(fmap_map_all_aln_t), "return");
}

static inline void
fmap_map_all_aln_destroy(fmap_map_all_aln_t *aln)
{
  int32_t i;
  for(i=0;i<aln->n;i++) {
      free(aln->hits[i].cigar);
  }
  free(aln->hits);
  free(aln);
}

static inline void
fmap_map_all_aln_realloc(fmap_map_all_aln_t *aln, int32_t n)
{
  int32_t i;
  for(i=n;i<aln->n;i++) {
      free(aln->hits[i].cigar);
  }
  aln->hits = fmap_realloc(aln->hits, sizeof(fmap_map_all_hit_t) * n, "aln->hits");
  aln->n = n;
}

static inline void
fmap_map_all_aln_mapq(fmap_map_all_aln_t *aln, fmap_map_all_opt_t *opt)
{
  int32_t i;
  int32_t n_best = 0, n_best_subo = 0;
  int32_t best_score, cur_score, best_subo;
  int32_t mapq;
  int32_t stage = -1;
  int32_t algo_id = FMAP_MAP_ALL_ALGO_NONE; 

  // estimate mapping quality TODO: this needs to be refined
  best_score = INT32_MIN;
  best_subo = INT32_MIN+1;
  n_best = n_best_subo = 0;
  for(i=0;i<aln->n;i++) {
      cur_score = aln->hits[i].score;
      if(best_score < cur_score) {
          // save sub-optimal
          best_subo = best_score;
          n_best_subo = n_best;
          // update
          best_score = cur_score;
          n_best = 1;
          stage = (algo_id == FMAP_MAP_ALL_ALGO_NONE) ? aln->hits[i].algo_stage : -1;
          algo_id = (algo_id == FMAP_MAP_ALL_ALGO_NONE) ? aln->hits[i].algo_id : -1;
      }
      else if(cur_score == best_score) { // qual
          n_best++;
      }
      else {
          if(best_subo < cur_score) {
              best_subo = cur_score;
              n_best_subo = 1;
          }
          else if(best_subo == cur_score) {
              n_best_subo++;
          }
      }
      if(FMAP_MAP_ALL_ALGO_MAP2 == aln->hits[i].algo_id
         || FMAP_MAP_ALL_ALGO_MAP3 == aln->hits[i].algo_id) {
          cur_score = aln->hits[i].score_subo;
          if(INT32_MIN == cur_score) {
              // ignore
          }
          else if(best_subo < cur_score) {
              best_subo = cur_score;
              n_best_subo=1;
          } 
          else if(best_subo == cur_score) {
              n_best_subo++;
          }
      }
  }
  if(1 < n_best || best_score <= best_subo) {
      mapq = 0;
  }
  else {
      if(0 == n_best_subo) {
          n_best_subo = 1;
          switch(algo_id) {
            case FMAP_MAP_ALL_ALGO_MAP1:
              // what is the best value for map1
              if(0 < opt->opt_map1[stage]->seed_length) {
                  best_subo = opt->score_match * (opt->opt_map1[stage]->seed_length - opt->opt_map1[stage]->seed_max_mm);
              }
              else {
                  best_subo = 0;
              }
              break;
            case FMAP_MAP_ALL_ALGO_MAP2:
              best_subo = opt->opt_map2[stage]->score_thr; break;
            case FMAP_MAP_ALL_ALGO_MAP3:
              best_subo = opt->opt_map3[stage]->score_thr; break;
            default:
              best_subo = 0; break;
          }
      }
      /*
      fprintf(stderr, "n_best=%d n_best_subo=%d\n",
              n_best, n_best_subo);
      fprintf(stderr, "best_score=%d best_subo=%d\n",
              best_score, best_subo);
      */
      mapq = (int32_t)((n_best / (1.0 * n_best_subo)) * (best_score - best_subo) * (250.0 / best_score + 0.03 / opt->score_match) + .499);
      if(mapq > 250) mapq = 250;
      if(mapq <= 0) mapq = 1;
  }
  for(i=0;i<aln->n;i++) {
      cur_score = aln->hits[i].score;
      if(cur_score == best_score) {
          aln->hits[i].mapq = mapq;
      }
      else {
          aln->hits[i].mapq = 0;
      }
  }
}

static inline void
fmap_map_all_aln_filter(fmap_map_all_aln_t *aln, fmap_map_all_opt_t *opt, int32_t algo_id)
{
  int32_t i, j, k;
  int32_t n_best = 0;
  int32_t best_score, cur_score;

  if(FMAP_MAP_UTIL_ALN_MODE_ALL == opt->aln_output_mode
     || aln->n <= 1) {
      return;
  }

  best_score = INT32_MIN;
  n_best = 0;
  for(i=0;i<aln->n;i++) {
      if(FMAP_MAP_ALL_ALGO_NONE == algo_id
         || aln->hits[i].algo_id == algo_id) {
          cur_score = aln->hits[i].score;
          if(best_score < cur_score) {
              best_score = cur_score;
              n_best = 1;
          }
          else if(!(cur_score < best_score)) { // equal
              n_best++;
          }
      }
  }

  // copy to the front
  if(n_best < aln->n) {
      for(i=j=0;i<aln->n;i++) {
          if(FMAP_MAP_ALL_ALGO_NONE == algo_id
             || aln->hits[i].algo_id == algo_id) {
              cur_score = aln->hits[i].score;
              if(cur_score < best_score) { // not the best
                  free(aln->hits[i].cigar);
                  aln->hits[i].cigar = NULL;
                  aln->hits[i].n_cigar = 0;
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      aln->hits[j] = aln->hits[i];
                      aln->hits[i].cigar = NULL;
                  }
                  j++;
              }
          }
          else {
              if(j < i) { // copy if we are not on the same index
                  aln->hits[j] = aln->hits[i];
                  aln->hits[i].cigar = NULL;
              }
              j++;
          }
      }
      // reallocate
      fmap_map_all_aln_realloc(aln, j);
  }

  if(FMAP_MAP_UTIL_ALN_MODE_UNIQ_BEST == opt->aln_output_mode) {
      if(1 < n_best) { // there can only be one
          if(FMAP_MAP_ALL_ALGO_NONE == algo_id) {
              fmap_map_all_aln_realloc(aln, 0);
          }
          else {
              // get rid of all of them
              for(i=j=0;i<aln->n;i++) {
                  if(aln->hits[i].algo_id == algo_id) {
                      free(aln->hits[i].cigar);
                      aln->hits[i].cigar=NULL;
                  }
                  else {
                      if(j < i) { // copy if we are not on the same index
                          aln->hits[j] = aln->hits[i];
                          aln->hits[i].cigar = NULL;
                      }
                      j++;
                  }
              }
              fmap_map_all_aln_realloc(aln, j);
          }
      }
  }
  else if(FMAP_MAP_UTIL_ALN_MODE_RAND_BEST == opt->aln_output_mode) { // get a random
      int32_t r = (int32_t)(drand48() * n_best);

      // keep the rth one
      if(FMAP_MAP_ALL_ALGO_NONE == algo_id) {
          if(0 != r) {
              free(aln->hits[0].cigar);
              aln->hits[0] = aln->hits[r];
              aln->hits[r].cigar = NULL;
          }
          // reallocate
          fmap_map_all_aln_realloc(aln, 1);
      }
      else {
          // keep the rth one
          for(i=j=k=0;i<aln->n;i++) {
              if(aln->hits[i].algo_id == algo_id) {
                  if(k == r) { // keep
                      if(j < i) { // copy if we are not on the same index
                          aln->hits[j] = aln->hits[i];
                          aln->hits[i].cigar = NULL;
                      }
                      j++;
                  }
                  else { // free
                      free(aln->hits[i].cigar);
                      aln->hits[i].cigar=NULL;
                  }
                  k++;
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      aln->hits[j] = aln->hits[i];
                      aln->hits[i].cigar = NULL;
                  }
                  j++;
              }
          }
          fmap_map_all_aln_realloc(aln, j);
      }
  }
  else if(FMAP_MAP_UTIL_ALN_MODE_ALL_BEST == opt->aln_output_mode) {
      // do nothing
  }
  else {
      fmap_error("bug encountered", Exit, OutOfRange);
  }
}

static fmap_map_all_aln_t *
fmap_map_all_aln_merge(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_bwt_t *bwt[2], fmap_sa_t *sa[2],
                       fmap_map1_aln_t *aln_map1, fmap_map2_sam_t *aln_map2, fmap_map3_aln_t *aln_map3,
                       int32_t stage, fmap_map_all_opt_t *opt)
{
  int32_t i, j, k, aln_i;
  uint32_t pacpos;
  uint32_t seqid, pos;
  fmap_map_all_aln_t *aln = NULL;
  int32_t seq_len = 0;

  aln = fmap_map_all_aln_init();

  // sequence length
  seq_len = fmap_seq_get_bases(seq)->l;
  if(FMAP_SEQ_TYPE_SFF == seq->type) {
      seq_len -= seq->data.sff->gheader->key_length; // soft clip the key sequence
  }

  // pre-allocate
  aln->n = 0;
  for(i=0;i<aln_map1->n;i++) {
      aln->n += aln_map1->hits[i].l - aln_map1->hits[i].k + 1;
  }
  aln->n += aln_map2->num_entries;
  aln->n += aln_map3->n;
  aln->hits = fmap_malloc(sizeof(fmap_map_all_hit_t) * aln->n, "aln->hits");

  aln_i = 0;

  // copy over map1
  for(i=0;i<aln_map1->n;i++) {
      int32_t aln_ref_l = 0;
      fmap_map1_hit_t *h = &aln_map1->hits[i];

      // get the reference length
      for(j=0;j<h->n_cigar;j++) {
          switch((h->cigar[j] & 0xf)) {
            case BAM_CMATCH:
            case BAM_CDEL:
              aln_ref_l += (h->cigar[j] >> 4); break;
            default:
              break;
          }
      }

      // go through all hits
      for(j=aln_map1->hits[i].k;j<=aln_map1->hits[i].l;j++) {
          pacpos = bwt[1]->seq_len - fmap_sa_pac_pos(sa[1], bwt[1], j) - aln_ref_l + 1;
          if(0 < fmap_refseq_pac2real(refseq, pacpos, seq_len, &seqid, &pos)) {
              // common info
              aln->hits[aln_i].algo_id = FMAP_MAP_ALL_ALGO_MAP1;
              aln->hits[aln_i].algo_stage = stage;
              aln->hits[aln_i].strand = h->strand;
              aln->hits[aln_i].seqid = seqid;
              aln->hits[aln_i].pos = pos - 1; // zero-based
              aln->hits[aln_i].score = h->score;
              // copy cigar - cannot shallow copy here
              aln->hits[aln_i].cigar = fmap_malloc(sizeof(uint32_t) * h->n_cigar, "aln->hits[aln_i]");
              for(k=0;k<h->n_cigar;k++) {
                  aln->hits[aln_i].cigar[k] = h->cigar[k];
              }
              aln->hits[aln_i].n_cigar = h->n_cigar;

              // map1 specific info
              aln->hits[aln_i].n_mm = h->n_mm;
              aln->hits[aln_i].n_gapo = h->n_gapo;
              aln->hits[aln_i].n_gape = h->n_gape;

              aln_i++;
          }
      }
  }

  // copy over map2
  for(i=0;i<aln_map2->num_entries;i++) {
      fmap_map2_sam_entry_t *e = &aln_map2->entries[i];

      // common info
      aln->hits[aln_i].algo_id = FMAP_MAP_ALL_ALGO_MAP2;
      aln->hits[aln_i].algo_stage = stage;
      aln->hits[aln_i].strand = e->strand;
      aln->hits[aln_i].seqid = e->seqid;
      aln->hits[aln_i].pos = e->pos;
      aln->hits[aln_i].score = e->AS;
      aln->hits[aln_i].cigar = e->cigar;
      aln->hits[aln_i].n_cigar = e->n_cigar;
      e->cigar = NULL;

      // map2 specific info
      aln->hits[aln_i].XF = e->XF;
      aln->hits[aln_i].XI = e->XI;
      aln->hits[aln_i].n_seeds = e->XE;
      aln->hits[aln_i].score_subo= e->XS;

      aln_i++;
  }

  // copy over map3
  for(i=0;i<aln_map3->n;i++) {
      fmap_map3_hit_t *h = &aln_map3->hits[i];

      // common info
      aln->hits[aln_i].algo_id = FMAP_MAP_ALL_ALGO_MAP3;
      aln->hits[aln_i].algo_stage = stage;
      aln->hits[aln_i].strand = h->strand;
      aln->hits[aln_i].seqid = h->seqid;
      aln->hits[aln_i].pos = h->pos ;
      aln->hits[aln_i].score = h->score;
      aln->hits[aln_i].cigar = h->cigar;
      aln->hits[aln_i].n_cigar = h->n_cigar;
      h->cigar = NULL;

      // map3 specific info
      aln->hits[aln_i].n_seeds = h->n_seeds;
      aln->hits[aln_i].score_subo= h->score_subo;

      aln_i++;
  }

  // reallocate/resize
  aln->hits = fmap_realloc(aln->hits, sizeof(fmap_map_all_hit_t) * aln_i, "aln->hits");
  aln->n = aln_i;

  // remove duplicates
  if(0 == opt->aln_output_mode_ind) {
      // sort
      fmap_sort_introsort(fmap_map_all_hit_t, aln->n, aln->hits);

      // remove duplicates within a window
      for(i=j=0;i<aln->n;) {
          int32_t end, best_score_i;

          // get the change
          end = best_score_i = i;
          while(end+1< aln->n) {
              if(aln->hits[end].seqid == aln->hits[end+1].seqid
                 && fabs(aln->hits[end].pos - aln->hits[end+1].pos) <= opt->dup_window) {
                  // track the best scoring
                  if(aln->hits[best_score_i].score < aln->hits[end].score) {
                      best_score_i = end+1;
                  }
                  end++;
              }
              else {
                  break;
              }
          }
          // TODO: randomize the best scoring

          // copy over the best
          if(j != best_score_i) {
              // free cigar
              free(aln->hits[j].cigar);
              // shallow copy
              aln->hits[j] = aln->hits[best_score_i];
              // destroy
              aln->hits[best_score_i].cigar = NULL;
          }

          // next
          i = end+1;
          j++;
      }

      // free cigars
      for(i=j;i<aln->n;i++) {
          free(aln->hits[i].cigar);
          aln->hits[i].cigar = NULL;
      }
      aln->n = j;
      aln->hits = fmap_realloc(aln->hits, sizeof(fmap_map_all_hit_t) * aln->n, "aln->hits");
  }

  // mapping quality
  fmap_map_all_aln_mapq(aln, opt);

  // choose alignment(s)
  if(0 == opt->aln_output_mode_ind) {
      // consider all algos together
      fmap_map_all_aln_filter(aln, opt, FMAP_MAP_ALL_ALGO_NONE);
  }
  else {
      // consider all algos independently
      fmap_map_all_aln_filter(aln, opt, FMAP_MAP_ALL_ALGO_MAP1);
      fmap_map_all_aln_filter(aln, opt, FMAP_MAP_ALL_ALGO_MAP2);
      fmap_map_all_aln_filter(aln, opt, FMAP_MAP_ALL_ALGO_MAP3);
  }

  return aln;
}

static inline void
fmap_map_all_print_sam(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_map_all_hit_t *hit)
{
  switch(hit->algo_id) {
    case FMAP_MAP_ALL_ALGO_MAP1:
      fmap_sam_print_mapped(fmap_file_stdout, seq, refseq,
                            hit->strand, hit->seqid, hit->pos,
                            hit->mapq, hit->cigar, hit->n_cigar,
                            "\tAS:i:%d\tNM:i:%d\tXM:i:%d\tXO:i:%d\tXG:i:%d\tXA:Z:%s-%d",
                            hit->score, 
                            (hit->n_mm + hit->n_gapo + hit->n_gape),
                            hit->n_mm, hit->n_gapo, hit->n_gape, algos[hit->algo_id], hit->algo_stage);
      break;
    case FMAP_MAP_ALL_ALGO_MAP2:
      if(0 < hit->XI) {
          fmap_sam_print_mapped(fmap_file_stdout, seq, refseq,
                                hit->strand, hit->seqid, hit->pos, hit->mapq,
                                hit->cigar, hit->n_cigar,
                                "\tAS:i:%d\tXS:i:%d\tXF:i:%d\tXE:i:%d\tXI:i:%d\tXA:Z:%s-%d",
                                hit->score, hit->score_subo, hit->XF, hit->n_seeds, hit->XI, 
                                algos[hit->algo_id], hit->algo_stage);
      }
      else {
          fmap_sam_print_mapped(fmap_file_stdout, seq, refseq,
                                hit->strand, hit->seqid, hit->pos, hit->mapq,
                                hit->cigar, hit->n_cigar,
                                "\tAS:i:%d\tXS:i:%d\tXF:i:%d\tXE:i:%d\tXA:Z:%s-%d",
                                hit->score, hit->score_subo, hit->XF, hit->n_seeds, 
                                algos[hit->algo_id], hit->algo_stage);
      }
      break;
    case FMAP_MAP_ALL_ALGO_MAP3:
      fmap_sam_print_mapped(fmap_file_stdout, seq, refseq,
                            hit->strand, hit->seqid, hit->pos,
                            hit->mapq, hit->cigar, hit->n_cigar,
                            "\tAS:i:%d\tXS:i:%d\tXE:i:%d\tXA:Z:%s-%d",
                            hit->score, hit->score_subo, hit->n_seeds, 
                            algos[hit->algo_id], hit->algo_stage);
      break;
    default:
      fmap_error("hit->algo_id", Exit, OutOfRange);
      break;
  }
}

static void
fmap_map_all_core_worker(fmap_seq_t **seq_buffer, fmap_map_all_aln_t **alns, int32_t seq_buffer_length, 
                         fmap_refseq_t *refseq, fmap_bwt_t *bwt[2], fmap_sa_t *sa[2],
                         int32_t tid, fmap_map_all_opt_t *opt)
{
  int32_t low = 0, high, i, j;
  // map1 
  fmap_bwt_match_width_t *width_map1[2][2], *seed_width_map1[2][2];
  int32_t width_length_map1[2] = {0, 0};
  fmap_map1_aux_stack_t *stack_map1=NULL;
  fmap_map1_aln_t *aln_map1;
  // map2
  fmap_map2_global_mempool_t *pool_map2 = NULL;
  fmap_map2_sam_t *aln_map2;
  // map3
  fmap_map3_aln_t *aln_map3;
  uint8_t *flow_map3[2][2];

  // map1
  for(i=0;i<2;i++) {
      if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP1) {
          width_map1[i][0] = width_map1[i][1] = NULL;
          width_length_map1[i] = 0;
          seed_width_map1[i][0] = fmap_calloc(opt->opt_map1[i]->seed_length, sizeof(fmap_bwt_match_width_t), "seed_width[0]");
          seed_width_map1[i][1] = fmap_calloc(opt->opt_map1[i]->seed_length, sizeof(fmap_bwt_match_width_t), "seed_width[1]");
          if(NULL == stack_map1) {
              stack_map1 = fmap_map1_aux_stack_init();
          }
      }
      // map2
      if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP2) {
          if(NULL == pool_map2) {
              pool_map2 = fmap_map2_global_mempool_init();
          }
      }
      // map3
      if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP3) {
          // set up the flow order
          if(0 < seq_buffer_length) {
              flow_map3[i][0] = fmap_malloc(sizeof(uint8_t) * 4, "flow[i][0]");
              flow_map3[i][1] = fmap_malloc(sizeof(uint8_t) * 4, "flow[i][1]");
              for(j=0;j<4;j++) {
                  if(FMAP_SEQ_TYPE_SFF == seq_buffer[0]->type) {
                      flow_map3[i][0][j] = fmap_nt_char_to_int[(int)seq_buffer[0]->data.sff->gheader->flow->s[j]];
                      flow_map3[i][1][j] = 3 - fmap_nt_char_to_int[(int)seq_buffer[0]->data.sff->gheader->flow->s[j]];
                  }
                  else {
                      flow_map3[i][0][j] = fmap_nt_char_to_int[(int)opt->opt_map3[i]->flow[j]];
                      flow_map3[i][1][j] = 3 - fmap_nt_char_to_int[(int)opt->opt_map3[i]->flow[j]];
                  }
              }
          }
      }
  }

  while(low < seq_buffer_length) {
#ifdef HAVE_LIBPTHREAD
      if(1 < opt->num_threads) {
          pthread_mutex_lock(&fmap_map_all_read_lock);

          // update bounds
          low = fmap_map_all_read_lock_low;
          fmap_map_all_read_lock_low += FMAP_MAP_ALL_THREAD_BLOCK_SIZE;
          high = low + FMAP_MAP_ALL_THREAD_BLOCK_SIZE;
          if(seq_buffer_length < high) {
              high = seq_buffer_length; 
          }

          pthread_mutex_unlock(&fmap_map_all_read_lock);
      }
      else {
          high = seq_buffer_length; // process all
      }
#else 
      high = seq_buffer_length; // process all
#endif
      while(low<high) {
          fmap_seq_t *seq[2]={NULL, NULL}, *orig_seq=NULL, *seq_char=NULL;
          orig_seq = seq_buffer[low];
          fmap_string_t *bases[2]={NULL, NULL};
          fmap_map1_opt_t opt_local_map1[2];

          // map1
          for(i=0;i<2;i++) {
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP1) {
                  opt_local_map1[i] = (*opt->opt_map1[i]); // copy over values
              }
              // map2
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP2) {
                  // - none
              }
              // map3
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP3) {
                  // - none
              }
          }

          // clone the sequence 
          seq[0] = fmap_seq_clone(orig_seq);
          seq[1] = fmap_seq_clone(orig_seq);
          seq_char = fmap_seq_clone(orig_seq);

          // Adjust for SFF
          fmap_seq_remove_key_sequence(seq[0]);
          fmap_seq_remove_key_sequence(seq[1]);
          fmap_seq_remove_key_sequence(seq_char);

          // reverse compliment
          fmap_seq_reverse_compliment(seq[1]);

          // convert to integers
          fmap_seq_to_int(seq[0]);
          fmap_seq_to_int(seq[1]);

          // get bases
          bases[0] = fmap_seq_get_bases(seq[0]);
          bases[1] = fmap_seq_get_bases(seq[1]);

          // map1
          for(i=0;i<2;i++) {
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP1) {
                  opt_local_map1[i].max_mm = (opt->opt_map1[i]->max_mm < 0) ? (int)(0.99 + opt->opt_map1[i]->max_mm_frac * bases[0]->l) : opt->opt_map1[i]->max_mm;
                  opt_local_map1[i].max_gape = (opt->opt_map1[i]->max_gape < 0) ? (int)(0.99 + opt->opt_map1[i]->max_gape_frac * bases[0]->l) : opt->opt_map1[i]->max_gape;
                  opt_local_map1[i].max_gapo = (opt->opt_map1[i]->max_gapo < 0) ? (int)(0.99 + opt->opt_map1[i]->max_gapo_frac * bases[0]->l) : opt->opt_map1[i]->max_gapo;
                  if(width_length_map1[i] < bases[0]->l) {
                      free(width_map1[i][0]); free(width_map1[i][1]);
                      width_length_map1[i] = bases[0]->l;
                      width_map1[i][0] = fmap_calloc(width_length_map1[i], sizeof(fmap_bwt_match_width_t), "width[0]");
                      width_map1[i][1] = fmap_calloc(width_length_map1[i], sizeof(fmap_bwt_match_width_t), "width[1]");
                  }
                  fmap_bwt_match_cal_width(bwt[0], bases[0]->l, bases[0]->s, width_map1[i][0]);
                  fmap_bwt_match_cal_width(bwt[0], bases[1]->l, bases[1]->s, width_map1[i][1]);

                  if(bases[0]->l < opt->opt_map1[i]->seed_length) {
                      opt_local_map1[i].seed_length = -1;
                  }
                  else {
                      fmap_bwt_match_cal_width(bwt[0], opt->opt_map1[i]->seed_length, bases[0]->s, seed_width_map1[i][0]);
                      fmap_bwt_match_cal_width(bwt[0], opt->opt_map1[i]->seed_length, bases[1]->s, seed_width_map1[i][1]);
                  }
              }
              // map2
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP2) {
                  // - none
              }
              // map3
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP3) {
                  // - none
              }
          }

          // run the core algorithms
          for(i=0;i<2;i++) {
              // fmap_map1_aux_core
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP1) {
                  aln_map1 = fmap_map1_aux_core(seq, bwt[1], width_map1[i], 
                                                (0 < opt_local_map1[i].seed_length) ? seed_width_map1[i] : NULL, 
                                                &opt_local_map1[i], stack_map1);
                  // adjust map1 scoring, since it does not consider opt->score_match
                  for(j=0;j<aln_map1->n;j++) {
                      fmap_map1_hit_t *h = &aln_map1->hits[j];
                      int32_t k, num_match;

                      num_match = 0 - h->n_mm; // # of mismatches
                      for(k=0;k<h->n_cigar;k++) {
                          switch(h->cigar[k] & 0xf) {
                            case BAM_CMATCH:
                              num_match += h->cigar[k] >> 4;
                              break;
                            default:
                              break;
                          }
                      }

                      // update the score
                      h->score = num_match * opt->score_match;
                      h->score -= h->n_mm * opt->pen_mm;
                      h->score -= h->n_gapo * opt->pen_gapo;
                      h->score -= h->n_gape * opt->pen_gape;
                  }
              }
              else {
                  // empty
                  aln_map1 = fmap_map1_aln_init();
              }
              // fmap_map2_aux_core
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP2) {
                  aln_map2 = fmap_map2_aux_core(opt->opt_map2[i], seq_char, refseq, bwt, sa, pool_map2);
              }
              else {
                  // empty
                  aln_map2 = fmap_map2_sam_init(0);
              }
              // fmap_map3_aux_core
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP3) {
                  aln_map3 = fmap_map3_aux_core(seq, flow_map3[i], refseq, bwt[1], sa[1], opt->opt_map3[i]);
              }
              else {
                  // empty
                  aln_map3 = fmap_map3_aln_init();
              }

              // consolidate mappings
              alns[low] = fmap_map_all_aln_merge(orig_seq, refseq, bwt, sa,
                                                 aln_map1, aln_map2, aln_map3, i, opt);

              // destroy
              // map1
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP1) {
                  fmap_map1_aln_destroy(aln_map1);
                  aln_map1 = NULL;
              }
              // map2
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP2) {
                  fmap_map2_sam_destroy(aln_map2);
                  aln_map2 = NULL;
              }
              // map3
              if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP3) {
                  fmap_map3_aln_destroy(aln_map3);
                  aln_map3 = NULL;
              }

              // check if we found any mappings
              if(0 < alns[low]->n) {
                  // yes we did
                  break;
              }
          }

          // destroy
          fmap_seq_destroy(seq[0]);
          fmap_seq_destroy(seq[1]);
          fmap_seq_destroy(seq_char);

          // next
          low++;
      }
  }

  // map1
  for(i=0;i<2;i++) {
      if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP1) {
          if(NULL != stack_map1) {
              fmap_map1_aux_stack_destroy(stack_map1);
              stack_map1 = NULL;
          }
          free(seed_width_map1[i][0]);
          free(seed_width_map1[i][1]);
          free(width_map1[i][0]);
          free(width_map1[i][1]);
      }
      // map2
      if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP2) {
          if(NULL != pool_map2) {
              fmap_map2_global_mempool_destroy(pool_map2);
              pool_map2 = NULL;
          }
      }
      // map3
      if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP3) {
          free(flow_map3[i][0]);
          free(flow_map3[i][1]);
      }
  }
}

static void *
fmap_map_all_core_thread_worker(void *arg)
{
  fmap_map_all_thread_data_t *thread_data = (fmap_map_all_thread_data_t*)arg;

  fmap_map_all_core_worker(thread_data->seq_buffer, thread_data->alns, thread_data->seq_buffer_length, 
                           thread_data->refseq, thread_data->bwt, thread_data->sa, 
                           thread_data->tid, thread_data->opt);

  return arg;
}

static void 
fmap_map_all_core(fmap_map_all_opt_t *opt)
{
  uint32_t i, j, n_reads_processed=0;
  int32_t seq_buffer_length;
  fmap_refseq_t *refseq=NULL;
  fmap_bwt_t *bwt[2]={NULL, NULL};
  fmap_sa_t *sa[2]={NULL, NULL};
  fmap_file_t *fp_reads=NULL;
  fmap_seq_io_t *seqio = NULL;
  fmap_seq_t **seq_buffer = NULL;
  fmap_map_all_aln_t **alns = NULL;
  fmap_shm_t *shm = NULL;
  int32_t reads_queue_size;
  
  if(NULL == opt->fn_reads) {
      fmap_progress_set_verbosity(0); 
  }

  if(0 == opt->shm_key) {
      fmap_progress_print("reading in reference data");
      refseq = fmap_refseq_read(opt->fn_fasta, 0);
      bwt[0] = fmap_bwt_read(opt->fn_fasta, 0);
      bwt[1] = fmap_bwt_read(opt->fn_fasta, 1);
      sa[0] = fmap_sa_read(opt->fn_fasta, 0);
      sa[1] = fmap_sa_read(opt->fn_fasta, 1);
      fmap_progress_print2("reference data read in");
  }
  else {
      fmap_progress_print("retrieving reference data from shared memory");
      shm = fmap_shm_init(opt->shm_key, 0, 0);
      if(NULL == (refseq = fmap_refseq_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_REFSEQ)))) {
          fmap_error("the packed reference sequence was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (bwt[0] = fmap_bwt_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_BWT)))) {
          fmap_error("the BWT string was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (bwt[1] = fmap_bwt_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_REV_BWT)))) {
          fmap_error("the reverse BWT string was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (sa[0] = fmap_sa_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_SA)))) {
          fmap_error("the SA was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (sa[1] = fmap_sa_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_REV_SA)))) {
          fmap_error("the reverse SA was not found in shared memory", Exit, SharedMemoryListing);
      }
      fmap_progress_print2("reference data retrieved from shared memory");

  }

  // TODO: this could be dangerous if algorithms change, needs to be refactored
  // adjust mapping algorithm specific options
  for(i=0;i<2;i++) {
      // map1
      if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP1) {
          // - none
      }
      // map2
      if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP2) {
          opt->opt_map2[i]->score_thr *= opt->score_match;
          opt->opt_map2[i]->length_coef *= opt->score_match;
      }
      // map3
      if(opt->algos[i] & FMAP_MAP_ALL_ALGO_MAP3) {
          opt->opt_map3[i]->score_thr *= opt->score_match;
          if(-1 == opt->opt_map3[i]->seed_length) {
              opt->opt_map3[i]->seed_length = fmap_map3_get_seed_length(refseq->len);
              fmap_progress_print("setting the seed length to %d for map3",
                                  opt->opt_map3[i]->seed_length);
          }
      }
  }
  // allocate the buffer
  if(-1 == opt->reads_queue_size) {
      reads_queue_size = 1;
  }
  else {
      reads_queue_size = opt->reads_queue_size;
  }
  seq_buffer = fmap_malloc(sizeof(fmap_seq_t*)*reads_queue_size, "seq_buffer");
  alns = fmap_malloc(sizeof(fmap_map_all_aln_t*)*reads_queue_size, "alns");

  if(NULL == opt->fn_reads) {
      fp_reads = fmap_file_fdopen(fileno(stdin), "rb", opt->input_compr);
      fmap_progress_set_verbosity(0); 
  }
  else {
      fp_reads = fmap_file_fopen(opt->fn_reads, "rb", opt->input_compr);
  }
  switch(opt->reads_format) {
    case FMAP_READS_FORMAT_FASTA:
    case FMAP_READS_FORMAT_FASTQ:
      seqio = fmap_seq_io_init(fp_reads, FMAP_SEQ_TYPE_FQ);
      for(i=0;i<reads_queue_size;i++) { // initialize the buffer
          seq_buffer[i] = fmap_seq_init(FMAP_SEQ_TYPE_FQ);
      }
      break;
    case FMAP_READS_FORMAT_SFF:
      seqio = fmap_seq_io_init(fp_reads, FMAP_SEQ_TYPE_SFF);
      for(i=0;i<reads_queue_size;i++) { // initialize the buffer
          seq_buffer[i] = fmap_seq_init(FMAP_SEQ_TYPE_SFF);
      }
      break;
    default:
      fmap_error("unrecognized input format", Exit, CommandLineArgument);
      break;
  }

  // Note: 'fmap_file_stdout' should not have been previously modified
  fmap_file_stdout = fmap_file_fdopen(fileno(stdout), "wb", opt->output_compr);

  // SAM header
  fmap_sam_print_header(fmap_file_stdout, refseq, seqio, opt->sam_rg, opt->argc, opt->argv);

  fmap_progress_print("processing reads");
  while(0 < (seq_buffer_length = fmap_seq_io_read_buffer(seqio, seq_buffer, reads_queue_size))) {

      // do alignment
#ifdef HAVE_LIBPTHREAD
      int32_t num_threads = opt->num_threads;
      if(seq_buffer_length < num_threads * FMAP_MAP_ALL_THREAD_BLOCK_SIZE) {
          num_threads = 1 + (seq_buffer_length / FMAP_MAP_ALL_THREAD_BLOCK_SIZE);
      }
      fmap_map_all_read_lock_low = 0; // ALWAYS set before running threads 
      if(1 == num_threads) {
          fmap_map_all_core_worker(seq_buffer, alns, seq_buffer_length, refseq, bwt, sa, 0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          fmap_map_all_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = fmap_calloc(num_threads, sizeof(pthread_t), "threads");
          thread_data = fmap_calloc(num_threads, sizeof(fmap_map_all_thread_data_t), "thread_data");

          for(i=0;i<num_threads;i++) {
              thread_data[i].seq_buffer = seq_buffer;
              thread_data[i].seq_buffer_length = seq_buffer_length;
              thread_data[i].alns = alns;
              thread_data[i].refseq = refseq;
              thread_data[i].bwt[0] = bwt[0];
              thread_data[i].bwt[1] = bwt[1];
              thread_data[i].sa[0] = sa[0];
              thread_data[i].sa[1] = sa[1];
              thread_data[i].tid = i;
              thread_data[i].opt = opt; 
              if(0 != pthread_create(&threads[i], &attr, fmap_map_all_core_thread_worker, &thread_data[i])) {
                  fmap_error("error creating threads", Exit, ThreadError);
              }
          }
          for(i=0;i<num_threads;i++) {
              if(0 != pthread_join(threads[i], NULL)) {
                  fmap_error("error joining threads", Exit, ThreadError);
              }
          }

          free(threads);
          free(thread_data);
      }
#else 
      fmap_map_all_core_worker(seq_buffer, alns, seq_buffer_length, refseq, bwt, sa, 0, opt);
#endif

      if(-1 != opt->reads_queue_size) {
          fmap_progress_print("writing alignments");
      }
      for(i=0;i<seq_buffer_length;i++) {
          if(0 < alns[i]->n) {
              for(j=0;j<alns[i]->n;j++) {
                  fmap_map_all_print_sam(seq_buffer[i], refseq, &alns[i]->hits[j]);
              }
          }
          else {
              fmap_sam_print_unmapped(fmap_file_stdout, seq_buffer[i]);
          }

          // free alignments
          fmap_map_all_aln_destroy(alns[i]);
          alns[i] = NULL;
      }

      if(-1 == opt->reads_queue_size) {
          fmap_file_fflush(fmap_file_stdout, 1);
      }

      n_reads_processed += seq_buffer_length;
      if(-1 != opt->reads_queue_size) {
          fmap_progress_print2("processed %d reads", n_reads_processed);
      }
  }
  if(-1 == opt->reads_queue_size) {
      fmap_progress_print2("processed %d reads", n_reads_processed);
  }

  // close the input/output
  fmap_file_fclose(fmap_file_stdout);
  fmap_file_fclose(fp_reads);

  // free memory
  for(i=0;i<reads_queue_size;i++) {
      fmap_seq_destroy(seq_buffer[i]);
  }
  free(seq_buffer);
  free(alns);
  fmap_refseq_destroy(refseq);
  fmap_bwt_destroy(bwt[0]);
  fmap_bwt_destroy(bwt[1]);
  fmap_sa_destroy(sa[0]);
  fmap_sa_destroy(sa[1]);
  fmap_seq_io_destroy(seqio);
  if(0 < opt->shm_key) {
      fmap_shm_destroy(shm, 0);
  }
}

int 
fmap_map_all_usage(fmap_map_all_opt_t *opt)
{
  char *reads_format = fmap_get_reads_file_format_string(opt->reads_format);

  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Usage: %s mapall [global options] (algorithm [options])+ (ALGORITHM [options])*", PACKAGE);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (required):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -f FILE     the FASTA reference file name [%s]\n", opt->fn_fasta);
  fmap_file_fprintf(fmap_file_stderr, "         -r FILE     the reads file name [%s]\n", (NULL == opt->fn_reads) ? "stdin" : opt->fn_reads);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (optional):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff) [%s]\n", reads_format);
  fmap_file_fprintf(fmap_file_stderr, "         -A INT      the match score [%d]\n", opt->score_match); 
  fmap_file_fprintf(fmap_file_stderr, "         -M INT      the mismatch penalty [%d]\n", opt->pen_mm); 
  fmap_file_fprintf(fmap_file_stderr, "         -O INT      the indel start penalty [%d]\n", opt->pen_gapo); 
  fmap_file_fprintf(fmap_file_stderr, "         -E INT      the indel extend penalty [%d]\n", opt->pen_gape); 
  fmap_file_fprintf(fmap_file_stderr, "         -X INT      the flow score penalty [%d]\n", opt->fscore);
  fmap_file_fprintf(fmap_file_stderr, "         -w INT      the extra bases to add before and after the target during Smith-Waterman [%d]\n", opt->bw);
  fmap_file_fprintf(fmap_file_stderr, "         -g          align the full read (global alignment) [%s]\n", (0 == opt->aln_global) ? "false" : "true");
  fmap_file_fprintf(fmap_file_stderr, "         -q INT      the queue size for the reads (-1 disables) [%d]\n", opt->reads_queue_size);
  fmap_file_fprintf(fmap_file_stderr, "         -n INT      the number of threads [%d]\n", opt->num_threads);
  fmap_file_fprintf(fmap_file_stderr, "         -a INT      output filter [%d]\n", opt->aln_output_mode);
  fmap_file_fprintf(fmap_file_stderr, "                             0 - unique best hits\n");
  fmap_file_fprintf(fmap_file_stderr, "                             1 - random best hit\n");
  fmap_file_fprintf(fmap_file_stderr, "                             2 - all best hits\n");
  fmap_file_fprintf(fmap_file_stderr, "                             3 - all alignments\n");
  fmap_file_fprintf(fmap_file_stderr, "         -R STRING   the RG line in the SAM header [%s]\n", opt->sam_rg);
  fmap_file_fprintf(fmap_file_stderr, "         -W INT      remove duplicate alignments from different algorithms within this bp window [%d]\n",
                    opt->dup_window);
  fmap_file_fprintf(fmap_file_stderr, "         -I          apply the output filter for each algorithm separately [%s]\n",
                    (1 == opt->aln_output_mode_ind) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -j          the input is bz2 compressed (bzip2) [%s]\n",
                    (FMAP_FILE_BZ2_COMPRESSION == opt->input_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -z          the input is gz compressed (gzip) [%s]\n",
                    (FMAP_FILE_GZ_COMPRESSION == opt->input_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -J          the output is bz2 compressed (bzip2) [%s]\n",
                    (FMAP_FILE_BZ2_COMPRESSION == opt->output_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -Z          the output is gz compressed (gzip) [%s]\n",
                    (FMAP_FILE_GZ_COMPRESSION == opt->output_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -s INT      use shared memory with the following key [%d]\n", opt->shm_key);
  fmap_file_fprintf(fmap_file_stderr, "         -v          print verbose progress information\n");
  fmap_file_fprintf(fmap_file_stderr, "         -h          print this message\n");
  fmap_file_fprintf(fmap_file_stderr, "\n");
  // TODO: more description of how to use map{1,2,3} and 2-stage algo

  free(reads_format);

  return 1;
}

fmap_map_all_opt_t *
fmap_map_all_opt_init()
{
  int32_t i;
  fmap_map_all_opt_t *opt = NULL;

  opt = fmap_calloc(1, sizeof(fmap_map_all_opt_t), "opt");

  // program defaults
  opt->argv = NULL;
  opt->argc = -1;
  opt->fn_fasta = opt->fn_reads = NULL;
  opt->reads_format = FMAP_READS_FORMAT_UNKNOWN;
  opt->score_match = FMAP_MAP_UTIL_SCORE_MATCH;
  opt->pen_mm = FMAP_MAP_UTIL_PEN_MM; 
  opt->pen_gapo = FMAP_MAP_UTIL_PEN_GAPO;
  opt->pen_gape = FMAP_MAP_UTIL_PEN_GAPE;
  opt->fscore = FMAP_MAP_UTIL_FSCORE;
  opt->bw = 10; 
  opt->aln_global = 0;
  opt->reads_queue_size = 65536; 
  opt->num_threads = 1;
  opt->aln_output_mode = FMAP_MAP_UTIL_ALN_MODE_RAND_BEST;
  opt->dup_window = 128;
  opt->aln_output_mode_ind = 0;
  opt->sam_rg = NULL;
  opt->input_compr = FMAP_FILE_NO_COMPRESSION;
  opt->output_compr = FMAP_FILE_NO_COMPRESSION;
  opt->shm_key = 0;

  // these must be checked for agreement above later
  for(i=0;i<2;i++) {
      opt->algos[i] = 0;
      opt->opt_map1[i] = fmap_map1_opt_init();
      opt->opt_map2[i] = fmap_map2_opt_init();
      opt->opt_map3[i] = fmap_map3_opt_init();
  }

  return opt;
}

void
fmap_map_all_opt_destroy(fmap_map_all_opt_t *opt)
{
  int32_t i;

  // free common options
  free(opt->fn_fasta);
  free(opt->fn_reads);
  free(opt->sam_rg);

  for(i=0;i<2;i++) {
      // since we shallow copied, do not free
      opt->opt_map1[i]->fn_fasta = opt->opt_map2[i]->fn_fasta = opt->opt_map3[i]->fn_fasta = NULL;
      opt->opt_map1[i]->fn_reads = opt->opt_map2[i]->fn_reads = opt->opt_map3[i]->fn_reads = NULL;
      opt->opt_map1[i]->sam_rg = opt->opt_map2[i]->sam_rg = opt->opt_map3[i]->sam_rg = NULL;

      // destroy other opts
      fmap_map1_opt_destroy(opt->opt_map1[i]);
      fmap_map2_opt_destroy(opt->opt_map2[i]);
      fmap_map3_opt_destroy(opt->opt_map3[i]);
  }

  // free the current opt
  free(opt);
}

// for map1/map2/map3
#define __fmap_map_all_opts_copy1(opt_map_all, opt_map_other) do { \
    (opt_map_other)->fn_fasta = (opt_map_all)->fn_fasta; \
    (opt_map_other)->fn_reads = (opt_map_all)->fn_reads; \
    (opt_map_other)->reads_format = (opt_map_all)->reads_format; \
    (opt_map_other)->pen_mm = (opt_map_all)->pen_mm; \
    (opt_map_other)->pen_gapo = (opt_map_all)->pen_gapo; \
    (opt_map_other)->pen_gape = (opt_map_all)->pen_gape; \
    (opt_map_other)->reads_queue_size = (opt_map_all)->reads_queue_size; \
    (opt_map_other)->num_threads = (opt_map_all)->num_threads; \
    (opt_map_other)->aln_output_mode = FMAP_MAP_UTIL_ALN_MODE_ALL; \
    (opt_map_other)->sam_rg = (opt_map_all)->sam_rg; \
    (opt_map_other)->input_compr = (opt_map_all)->input_compr; \
    (opt_map_other)->output_compr = (opt_map_all)->output_compr; \
    (opt_map_other)->shm_key = (opt_map_all)->shm_key; \
} while(0)

// for map2 and map3
#define __fmap_map_all_opts_copy2(opt_map_all, opt_map_other) do { \
    __fmap_map_all_opts_copy1(opt_map_all, opt_map_other); \
    (opt_map_other)->score_match = (opt_map_all)->score_match; \
    (opt_map_other)->fscore = (opt_map_all)->fscore; \
    (opt_map_other)->bw = (opt_map_all)->bw; \
    (opt_map_other)->aln_global = (opt_map_all)->aln_global; \
} while(0)

static int32_t
fmap_map_all_opt_parse_common(int argc, char *argv[], fmap_map_all_opt_t *opt)
{
  int c;

  while((c = getopt(argc, argv, "f:r:F:A:M:O:E:X:w:gq:n:a:R:W:IjzJZs:vh")) >= 0) {
      switch(c) {
        case 'f':
          opt->fn_fasta = fmap_strdup(optarg); break;
        case 'r':
          opt->fn_reads = fmap_strdup(optarg); 
          fmap_get_reads_file_format_from_fn_int(opt->fn_reads, &opt->reads_format, &opt->input_compr);
          break;
        case 'F':
          opt->reads_format = fmap_get_reads_file_format_int(optarg); break;
        case 'A':
          opt->score_match = atoi(optarg); break;
        case 'M':
          opt->pen_mm = atoi(optarg); break;
        case 'O':
          opt->pen_gapo = atoi(optarg); break;
        case 'E':
          opt->pen_gape = atoi(optarg); break;
        case 'X':
          opt->fscore = atoi(optarg); break;
        case 'w':
          opt->bw = atoi(optarg); break;
        case 'g':
          opt->aln_global = 1; break;
        case 'q': 
          opt->reads_queue_size = atoi(optarg); break;
        case 'n':
          opt->num_threads = atoi(optarg); break;
        case 'a':
          opt->aln_output_mode = atoi(optarg); break;
        case 'R':
          opt->sam_rg = fmap_strdup(optarg); break;
        case 'W':
          opt->dup_window = atoi(optarg); break;
        case 'I':
          opt->aln_output_mode_ind = 1; break;
        case 'j':
          opt->input_compr = FMAP_FILE_BZ2_COMPRESSION; 
          fmap_get_reads_file_format_from_fn_int(opt->fn_reads, &opt->reads_format, &opt->input_compr);
          break;
        case 'z':
          opt->input_compr = FMAP_FILE_GZ_COMPRESSION; 
          fmap_get_reads_file_format_from_fn_int(opt->fn_reads, &opt->reads_format, &opt->input_compr);
          break;
        case 'J':
          opt->output_compr = FMAP_FILE_BZ2_COMPRESSION; break;
        case 'Z':
          opt->output_compr = FMAP_FILE_GZ_COMPRESSION; break;
        case 's':
          opt->shm_key = atoi(optarg); break;
        case 'v':
          fmap_progress_set_verbosity(1); break;
        case 'h':
        default:
          return 0;
      }
  }
  return 1;
}

int32_t
fmap_map_all_opt_parse(int argc, char *argv[], fmap_map_all_opt_t *opt)
{
  int32_t i, j, start, opt_type, opt_type_next, opt_stage, opt_stage_next;

  opt->argc = argc; opt->argv = argv;

  // parse common options as well as map1/map2/map3 commands
  start = 0;
  i = 1;
  opt_type = opt_type_next = FMAP_MAP_ALL_ALGO_NONE;
  opt_stage = opt_stage_next = 0;
  while(i<argc) {
      if(0 == strcmp("map1", argv[i])) {
          opt_type_next = FMAP_MAP_ALL_ALGO_MAP1;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
      }
      else if(0 == strcmp("map2", argv[i])) {
          opt_type_next = FMAP_MAP_ALL_ALGO_MAP2;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
      }
      else if(0 == strcmp("map3", argv[i])) {
          opt_type_next = FMAP_MAP_ALL_ALGO_MAP3;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
      }
      else if(0 == strcmp("MAP1", argv[i])) {
          opt_type_next = FMAP_MAP_ALL_ALGO_MAP1;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
      }
      else if(0 == strcmp("MAP2", argv[i])) {
          opt_type_next = FMAP_MAP_ALL_ALGO_MAP2;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
      }
      else if(0 == strcmp("MAP3", argv[i])) {
          opt_type_next = FMAP_MAP_ALL_ALGO_MAP3;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
      }

      /*
         fprintf(stderr, "i=%d start=%d argc=%d opt_type=%d opt_type_next=%d argv[%d]=%s\n",
         i, start, argc, opt_type, opt_type_next, i, argv[i]);
         */

      if(opt_type != opt_type_next
         || i == argc-1) {
          if(i == argc-1) {
              i++;
          }
          optind=1; // needed for getopt
          switch(opt_type) {
            case FMAP_MAP_ALL_ALGO_NONE:
              // parse common options
              if(0 == fmap_map_all_opt_parse_common(i-start, argv+start, opt)) {
                  return 0;
              }
              // copy over common values into the other opts
              for(j=0;j<2;j++) {
                  __fmap_map_all_opts_copy1(opt, opt->opt_map1[j]);
                  __fmap_map_all_opts_copy2(opt, opt->opt_map2[j]);
                  __fmap_map_all_opts_copy2(opt, opt->opt_map3[j]);
              }
              break;
            case FMAP_MAP_ALL_ALGO_MAP1:
              // parse map1 options
              if(0 < i - start) {
                  if(0 == fmap_map1_opt_parse(i-start, argv+start, opt->opt_map1[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            case FMAP_MAP_ALL_ALGO_MAP2:
              // parse map2 options
              if(0 < i - start) {
                  if(0 == fmap_map2_opt_parse(i-start, argv+start, opt->opt_map2[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            case FMAP_MAP_ALL_ALGO_MAP3:
              // parse map3 options
              if(0 < i - start) {
                  if(0 == fmap_map3_opt_parse(i-start, argv+start, opt->opt_map3[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            default:
              fmap_error("bug encountered", Exit, OutOfRange);
          }
          opt_type = opt_type_next;
          opt_stage = opt_stage_next;
          start = i;
      }
      i++;
      if(argc < i) {
          i = argc;
      }
  }
  optind = i;

  return 1;
}

static int32_t 
fmap_map_all_file_check_with_null(char *fn1, char *fn2)
{
    if(NULL == fn1 && NULL == fn2) { 
        return 0;
    }
    else if((NULL == fn1 && NULL != fn2) 
       || (NULL != fn1 && NULL == fn2)) {
        return 1;
    }
    else if(0 != strcmp(fn1, fn2)) {
        return 1;
    }
    return 0;
}

// for map1
#define __fmap_map_all_opts_check_common1(opt_map_all, opt_map_other) do { \
    if(0 != fmap_map_all_file_check_with_null((opt_map_other)->fn_fasta, (opt_map_all)->fn_fasta)) { \
        fmap_error("option -f was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if(0 != fmap_map_all_file_check_with_null((opt_map_other)->fn_reads, (opt_map_all)->fn_reads)) { \
        fmap_error("option -r was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->reads_format != (opt_map_all)->reads_format) { \
        fmap_error("option -F was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->pen_mm != (opt_map_all)->pen_mm) { \
        fmap_error("option -M was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->pen_gapo != (opt_map_all)->pen_gapo) { \
        fmap_error("option -O was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->pen_gape != (opt_map_all)->pen_gape) { \
        fmap_error("option -E was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->reads_queue_size != (opt_map_all)->reads_queue_size) { \
        fmap_error("option -q was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->num_threads != (opt_map_all)->num_threads) { \
        fmap_error("option -n was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if(0 != fmap_map_all_file_check_with_null((opt_map_other)->sam_rg, (opt_map_all)->sam_rg)) { \
        fmap_error("option -R was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->input_compr != (opt_map_all)->input_compr) { \
        fmap_error("option -j or -z was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->output_compr != (opt_map_all)->output_compr) { \
        fmap_error("option -J or -Z was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->shm_key != (opt_map_all)->shm_key) { \
        fmap_error("option -s was specified outside of the common options", Exit, CommandLineArgument); \
    } \
} while(0)

// for map2/map3
#define __fmap_map_all_opts_check_common2(opt_map_all, opt_map_other) do { \
    __fmap_map_all_opts_check_common1(opt_map_all, opt_map_other); \
    if((opt_map_other)->score_match != (opt_map_all)->score_match) { \
        fmap_error("option -A was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->fscore != (opt_map_all)->fscore) { \
        fmap_error("option -X was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->bw != (opt_map_all)->bw) { \
        fmap_error("option -w was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->aln_global != (opt_map_all)->aln_global) { \
        fmap_error("option -g was specified outside of the common options", Exit, CommandLineArgument); \
    } \
} while(0)

void
fmap_map_all_opt_check(fmap_map_all_opt_t *opt)
{
  int32_t i;
  // common options
  if(NULL == opt->fn_fasta && 0 == opt->shm_key) {
      fmap_error("option -f or option -s must be specified", Exit, CommandLineArgument);
  }
  else if(NULL != opt->fn_fasta && 0 < opt->shm_key) {
      fmap_error("option -f and option -s may not be specified together", Exit, CommandLineArgument);
  }
  if(NULL == opt->fn_reads && FMAP_READS_FORMAT_UNKNOWN == opt->reads_format) {
      fmap_error("option -F or option -r must be specified", Exit, CommandLineArgument);
  }
  if(FMAP_READS_FORMAT_UNKNOWN == opt->reads_format) {
      fmap_error("the reads format (-r) was unrecognized", Exit, CommandLineArgument);
  }
  fmap_error_cmd_check_int(opt->score_match, 0, INT32_MAX, "-A");
  fmap_error_cmd_check_int(opt->pen_mm, 0, INT32_MAX, "-M");
  fmap_error_cmd_check_int(opt->pen_gapo, 0, INT32_MAX, "-O");
  fmap_error_cmd_check_int(opt->pen_gape, 0, INT32_MAX, "-E");
  fmap_error_cmd_check_int(opt->fscore, 0, INT32_MAX, "-X");
  fmap_error_cmd_check_int(opt->bw, 1, INT32_MAX, "-w");
  if(-1 != opt->reads_queue_size) fmap_error_cmd_check_int(opt->reads_queue_size, 1, INT32_MAX, "-q");
  fmap_error_cmd_check_int(opt->num_threads, 1, INT32_MAX, "-n");
  fmap_error_cmd_check_int(opt->aln_output_mode, 0, 3, "-a");
  fmap_error_cmd_check_int(opt->dup_window, 0, INT32_MAX, "-W");
  fmap_error_cmd_check_int(opt->aln_output_mode_ind, 0, 1, "-I");

  if(0 == opt->algos[0]) {
      fmap_error("no algorithms given for stage 1", Exit, CommandLineArgument);
  }

  if(FMAP_FILE_BZ2_COMPRESSION == opt->output_compr 
     && -1 == opt->reads_queue_size) {
      fmap_error("cannot buffer reads with bzip2 output (options \"-q 1 -J\")", Exit, OutOfRange);
  }

  for(i=0;i<2;i++) {
      // check mapping algorithm specific options
      fmap_map1_opt_check(opt->opt_map1[i]);
      fmap_map2_opt_check(opt->opt_map2[i]);
      fmap_map3_opt_check(opt->opt_map3[i]);

      // check that common values match other opt values
      __fmap_map_all_opts_check_common1(opt, opt->opt_map1[i]);
      __fmap_map_all_opts_check_common2(opt, opt->opt_map2[i]);
      __fmap_map_all_opts_check_common2(opt, opt->opt_map3[i]);
  }
}

int 
fmap_map_all_main(int argc, char *argv[])
{
  fmap_map_all_opt_t *opt = NULL;

  // random seed
  srand48(0); 

  // init opt
  opt = fmap_map_all_opt_init();
      
  // get options
  if(1 != fmap_map_all_opt_parse(argc, argv, opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      fprintf(stderr, "argc=%d optind=%d\n", argc, optind);
      return fmap_map_all_usage(opt);
  }
  else { 
      // check command line arguments
      fmap_map_all_opt_check(opt);
  }

  // run map_all
  fmap_map_all_core(opt);

  // destroy opt
  fmap_map_all_opt_destroy(opt);

  fmap_progress_print2("terminating successfully");

  return 0;
}
