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
#include "fmap_map_util.h"
#include "fmap_map1.h"
#include "fmap_map1_aux.h"
#include "fmap_map2.h"
#include "fmap_map2_aux.h"
#include "fmap_map3.h"
#include "fmap_map3_aux.h"
#include "fmap_map_all.h"

/* Notes:
   We could avoid some computation give we know which algorithms to 
   run. This includes stacks, memory pools, as well as not loading 
   in all reference data.
   */

// sort by min-seqid, min-position, max-score
#define __fmap_map_sam_sort_lt(a, b) ( ((a).seqid < (b).seqid \
                                            || ( (a).seqid == (b).seqid && (a).pos < (b).pos ) \
                                            || ( (a).seqid == (b).seqid && (a).pos == (b).pos && (a).score < (b).score )) \
                                          ? 1 : 0 )

FMAP_SORT_INIT(fmap_map_sam_t, fmap_map_sam_t, __fmap_map_sam_sort_lt)

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t fmap_map_all_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t fmap_map_all_read_lock_low = 0;
#define FMAP_MAP_ALL_THREAD_BLOCK_SIZE 512
#endif

static inline void
fmap_map_all_mapq(fmap_map_sams_t *sams, fmap_map_opt_t *opt)
{
  int32_t i;
  int32_t n_best = 0, n_best_subo = 0;
  int32_t best_score, cur_score, best_subo;
  int32_t mapq;
  int32_t stage = -1;
  int32_t algo_id = FMAP_MAP_ALGO_NONE; 

  // estimate mapping quality TODO: this needs to be refined
  best_score = INT32_MIN;
  best_subo = INT32_MIN+1;
  n_best = n_best_subo = 0;
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      if(best_score < cur_score) {
          // save sub-optimal
          best_subo = best_score;
          n_best_subo = n_best;
          // update
          best_score = cur_score;
          n_best = 1;
          stage = (algo_id == FMAP_MAP_ALGO_NONE) ? sams->sams[i].algo_stage : -1;
          algo_id = (algo_id == FMAP_MAP_ALGO_NONE) ? sams->sams[i].algo_id : -1;
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
      if(FMAP_MAP_ALGO_MAP2 == sams->sams[i].algo_id
         || FMAP_MAP_ALGO_MAP3 == sams->sams[i].algo_id) {
          cur_score = sams->sams[i].score_subo;
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
            case FMAP_MAP_ALGO_MAP1:
              // what is the best value for map1
              if(0 < opt->opt_map1[stage]->seed_length) {
                  best_subo = opt->score_match * (opt->opt_map1[stage]->seed_length - opt->opt_map1[stage]->seed_max_mm);
              }
              else {
                  best_subo = 0;
              }
              break;
            case FMAP_MAP_ALGO_MAP2:
              best_subo = opt->opt_map2[stage]->score_thr; break;
            case FMAP_MAP_ALGO_MAP3:
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
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      if(cur_score == best_score) {
          sams->sams[i].mapq = mapq;
      }
      else {
          sams->sams[i].mapq = 0;
      }
  }
}

static void
fmap_map_all_sams_merge_helper(fmap_map_sams_t *dest, fmap_map_sams_t *src, int32_t stage)
{
  int32_t i, n;

  // make room
  n = dest->n;
  fmap_map_sams_realloc(dest, dest->n + src->n);
  // copy over
  for(i=0;i<src->n;i++) {
      // nullify
      fmap_map_sam_copy_and_nullify(&dest->sams[i+n], &src->sams[i]);
      // copy over the stage #
      dest->sams[i+n].algo_stage = stage;
  }
}

static void
fmap_map_all_remove_duplicates(fmap_map_sams_t *sams, int32_t dup_window)
{
  int32_t i, j, end, best_score_i;
  // sort
  fmap_sort_introsort(fmap_map_sam_t, sams->n, sams->sams);
  
  // remove duplicates within a window
  for(i=j=0;i<sams->n;) {

      // get the change
      end = best_score_i = i;
      while(end+1 < sams->n) {
          if(sams->sams[end].seqid == sams->sams[end+1].seqid
             && fabs(sams->sams[end].pos - sams->sams[end+1].pos) <= dup_window) {
              // track the best scoring
              if(sams->sams[best_score_i].score < sams->sams[end].score) {
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
          // destroy
          fmap_map_sam_destroy(&sams->sams[j]);
          // nullify
          fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[best_score_i]);
      }

      // next
      i = end+1;
      j++;
  }

  // destroy the sams
  for(i=j;i<sams->n;i++) {
      fmap_map_sam_destroy(&sams->sams[i]);
  }
  fmap_map_sams_realloc(sams, j);
}

static fmap_map_sams_t *
fmap_map_all_sams_merge(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_bwt_t *bwt[2], fmap_sa_t *sa[2],
                       fmap_map_sams_t *sams_map1, fmap_map_sams_t *sams_map2, fmap_map_sams_t *sams_map3,
                       int32_t stage, fmap_map_opt_t *opt)
{
  fmap_map_sams_t *sams = NULL;

  // init
  sams = fmap_map_sams_init();

  // merge
  fmap_map_all_sams_merge_helper(sams, sams_map1, stage);
  fmap_map_all_sams_merge_helper(sams, sams_map2, stage);
  fmap_map_all_sams_merge_helper(sams, sams_map3, stage);

  // remove duplicates
  fmap_map_all_remove_duplicates(sams, opt->dup_window);

  // mapping quality
  fmap_map_all_mapq(sams, opt);

  // choose alignment(s)
  if(0 == opt->aln_output_mode_ind) {
      // consider all algos together
      fmap_map_sams_filter1(sams, opt->aln_output_mode, FMAP_MAP_ALGO_NONE);
  }
  else {
      // consider all algos independently
      fmap_map_sams_filter1(sams, opt->aln_output_mode, FMAP_MAP_ALGO_MAP1);
      fmap_map_sams_filter1(sams, opt->aln_output_mode, FMAP_MAP_ALGO_MAP2);
      fmap_map_sams_filter1(sams, opt->aln_output_mode, FMAP_MAP_ALGO_MAP3);
  }

  return sams;
}

static void
fmap_map_all_core_worker(fmap_seq_t **seq_buffer, fmap_map_sams_t **sams, int32_t seq_buffer_length, 
                         fmap_refseq_t *refseq, fmap_bwt_t *bwt[2], fmap_sa_t *sa[2],
                         int32_t tid, fmap_map_opt_t *opt)
{
  int32_t low = 0, high, i, j;
  // map1 
  fmap_bwt_match_width_t *width_map1[2][2], *seed_width_map1[2][2];
  int32_t width_length_map1[2] = {0, 0};
  fmap_map1_aux_stack_t *stack_map1=NULL;
  fmap_map_sams_t *sams_map1;
  // map2
  fmap_map2_global_mempool_t *pool_map2 = NULL;
  fmap_map_sams_t *sams_map2;
  // map3
  fmap_map_sams_t *sams_map3;
  uint8_t *flow_map3[2][2];

  // map1
  for(i=0;i<2;i++) {
      if(opt->algos[i] & FMAP_MAP_ALGO_MAP1) {
          width_map1[i][0] = width_map1[i][1] = NULL;
          width_length_map1[i] = 0;
          seed_width_map1[i][0] = fmap_calloc(opt->opt_map1[i]->seed_length, sizeof(fmap_bwt_match_width_t), "seed_width[0]");
          seed_width_map1[i][1] = fmap_calloc(opt->opt_map1[i]->seed_length, sizeof(fmap_bwt_match_width_t), "seed_width[1]");
          if(NULL == stack_map1) {
              stack_map1 = fmap_map1_aux_stack_init();
          }
      }
      // map2
      if(opt->algos[i] & FMAP_MAP_ALGO_MAP2) {
          if(NULL == pool_map2) {
              pool_map2 = fmap_map2_global_mempool_init();
          }
      }
      // map3
      if(opt->algos[i] & FMAP_MAP_ALGO_MAP3) {
          // set up the flow order
          if(0 < seq_buffer_length && 0 < opt->opt_map3[i]->hp_diff) {
              if(FMAP_SEQ_TYPE_SFF == seq_buffer[0]->type) {
                  flow_map3[i][0] = fmap_malloc(sizeof(uint8_t) * 4, "flow[i][0]");
                  flow_map3[i][1] = fmap_malloc(sizeof(uint8_t) * 4, "flow[i][1]");
                  for(j=0;j<4;j++) {
                      flow_map3[i][0][j] = fmap_nt_char_to_int[(int)seq_buffer[0]->data.sff->gheader->flow->s[j]];
                      flow_map3[i][1][j] = 3 - fmap_nt_char_to_int[(int)seq_buffer[0]->data.sff->gheader->flow->s[j]];
                  }
              }
              else {
                  fmap_error("bug encountered", Exit, OutOfRange);
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
          fmap_map_opt_t opt_local_map1[2];

          // map1
          for(i=0;i<2;i++) {
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP1) {
                  opt_local_map1[i] = (*opt->opt_map1[i]); // copy over values
              }
              // map2
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP2) {
                  // - none
              }
              // map3
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP3) {
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
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP1) {
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
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP2) {
                  // - none
              }
              // map3
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP3) {
                  // - none
              }
          }

          // run the core algorithms
          for(i=0;i<2;i++) {
              // fmap_map1_aux_core
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP1) {
                  sams_map1 = fmap_map1_aux_core(seq, refseq, bwt[1], sa[1], width_map1[i], 
                                                (0 < opt_local_map1[i].seed_length) ? seed_width_map1[i] : NULL, 
                                                &opt_local_map1[i], stack_map1);
                  // adjust map1 scoring, since it does not consider opt->score_match
                  fmap_map_util_map1_adjust_score(sams_map1, opt->score_match, opt->pen_mm, opt->pen_gapo, opt->pen_gape);
              }
              else {
                  // empty
                  sams_map1 = fmap_map_sams_init();
              }
              // fmap_map2_aux_core
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP2) {
                  sams_map2 = fmap_map2_aux_core(opt->opt_map2[i], seq_char, refseq, bwt, sa, pool_map2);
              }
              else {
                  // empty
                  sams_map2 = fmap_map_sams_init();
              }
              // fmap_map3_aux_core
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP3) {
                  sams_map3 = fmap_map3_aux_core(seq, flow_map3[i], refseq, bwt[1], sa[1], opt->opt_map3[i]);
              }
              else {
                  // empty
                  sams_map3 = fmap_map_sams_init();
              }

              // consolidate mappings
              sams[low] = fmap_map_all_sams_merge(orig_seq, refseq, bwt, sa,
                                                 sams_map1, sams_map2, sams_map3, 
                                                 i+1, opt);

              // re-align the alignments in flow-space
              if(FMAP_SEQ_TYPE_SFF == seq_buffer[low]->type) {
                  fmap_map_util_fsw(seq_buffer[low]->data.sff,
                                    sams[low], refseq,
                                    opt->bw, opt->aln_global, INT32_MIN,
                                    opt->score_match, opt->pen_mm, opt->pen_gapo,
                                    opt->pen_gape, opt->fscore);
              }


              // destroy
              // map1
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP1) {
                  fmap_map_sams_destroy(sams_map1);
                  sams_map1 = NULL;
              }
              // map2
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP2) {
                  fmap_map_sams_destroy(sams_map2);
                  sams_map2 = NULL;
              }
              // map3
              if(opt->algos[i] & FMAP_MAP_ALGO_MAP3) {
                  fmap_map_sams_destroy(sams_map3);
                  sams_map3 = NULL;
              }

              // check if we found any mappings
              if(0 < sams[low]->n) {
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
      if(opt->algos[i] & FMAP_MAP_ALGO_MAP1) {
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
      if(opt->algos[i] & FMAP_MAP_ALGO_MAP2) {
          if(NULL != pool_map2) {
              fmap_map2_global_mempool_destroy(pool_map2);
              pool_map2 = NULL;
          }
      }
      // map3
      if(opt->algos[i] & FMAP_MAP_ALGO_MAP3) {
          free(flow_map3[i][0]);
          free(flow_map3[i][1]);
      }
  }
}

static void *
fmap_map_all_core_thread_worker(void *arg)
{
  fmap_map_all_thread_data_t *thread_data = (fmap_map_all_thread_data_t*)arg;

  fmap_map_all_core_worker(thread_data->seq_buffer, thread_data->sams, thread_data->seq_buffer_length, 
                           thread_data->refseq, thread_data->bwt, thread_data->sa, 
                           thread_data->tid, thread_data->opt);

  return arg;
}

static void 
fmap_map_all_core(fmap_map_opt_t *opt)
{
  uint32_t i, n_reads_processed=0;
  int32_t seq_buffer_length;
  fmap_refseq_t *refseq=NULL;
  fmap_bwt_t *bwt[2]={NULL, NULL};
  fmap_sa_t *sa[2]={NULL, NULL};
  fmap_file_t *fp_reads=NULL;
  fmap_seq_io_t *seqio = NULL;
  fmap_seq_t **seq_buffer = NULL;
  fmap_map_sams_t **sams = NULL;
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
      if(opt->algos[i] & FMAP_MAP_ALGO_MAP1) {
          // - none
      }
      // map2
      if(opt->algos[i] & FMAP_MAP_ALGO_MAP2) {
          opt->opt_map2[i]->score_thr *= opt->score_match;
          opt->opt_map2[i]->length_coef *= opt->score_match;
      }
      // map3
      if(opt->algos[i] & FMAP_MAP_ALGO_MAP3) {
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
  sams = fmap_malloc(sizeof(fmap_map_sams_t*)*reads_queue_size, "sams");

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
  fmap_sam_print_header(fmap_file_stdout, refseq, seqio, opt->sam_rg, opt->sam_sff_tags, opt->argc, opt->argv);

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
          fmap_map_all_core_worker(seq_buffer, sams, seq_buffer_length, refseq, bwt, sa, 0, opt);
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
              thread_data[i].sams = sams;
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
      fmap_map_all_core_worker(seq_buffer, sams, seq_buffer_length, refseq, bwt, sa, 0, opt);
#endif

      if(-1 != opt->reads_queue_size) {
          fmap_progress_print("writing alignments");
      }
      for(i=0;i<seq_buffer_length;i++) {
          // write
          fmap_map_sams_print(seq_buffer[i], refseq, sams[i], opt->sam_sff_tags);

          // free alignments
          fmap_map_sams_destroy(sams[i]);
          sams[i] = NULL;
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
  free(sams);
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

// for map1/map2/map3
#define __fmap_map_all_opts_copy1(opt_map_all, opt_map_other) do { \
    (opt_map_other)->fn_fasta = fmap_strdup((opt_map_all)->fn_fasta); \
    (opt_map_other)->fn_reads = fmap_strdup((opt_map_all)->fn_reads); \
    (opt_map_other)->reads_format = (opt_map_all)->reads_format; \
    (opt_map_other)->pen_mm = (opt_map_all)->pen_mm; \
    (opt_map_other)->pen_gapo = (opt_map_all)->pen_gapo; \
    (opt_map_other)->pen_gape = (opt_map_all)->pen_gape; \
    (opt_map_other)->reads_queue_size = (opt_map_all)->reads_queue_size; \
    (opt_map_other)->num_threads = (opt_map_all)->num_threads; \
    (opt_map_other)->aln_output_mode = FMAP_MAP_UTIL_ALN_MODE_ALL; \
    (opt_map_other)->sam_rg = fmap_strdup((opt_map_all)->sam_rg); \
    (opt_map_other)->input_compr = (opt_map_all)->input_compr; \
    (opt_map_other)->output_compr = (opt_map_all)->output_compr; \
    (opt_map_other)->shm_key = (opt_map_all)->shm_key; \
} while(0)

// for map2 and map3
#define __fmap_map_all_opts_copy2(opt_map_all, opt_map_other) do { \
    __fmap_map_all_opts_copy1(opt_map_all, opt_map_other); \
    (opt_map_other)->score_match = (opt_map_all)->score_match; \
    (opt_map_other)->bw = (opt_map_all)->bw; \
    (opt_map_other)->aln_global = (opt_map_all)->aln_global; \
} while(0)

int32_t
fmap_map_all_opt_parse(int argc, char *argv[], fmap_map_opt_t *opt)
{
  int32_t i, j, start, opt_type, opt_type_next, opt_stage, opt_stage_next;

  opt->argc = argc; opt->argv = argv;

  // parse common options as well as map1/map2/map3 commands
  start = 0;
  i = 1;
  opt_type = opt_type_next = FMAP_MAP_ALGO_NONE;
  opt_stage = opt_stage_next = 0;
  while(i<argc) {
      if(0 == strcmp("map1", argv[i])) {
          opt_type_next = FMAP_MAP_ALGO_MAP1;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
      }
      else if(0 == strcmp("map2", argv[i])) {
          opt_type_next = FMAP_MAP_ALGO_MAP2;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
      }
      else if(0 == strcmp("map3", argv[i])) {
          opt_type_next = FMAP_MAP_ALGO_MAP3;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
      }
      else if(0 == strcmp("MAP1", argv[i])) {
          opt_type_next = FMAP_MAP_ALGO_MAP1;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
      }
      else if(0 == strcmp("MAP2", argv[i])) {
          opt_type_next = FMAP_MAP_ALGO_MAP2;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
      }
      else if(0 == strcmp("MAP3", argv[i])) {
          opt_type_next = FMAP_MAP_ALGO_MAP3;
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
            case FMAP_MAP_ALGO_NONE:
              // parse common options
              if(0 == fmap_map_opt_parse(i-start, argv+start, opt)) {
                  return 0;
              }
              // copy over common values into the other opts
              for(j=0;j<2;j++) {
                  __fmap_map_all_opts_copy1(opt, opt->opt_map1[j]);
                  __fmap_map_all_opts_copy2(opt, opt->opt_map2[j]);
                  __fmap_map_all_opts_copy2(opt, opt->opt_map3[j]);
              }
              break;
            case FMAP_MAP_ALGO_MAP1:
              if(0 < i - start) {
                  if(0 == fmap_map_opt_parse(i-start, argv+start, opt->opt_map1[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            case FMAP_MAP_ALGO_MAP2:
              // parse map2 options
              if(0 < i - start) {
                  if(0 == fmap_map_opt_parse(i-start, argv+start, opt->opt_map2[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            case FMAP_MAP_ALGO_MAP3:
              // parse map3 options
              if(0 < i - start) {
                  if(0 == fmap_map_opt_parse(i-start, argv+start, opt->opt_map3[opt_stage-1])) {
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

int 
fmap_map_all_main(int argc, char *argv[])
{
  fmap_map_opt_t *opt = NULL;

  // random seed
  srand48(0); 

  // init opt
  opt = fmap_map_opt_init(FMAP_MAP_ALGO_MAPALL);
      
  // get options
  if(1 != fmap_map_all_opt_parse(argc, argv, opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return fmap_map_opt_usage(opt);
  }
  else { 
      // check command line arguments
      fmap_map_opt_check(opt);
  }

  // run map_all
  fmap_map_all_core(opt);

  // destroy opt
  fmap_map_opt_destroy(opt);

  fmap_progress_print2("terminating successfully");

  return 0;
}
