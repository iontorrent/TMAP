/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <config.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#include <unistd.h>
#endif
#include <unistd.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sam.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt_gen.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_bwt_match.h"
#include "../index/tmap_sa.h"
#include "../io/tmap_seq_io.h"
#include "../server/tmap_shm.h"
#include "tmap_map_util.h"
#include "tmap_map1.h"
#include "tmap_map1_aux.h"
#include "tmap_map2.h"
#include "tmap_map2_aux.h"
#include "tmap_map3.h"
#include "tmap_map3_aux.h"
#include "tmap_map_all.h"

/* Notes:
   We could avoid some computation give we know which algorithms to 
   run. This includes stacks, memory pools, as well as not loading 
   in all reference data.
   */

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t tmap_map_all_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t tmap_map_all_read_lock_low = 0;
#define TMAP_MAP_ALL_THREAD_BLOCK_SIZE 512
#endif

static inline void
tmap_map_all_mapq(tmap_map_sams_t *sams, tmap_map_opt_t *opt)
{
  int32_t i;
  int32_t n_best = 0, n_best_subo = 0;
  int32_t best_score, cur_score, best_subo;
  int32_t mapq;
  int32_t stage = -1;
  int32_t algo_id = TMAP_MAP_ALGO_NONE; 

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
          stage = (algo_id == TMAP_MAP_ALGO_NONE) ? sams->sams[i].algo_stage-1 : -1;
          algo_id = (algo_id == TMAP_MAP_ALGO_NONE) ? sams->sams[i].algo_id : -1;
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
      if(TMAP_MAP_ALGO_MAP2 == sams->sams[i].algo_id
         || TMAP_MAP_ALGO_MAP3 == sams->sams[i].algo_id) {
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
            case TMAP_MAP_ALGO_MAP1:
              // what is the best value for map1
              if(0 < opt->opt_map1[stage]->seed_length) {
                  best_subo = opt->score_match * (opt->opt_map1[stage]->seed_length - opt->opt_map1[stage]->seed_max_mm);
              }
              else {
                  best_subo = 0;
              }
              break;
            case TMAP_MAP_ALGO_MAP2:
              best_subo = opt->opt_map2[stage]->score_thr; break;
            case TMAP_MAP_ALGO_MAP3:
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
tmap_map_all_sams_merge_helper(tmap_map_sams_t *dest, tmap_map_sams_t *src, int32_t stage)
{
  int32_t i, n;

  if(0 == src->n) return;

  // make room
  n = dest->n;
  tmap_map_sams_realloc(dest, dest->n + src->n);
  // copy over
  for(i=0;i<src->n;i++) {
      // nullify
      tmap_map_sam_copy_and_nullify(&dest->sams[i+n], &src->sams[i]);
      // copy over the stage #
      dest->sams[i+n].algo_stage = stage;
  }
}

static tmap_map_sams_t *
tmap_map_all_sams_merge(tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2],
                       tmap_map_sams_t *sams_map1, tmap_map_sams_t *sams_map2, tmap_map_sams_t *sams_map3,
                       int32_t stage, tmap_map_opt_t *opt)
{
  tmap_map_sams_t *sams = NULL;

  // init
  sams = tmap_map_sams_init();

  // merge
  tmap_map_all_sams_merge_helper(sams, sams_map1, stage);
  tmap_map_all_sams_merge_helper(sams, sams_map2, stage);
  tmap_map_all_sams_merge_helper(sams, sams_map3, stage);

  // no alignments
  if(0 == sams->n) return sams;

  // remove duplicates
  tmap_map_util_remove_duplicates(sams, opt->dup_window);

  // mapping quality
  tmap_map_all_mapq(sams, opt);

  // choose alignment(s)
  if(0 == opt->aln_output_mode_ind) {
      // consider all algos together
      tmap_map_sams_filter1(sams, opt->aln_output_mode, TMAP_MAP_ALGO_NONE);
  }
  else {
      // consider all algos independently
      tmap_map_sams_filter1(sams, opt->aln_output_mode, TMAP_MAP_ALGO_MAP1);
      tmap_map_sams_filter1(sams, opt->aln_output_mode, TMAP_MAP_ALGO_MAP2);
      tmap_map_sams_filter1(sams, opt->aln_output_mode, TMAP_MAP_ALGO_MAP3);
  }

  return sams;
}

static void
tmap_map_all_core_worker(tmap_seq_t **seq_buffer, tmap_map_sams_t **sams, int32_t seq_buffer_length, 
                         tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2],
                         int32_t tid, tmap_map_opt_t *opt)
{
  int32_t low = 0, high, i, j;
  // map1 
  tmap_bwt_match_width_t *width_map1[2][2], *seed_width_map1[2][2];
  int32_t width_length_map1[2] = {0, 0}, seed2_len_map1[2];
  tmap_map1_aux_stack_t *stack_map1=NULL;
  tmap_map_sams_t *sams_map1 = NULL;
  // map2
  tmap_map2_global_mempool_t *pool_map2 = NULL;
  tmap_map_sams_t *sams_map2 = NULL;
  // map3
  tmap_map_sams_t *sams_map3 = NULL;
  uint8_t *flow_map3[2][2];

  // map1
  for(i=0;i<opt->num_stages;i++) {
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP1) {
          width_map1[i][0] = width_map1[i][1] = NULL;
          width_length_map1[i] = 0;
          seed_width_map1[i][0] = tmap_calloc(opt->opt_map1[i]->seed_length, sizeof(tmap_bwt_match_width_t), "seed_width[0]");
          seed_width_map1[i][1] = tmap_calloc(opt->opt_map1[i]->seed_length, sizeof(tmap_bwt_match_width_t), "seed_width[1]");
          if(NULL == stack_map1) {
              stack_map1 = tmap_map1_aux_stack_init();
          }
      }
      // map2
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP2) {
          if(NULL == pool_map2) {
              pool_map2 = tmap_map2_global_mempool_init();
          }
      }
      // map3
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP3) {
          // set up the flow order
          if(0 < seq_buffer_length && 0 < opt->opt_map3[i]->hp_diff) {
              if(TMAP_SEQ_TYPE_SFF == seq_buffer[0]->type) {
                  flow_map3[i][0] = tmap_malloc(sizeof(uint8_t) * 4, "flow[i][0]");
                  flow_map3[i][1] = tmap_malloc(sizeof(uint8_t) * 4, "flow[i][1]");
                  for(j=0;j<4;j++) {
                      flow_map3[i][0][j] = tmap_nt_char_to_int[(int)seq_buffer[0]->data.sff->gheader->flow->s[j]];
                      flow_map3[i][1][j] = 3 - tmap_nt_char_to_int[(int)seq_buffer[0]->data.sff->gheader->flow->s[j]];
                  }
              }
              else {
                  tmap_error("bug encountered", Exit, OutOfRange);
              }
          }
      }
  }

  while(low < seq_buffer_length) {
#ifdef HAVE_LIBPTHREAD
      if(1 < opt->num_threads) {
          pthread_mutex_lock(&tmap_map_all_read_lock);

          // update bounds
          low = tmap_map_all_read_lock_low;
          tmap_map_all_read_lock_low += TMAP_MAP_ALL_THREAD_BLOCK_SIZE;
          high = low + TMAP_MAP_ALL_THREAD_BLOCK_SIZE;
          if(seq_buffer_length < high) {
              high = seq_buffer_length; 
          }

          pthread_mutex_unlock(&tmap_map_all_read_lock);
      }
      else {
          high = seq_buffer_length; // process all
      }
#else 
      high = seq_buffer_length; // process all
#endif
      while(low<high) {
          tmap_seq_t *seq[2]={NULL, NULL}, *orig_seq=NULL, *seq_char=NULL;
          orig_seq = seq_buffer[low];
          tmap_string_t *bases[2]={NULL, NULL};
          tmap_map_opt_t opt_local_map1[2];

          // map1
          for(i=0;i<opt->num_stages;i++) {
              if(opt->algos[i] & TMAP_MAP_ALGO_MAP1) {
                  opt_local_map1[i] = (*opt->opt_map1[i]); // copy over values
              }
              // map2
              if(opt->algos[i] & TMAP_MAP_ALGO_MAP2) {
                  // - none
              }
              // map3
              if(opt->algos[i] & TMAP_MAP_ALGO_MAP3) {
                  // - none
              }
          }

          // clone the sequence 
          seq[0] = tmap_seq_clone(orig_seq);
          seq[1] = tmap_seq_clone(orig_seq);
          seq_char = tmap_seq_clone(orig_seq);

          // Adjust for SFF
          tmap_seq_remove_key_sequence(seq[0]);
          tmap_seq_remove_key_sequence(seq[1]);
          tmap_seq_remove_key_sequence(seq_char);

          // reverse compliment
          tmap_seq_reverse_compliment(seq[1]);

          // convert to integers
          tmap_seq_to_int(seq[0]);
          tmap_seq_to_int(seq[1]);

          // get bases
          bases[0] = tmap_seq_get_bases(seq[0]);
          bases[1] = tmap_seq_get_bases(seq[1]);

          // map1
          for(i=0;i<opt->num_stages;i++) {
              if((opt->algos[i] & TMAP_MAP_ALGO_MAP1)
                 && (opt->opt_map1[i]->seed_length < 0 || opt->opt_map1[i]->seed_length <= bases[0]->l)) {
                  if(opt_local_map1[i].seed2_length < 0 || bases[0]->l < opt_local_map1[i].seed2_length) {
                      seed2_len_map1[i] = bases[0]->l;
                  }
                  else {
                      seed2_len_map1[i] = opt_local_map1[i].seed2_length;
                  }

                  opt_local_map1[i].max_mm = (opt->opt_map1[i]->max_mm < 0) ? (int)(0.99 + opt->opt_map1[i]->max_mm_frac * seed2_len_map1[i]) : opt->opt_map1[i]->max_mm;
                  opt_local_map1[i].max_gape = (opt->opt_map1[i]->max_gape < 0) ? (int)(0.99 + opt->opt_map1[i]->max_gape_frac * seed2_len_map1[i]) : opt->opt_map1[i]->max_gape;
                  opt_local_map1[i].max_gapo = (opt->opt_map1[i]->max_gapo < 0) ? (int)(0.99 + opt->opt_map1[i]->max_gapo_frac * seed2_len_map1[i]) : opt->opt_map1[i]->max_gapo;
                  if(width_length_map1[i] < seed2_len_map1[i]) {
                      free(width_map1[i][0]); free(width_map1[i][1]);
                      width_length_map1[i] = seed2_len_map1[i];
                      width_map1[i][0] = tmap_calloc(width_length_map1[i], sizeof(tmap_bwt_match_width_t), "width[0]");
                      width_map1[i][1] = tmap_calloc(width_length_map1[i], sizeof(tmap_bwt_match_width_t), "width[1]");
                  }
                  tmap_bwt_match_cal_width(bwt[0], seed2_len_map1[i], bases[0]->s, width_map1[i][0]);
                  tmap_bwt_match_cal_width(bwt[0], seed2_len_map1[i], bases[1]->s, width_map1[i][1]);

                  if(seed2_len_map1[i] < opt->opt_map1[i]->seed_length) {
                      opt_local_map1[i].seed_length = -1;
                  }
                  else {
                      tmap_bwt_match_cal_width(bwt[0], opt->opt_map1[i]->seed_length, bases[0]->s, seed_width_map1[i][0]);
                      tmap_bwt_match_cal_width(bwt[0], opt->opt_map1[i]->seed_length, bases[1]->s, seed_width_map1[i][1]);
                  }
              }
              // map2
              if(opt->algos[i] & TMAP_MAP_ALGO_MAP2) {
                  // - none
              }
              // map3
              if(opt->algos[i] & TMAP_MAP_ALGO_MAP3) {
                  // - none
              }
          }

          // run the core algorithms
          for(i=0;i<opt->num_stages;i++) {
              // tmap_map1_aux_core
              if((opt->algos[i] & TMAP_MAP_ALGO_MAP1)
                 && (opt->opt_map1[i]->seed_length < 0 || opt->opt_map1[i]->seed_length <= bases[0]->l)) {
                  // TODO: seed2_length
                  sams_map1 = tmap_map1_aux_core(seq, refseq, bwt[1], sa[1], width_map1[i], 
                                                (0 < opt_local_map1[i].seed_length) ? seed_width_map1[i] : NULL, 
                                                &opt_local_map1[i], stack_map1, seed2_len_map1[i]);
                  // adjust map1 scoring, since it does not consider opt->score_match
                  tmap_map_util_map1_adjust_score(sams_map1, opt->score_match, opt->pen_mm, opt->pen_gapo, opt->pen_gape);
              }
              else {
                  // empty
                  sams_map1 = tmap_map_sams_init();
              }
              // tmap_map2_aux_core
              if(opt->algos[i] & TMAP_MAP_ALGO_MAP2) {
                  sams_map2 = tmap_map2_aux_core(opt->opt_map2[i], seq_char, refseq, bwt, sa, pool_map2);
              }
              else {
                  // empty
                  sams_map2 = tmap_map_sams_init();
              }
              // tmap_map3_aux_core
              if(opt->algos[i] & TMAP_MAP_ALGO_MAP3) {
                  sams_map3 = tmap_map3_aux_core(seq, flow_map3[i], refseq, bwt[1], sa[1], opt->opt_map3[i]);
              }
              else {
                  // empty
                  sams_map3 = tmap_map_sams_init();
              }

              // consolidate mappings
              sams[low] = tmap_map_all_sams_merge(orig_seq, refseq, bwt, sa,
                                                 sams_map1, sams_map2, sams_map3, 
                                                 i+1, opt);

              // re-align the alignments in flow-space
              if(TMAP_SEQ_TYPE_SFF == seq_buffer[low]->type) {
                  tmap_map_util_fsw(seq_buffer[low]->data.sff,
                                    sams[low], refseq,
                                    opt->bw, opt->aln_global, INT32_MIN,
                                    opt->score_match, opt->pen_mm, opt->pen_gapo,
                                    opt->pen_gape, opt->fscore);
              }
              // destroy
              // map1
              tmap_map_sams_destroy(sams_map1);
              sams_map1 = NULL;
              // map2
              tmap_map_sams_destroy(sams_map2);
              sams_map2 = NULL;
              // map3
              tmap_map_sams_destroy(sams_map3);
              sams_map3 = NULL;

              // check if we found any mappings
              if(0 < sams[low]->n) {
                  // yes we did
                  break;
              }
              else if(0 == i && 2 == opt->num_stages) {
                  tmap_map_sams_destroy(sams[low]);
              }
          }

          // destroy
          tmap_seq_destroy(seq[0]);
          tmap_seq_destroy(seq[1]);
          tmap_seq_destroy(seq_char);

          // next
          low++;
      }
  }

  // map1
  for(i=0;i<opt->num_stages;i++) {
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP1) {
          if(NULL != stack_map1) {
              tmap_map1_aux_stack_destroy(stack_map1);
              stack_map1 = NULL;
          }
          free(seed_width_map1[i][0]);
          free(seed_width_map1[i][1]);
          free(width_map1[i][0]);
          free(width_map1[i][1]);
      }
      // map2
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP2) {
          if(NULL != pool_map2) {
              tmap_map2_global_mempool_destroy(pool_map2);
              pool_map2 = NULL;
          }
      }
      // map3
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP3) {
          if(0 < seq_buffer_length && 0 < opt->opt_map3[i]->hp_diff) {
              if(TMAP_SEQ_TYPE_SFF == seq_buffer[0]->type) {
                  free(flow_map3[i][0]);
                  free(flow_map3[i][1]);
              }
          }
      }
  }
}

static void *
tmap_map_all_core_thread_worker(void *arg)
{
  tmap_map_all_thread_data_t *thread_data = (tmap_map_all_thread_data_t*)arg;

  tmap_map_all_core_worker(thread_data->seq_buffer, thread_data->sams, thread_data->seq_buffer_length, 
                           thread_data->refseq, thread_data->bwt, thread_data->sa, 
                           thread_data->tid, thread_data->opt);

  return arg;
}

static void 
tmap_map_all_core(tmap_map_opt_t *opt)
{
  uint32_t i, n_reads_processed=0;
  int32_t seq_buffer_length;
  tmap_refseq_t *refseq=NULL;
  tmap_bwt_t *bwt[2]={NULL, NULL};
  tmap_sa_t *sa[2]={NULL, NULL};
  tmap_file_t *fp_reads=NULL;
  tmap_seq_io_t *seqio = NULL;
  tmap_seq_t **seq_buffer = NULL;
  tmap_map_sams_t **sams = NULL;
  tmap_shm_t *shm = NULL;
  int32_t reads_queue_size;
  
  if(NULL == opt->fn_reads) {
      tmap_progress_set_verbosity(0); 
  }

  if(0 == opt->shm_key) {
      tmap_progress_print("reading in reference data");
      refseq = tmap_refseq_read(opt->fn_fasta, 0);
      bwt[0] = tmap_bwt_read(opt->fn_fasta, 0);
      bwt[1] = tmap_bwt_read(opt->fn_fasta, 1);
      sa[0] = tmap_sa_read(opt->fn_fasta, 0);
      sa[1] = tmap_sa_read(opt->fn_fasta, 1);
      tmap_progress_print2("reference data read in");
  }
  else {
      tmap_progress_print("retrieving reference data from shared memory");
      shm = tmap_shm_init(opt->shm_key, 0, 0);
      if(NULL == (refseq = tmap_refseq_shm_unpack(tmap_shm_get_buffer(shm, TMAP_SHM_LISTING_REFSEQ)))) {
          tmap_error("the packed reference sequence was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (bwt[0] = tmap_bwt_shm_unpack(tmap_shm_get_buffer(shm, TMAP_SHM_LISTING_BWT)))) {
          tmap_error("the BWT string was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (bwt[1] = tmap_bwt_shm_unpack(tmap_shm_get_buffer(shm, TMAP_SHM_LISTING_REV_BWT)))) {
          tmap_error("the reverse BWT string was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (sa[0] = tmap_sa_shm_unpack(tmap_shm_get_buffer(shm, TMAP_SHM_LISTING_SA)))) {
          tmap_error("the SA was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (sa[1] = tmap_sa_shm_unpack(tmap_shm_get_buffer(shm, TMAP_SHM_LISTING_REV_SA)))) {
          tmap_error("the reverse SA was not found in shared memory", Exit, SharedMemoryListing);
      }
      tmap_progress_print2("reference data retrieved from shared memory");

  }

  // TODO: this could be dangerous if algorithms change, needs to be refactored
  // adjust mapping algorithm specific options
  for(i=0;i<opt->num_stages;i++) {
      // map1
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP1) {
          // - none
      }
      // map2
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP2) {
          opt->opt_map2[i]->score_thr *= opt->score_match;
          opt->opt_map2[i]->length_coef *= opt->score_match;
      }
      // map3
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP3) {
          opt->opt_map3[i]->score_thr *= opt->score_match;
          if(-1 == opt->opt_map3[i]->seed_length) {
              opt->opt_map3[i]->seed_length = tmap_map3_get_seed_length(refseq->len);
              tmap_progress_print("setting the seed length to %d for map3 stage %d",
                                  opt->opt_map3[i]->seed_length,
                                  i+1);
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
  seq_buffer = tmap_malloc(sizeof(tmap_seq_t*)*reads_queue_size, "seq_buffer");
  sams = tmap_malloc(sizeof(tmap_map_sams_t*)*reads_queue_size, "sams");

  if(NULL == opt->fn_reads) {
      fp_reads = tmap_file_fdopen(fileno(stdin), "rb", opt->input_compr);
      tmap_progress_set_verbosity(0); 
  }
  else {
      fp_reads = tmap_file_fopen(opt->fn_reads, "rb", opt->input_compr);
  }
  switch(opt->reads_format) {
    case TMAP_READS_FORMAT_FASTA:
    case TMAP_READS_FORMAT_FASTQ:
      seqio = tmap_seq_io_init(fp_reads, TMAP_SEQ_TYPE_FQ);
      for(i=0;i<reads_queue_size;i++) { // initialize the buffer
          seq_buffer[i] = tmap_seq_init(TMAP_SEQ_TYPE_FQ);
      }
      break;
    case TMAP_READS_FORMAT_SFF:
      seqio = tmap_seq_io_init(fp_reads, TMAP_SEQ_TYPE_SFF);
      for(i=0;i<reads_queue_size;i++) { // initialize the buffer
          seq_buffer[i] = tmap_seq_init(TMAP_SEQ_TYPE_SFF);
      }
      break;
    default:
      tmap_error("unrecognized input format", Exit, CommandLineArgument);
      break;
  }

  // Note: 'tmap_file_stdout' should not have been previously modified
  tmap_file_stdout = tmap_file_fdopen(fileno(stdout), "wb", opt->output_compr);

  // SAM header
  tmap_sam_print_header(tmap_file_stdout, refseq, seqio, opt->sam_rg, opt->sam_sff_tags, opt->argc, opt->argv);

  tmap_progress_print("processing reads");
  while(0 < (seq_buffer_length = tmap_seq_io_read_buffer(seqio, seq_buffer, reads_queue_size))) {

      // do alignment
#ifdef HAVE_LIBPTHREAD
      int32_t num_threads = opt->num_threads;
      if(seq_buffer_length < num_threads * TMAP_MAP_ALL_THREAD_BLOCK_SIZE) {
          num_threads = 1 + (seq_buffer_length / TMAP_MAP_ALL_THREAD_BLOCK_SIZE);
      }
      tmap_map_all_read_lock_low = 0; // ALWAYS set before running threads 
      if(1 == num_threads) {
          tmap_map_all_core_worker(seq_buffer, sams, seq_buffer_length, refseq, bwt, sa, 0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          tmap_map_all_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = tmap_calloc(num_threads, sizeof(pthread_t), "threads");
          thread_data = tmap_calloc(num_threads, sizeof(tmap_map_all_thread_data_t), "thread_data");

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
              if(0 != pthread_create(&threads[i], &attr, tmap_map_all_core_thread_worker, &thread_data[i])) {
                  tmap_error("error creating threads", Exit, ThreadError);
              }
          }
          for(i=0;i<num_threads;i++) {
              if(0 != pthread_join(threads[i], NULL)) {
                  tmap_error("error joining threads", Exit, ThreadError);
              }
          }

          free(threads);
          free(thread_data);
      }
#else 
      tmap_map_all_core_worker(seq_buffer, sams, seq_buffer_length, refseq, bwt, sa, 0, opt);
#endif

      if(-1 != opt->reads_queue_size) {
          tmap_progress_print("writing alignments");
      }
      for(i=0;i<seq_buffer_length;i++) {
          // write
          tmap_map_sams_print(seq_buffer[i], refseq, sams[i], opt->sam_sff_tags);

          // free alignments
          tmap_map_sams_destroy(sams[i]);
          sams[i] = NULL;
      }

      if(-1 == opt->reads_queue_size) {
          tmap_file_fflush(tmap_file_stdout, 1);
      }

      n_reads_processed += seq_buffer_length;
      if(-1 != opt->reads_queue_size) {
          tmap_progress_print2("processed %d reads", n_reads_processed);
      }
  }
  if(-1 == opt->reads_queue_size) {
      tmap_progress_print2("processed %d reads", n_reads_processed);
  }

  // close the input/output
  tmap_file_fclose(tmap_file_stdout);
  tmap_file_fclose(fp_reads);

  // free memory
  for(i=0;i<reads_queue_size;i++) {
      tmap_seq_destroy(seq_buffer[i]);
  }
  free(seq_buffer);
  free(sams);
  tmap_refseq_destroy(refseq);
  tmap_bwt_destroy(bwt[0]);
  tmap_bwt_destroy(bwt[1]);
  tmap_sa_destroy(sa[0]);
  tmap_sa_destroy(sa[1]);
  tmap_seq_io_destroy(seqio);
  if(0 < opt->shm_key) {
      tmap_shm_destroy(shm, 0);
  }
}

// for map1/map2/map3
#define __tmap_map_all_opts_copy1(opt_map_all, opt_map_other) do { \
    (opt_map_other)->fn_fasta = tmap_strdup((opt_map_all)->fn_fasta); \
    (opt_map_other)->fn_reads = tmap_strdup((opt_map_all)->fn_reads); \
    (opt_map_other)->reads_format = (opt_map_all)->reads_format; \
    (opt_map_other)->score_match = (opt_map_all)->score_match; \
    (opt_map_other)->pen_mm = (opt_map_all)->pen_mm; \
    (opt_map_other)->pen_gapo = (opt_map_all)->pen_gapo; \
    (opt_map_other)->pen_gape = (opt_map_all)->pen_gape; \
    (opt_map_other)->fscore = (opt_map_all)->fscore; \
    (opt_map_other)->bw = (opt_map_all)->bw; \
    (opt_map_other)->aln_global = (opt_map_all)->aln_global; \
    (opt_map_other)->dup_window = -1; \
    (opt_map_other)->reads_queue_size = (opt_map_all)->reads_queue_size; \
    (opt_map_other)->num_threads = (opt_map_all)->num_threads; \
    (opt_map_other)->aln_output_mode = TMAP_MAP_UTIL_ALN_MODE_ALL; \
    (opt_map_other)->sam_rg = tmap_strdup((opt_map_all)->sam_rg); \
    (opt_map_other)->input_compr = (opt_map_all)->input_compr; \
    (opt_map_other)->output_compr = (opt_map_all)->output_compr; \
    (opt_map_other)->shm_key = (opt_map_all)->shm_key; \
} while(0)

int32_t
tmap_map_all_opt_parse(int argc, char *argv[], tmap_map_opt_t *opt)
{
  int32_t i, j, start, opt_type, opt_type_next, opt_stage, opt_stage_next;

  opt->argc = argc; opt->argv = argv;

  // parse common options as well as map1/map2/map3 commands
  start = 0;
  i = 1;
  opt_type = opt_type_next = TMAP_MAP_ALGO_NONE;
  opt_stage = opt_stage_next = 0;
  opt->num_stages = 0;
  while(i<argc) {
      if(0 == strcmp("map1", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP1;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
          if(0 == opt->num_stages) opt->num_stages++;
      }
      else if(0 == strcmp("map2", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP2;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
          if(0 == opt->num_stages) opt->num_stages++;
      }
      else if(0 == strcmp("map3", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP3;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
          if(0 == opt->num_stages) opt->num_stages++;
      }
      else if(0 == strcmp("MAP1", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP1;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
          if(1 == opt->num_stages) opt->num_stages++;
      }
      else if(0 == strcmp("MAP2", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP2;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
          if(1 == opt->num_stages) opt->num_stages++;
      }
      else if(0 == strcmp("MAP3", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP3;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
          if(1 == opt->num_stages) opt->num_stages++;
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
            case TMAP_MAP_ALGO_NONE:
              // parse common options
              if(0 == tmap_map_opt_parse(i-start, argv+start, opt)) {
                  return 0;
              }
              // copy over common values into the other opts
              for(j=0;j<2;j++) {
                  __tmap_map_all_opts_copy1(opt, opt->opt_map1[j]);
                  __tmap_map_all_opts_copy1(opt, opt->opt_map2[j]);
                  __tmap_map_all_opts_copy1(opt, opt->opt_map3[j]);
              }
              break;
            case TMAP_MAP_ALGO_MAP1:
              if(0 < i - start) {
                  if(0 == tmap_map_opt_parse(i-start, argv+start, opt->opt_map1[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            case TMAP_MAP_ALGO_MAP2:
              // parse map2 options
              if(0 < i - start) {
                  if(0 == tmap_map_opt_parse(i-start, argv+start, opt->opt_map2[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            case TMAP_MAP_ALGO_MAP3:
              // parse map3 options
              if(0 < i - start) {
                  if(0 == tmap_map_opt_parse(i-start, argv+start, opt->opt_map3[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            default:
              tmap_error("bug encountered", Exit, OutOfRange);
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
tmap_map_all_main(int argc, char *argv[])
{
  tmap_map_opt_t *opt = NULL;

  // random seed
  srand48(0); 

  // init opt
  opt = tmap_map_opt_init(TMAP_MAP_ALGO_MAPALL);
      
  // get options
  if(1 != tmap_map_all_opt_parse(argc, argv, opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return tmap_map_opt_usage(opt);
  }
  else { 
      // check command line arguments
      tmap_map_opt_check(opt);
  }

  // run map_all
  tmap_map_all_core(opt);

  // destroy opt
  tmap_map_opt_destroy(opt);

  tmap_progress_print2("terminating successfully");

  return 0;
}
