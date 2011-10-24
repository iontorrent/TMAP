/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <config.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif
#include <unistd.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sam_print.h"
#include "../util/tmap_sort.h"
#include "../util/tmap_rand.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt_gen.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_bwt_match.h"
#include "../index/tmap_sa.h"
#include "../index/tmap_index.h"
#include "../io/tmap_seq_io.h"
#include "../server/tmap_shm.h"
#include "../sw/tmap_fsw.h"
#include "../sw/tmap_sw.h"
#include "tmap_map_stats.h"
#include "tmap_map_util.h"
#include "tmap_map1.h"
#include "tmap_map1_aux.h"
#include "tmap_map2.h"
#include "tmap_map2_aux.h"
#include "tmap_map3.h"
#include "tmap_map3_aux.h"
#include "tmap_map_driver.h"

#define __tmap_map_sam_sort_score_lt(a, b) ((a).score > (b).score)
TMAP_SORT_INIT(tmap_map_sam_sort_score, tmap_map_sam_t, __tmap_map_sam_sort_score_lt)

#define __tmap_map_driver_check_func(_func_init, _func_thread_init, _func_thread_map, _func_thread_cleanup, _func_mapq, _opt) do { \
  if(NULL == _func_init) { \
      tmap_error("func_init == NULL", Exit, OutOfRange); \
  } \
  if(NULL == _func_thread_init) { \
      tmap_error("func_thread_init == NULL", Exit, OutOfRange); \
  } \
  if(NULL == _func_thread_map) { \
      tmap_error("func_thread_map == NULL", Exit, OutOfRange); \
  } \
  if(NULL == _func_thread_cleanup) { \
      tmap_error("func_thread_cleanup == NULL", Exit, OutOfRange); \
  } \
  if(NULL == _func_mapq && _opt->algo_id != TMAP_MAP_ALGO_MAPALL) { \
      tmap_error("func_mapq == NULL", Exit, OutOfRange); \
  } \
} while(0)

static tmap_fsw_flowseq_t *
tmap_map_driver_get_flow_info(tmap_seq_t *seq, tmap_map_opt_t *opt, 
                              uint8_t **flow_order, int32_t *flow_order_len, 
                              uint8_t **key_seq, int32_t *key_seq_len)
{
  int32_t i;
  static int32_t warned = 0;

  // NB: SAM/BAM could be supported, but is not yet
  if(1 == opt->flow_order_use_file) {
      (*flow_order_len) = tmap_seq_get_flow_order_int(seq, flow_order); 
      if(0 == (*flow_order_len)) {
          tmap_error("could not retrieve the flow order from the input file", Exit, OutOfRange);
      }
  }
  else if(NULL != opt->flow_order) {
      (*flow_order_len) = strlen(opt->flow_order);
      (*flow_order) = tmap_malloc(sizeof(uint8_t) * (*flow_order_len), "flow_order");
      for(i=0;i<(*flow_order_len);i++) {
          (*flow_order)[i] = tmap_nt_char_to_int[(int)opt->flow_order[i]];
      }
  }
  else {
      (*flow_order) = NULL;
      (*flow_order_len) = 0;
  }
  if(1 == opt->key_seq_use_file) {
      (*key_seq_len) = tmap_seq_get_key_seq_int(seq, key_seq);
      if(0 == (*key_seq_len)) {
          tmap_error("could not retrieve the key sequence from the input file", Exit, OutOfRange);
      }
  }
  else if(NULL != opt->key_seq) {
      (*key_seq_len) = strlen(opt->key_seq);
      (*key_seq) = tmap_malloc(sizeof(uint8_t) * (*key_seq_len), "key_seq");
      for(i=0;i<(*key_seq_len);i++) {
          (*key_seq)[i] = tmap_nt_char_to_int[(int)opt->key_seq[i]];
      }
  }
  else {
      (*key_seq) = NULL;
      (*key_seq_len) = 0;
  }

  if(0 < (*flow_order_len) && 0 < (*key_seq_len)) {
      return tmap_fsw_flowseq_init(NULL, 0, NULL, NULL, 0, 0, 0);
  }
  else if(0 != (*flow_order_len) && 0 == (*key_seq_len)) {
      if(0 == warned) {
          tmap_error("the flow order was specified but not the key sequence", Warn, OutOfRange);
      }
      warned++;
      return tmap_fsw_flowseq_init(NULL, 0, NULL, NULL, 0, 0, 0);
  }
  else if(0 == (*flow_order_len) && 0 != (*key_seq_len)) {
      tmap_error("the key sequence was specified but not the flow order", Exit, OutOfRange);
  }

  return NULL;
  
}

void
tmap_map_driver_core_worker(int32_t num_ends, tmap_seq_t ***seq_buffer, tmap_map_sams_t ***sams, int32_t seq_buffer_length,
                         tmap_index_t *index,
                         tmap_map_driver_func_thread_init func_thread_init, 
                         tmap_map_driver_func_thread_map func_thread_map, 
                         tmap_map_driver_func_mapq func_mapq,
                         tmap_map_driver_func_thread_cleanup func_thread_cleanup,
                         tmap_map_stats_t *stat,
                         tmap_rand_t *rand,
                         int32_t tid, tmap_map_opt_t *opt)
{
  int32_t i, low = 0;
  void *data = NULL;
  int32_t flow_order_len = 0, key_seq_len = 0;
  uint8_t *flow_order = NULL, *key_seq = NULL;
  tmap_fsw_flowseq_t *fs = NULL;

  // initialize thread data
  if(0 != func_thread_init(&data, opt)) {
      tmap_error("the thread function could not be initialized", Exit, OutOfRange);
  }

  if(0 < seq_buffer_length) {
      fs = tmap_map_driver_get_flow_info(seq_buffer[0][0], opt, &flow_order, &flow_order_len, &key_seq, &key_seq_len);
  }

  // Go through the buffer
  while(low < seq_buffer_length) {
      if(tid == (low % opt->num_threads)) {
          for(i=0;i<num_ends;i++) {
              tmap_seq_t *seq = NULL;

              seq = seq_buffer[i][low];
              
              // remove key sequence for seeding
              if(0 == tmap_seq_remove_key_sequence(seq, opt->remove_sff_clipping, key_seq, key_seq_len)) {
                  // key sequence did not match
                  continue;
              }
              
              stat->num_reads++;

              // map thread data,
              sams[i][low] = func_thread_map(&data, seq, index, stat, rand, opt);
              if(sams[i][low] == NULL) {
                  tmap_error("the thread function did not return a mapping", Exit, OutOfRange);
              }

              if(0 < sams[i][low]->n) {
                  // mapall should have already done this!
                  if(TMAP_MAP_ALGO_MAPALL != opt->algo_id) {
                      stat->num_after_seeding += sams[i][low]->n;

                      // smith waterman (score only)
                      sams[i][low] = tmap_map_util_sw_gen_score(index->refseq, sams[i][low], seq_buffer[i][low], rand, opt);
                      stat->num_after_scoring += sams[i][low]->n;

                      // remove duplicates
                      tmap_map_util_remove_duplicates(sams[i][low], opt->dup_window, rand);
                      stat->num_after_rmdup += sams[i][low]->n;

                      // mapping quality
                      func_mapq(sams[i][low], tmap_seq_get_bases(seq)->l, opt);

                      // set the number of hits before filtering
                      sams[i][low]->max = sams[i][low]->n;

                      // filter alignments
                      tmap_map_sams_filter(sams[i][low], opt->aln_output_mode, rand);
                      stat->num_after_filter += sams[i][low]->n;

                      // smith waterman - generate cigars
                      sams[i][low] = tmap_map_util_sw_gen_cigar(index->refseq, sams[i][low], seq_buffer[i][low], opt);
                  }

                  // re-align the alignments in flow-space
                  if(NULL != fs) {
                      // TODO: if this is run, we do not need to run
                      // tmap_sw_global_banded_core...
                      // NB: seq_buffer should have its key sequence if 0 <
                      // key_seq_len
                      tmap_map_util_fsw(fs, seq,
                                        flow_order, flow_order_len,
                                        key_seq, key_seq_len,
                                        sams[i][low], index->refseq, 
                                        opt->bw, opt->softclip_type, opt->score_thr,
                                        opt->score_match, opt->pen_mm, opt->pen_gapo,
                                        opt->pen_gape, opt->fscore, 1-opt->ignore_flowgram);
                  }

                  // sort by alignment score
                  if(1 < sams[i][low]->n) {
                      tmap_sort_introsort(tmap_map_sam_sort_score,
                                          sams[i][low]->n, sams[i][low]->sams);
                  }
              }
              if(NULL == seq_buffer[i][low]) {
                  tmap_error("bug encoutereed", Exit, OutOfRange);
              }
          }
      }
      // next
      low++;
  }
                  
  // free thread variables
  if(NULL != fs) {
      tmap_fsw_flowseq_destroy(fs);
  }
  free(flow_order);
  free(key_seq);

  // cleanup
  if(0 != func_thread_cleanup(&data, opt)) {
      tmap_error("the thread function could not be initialized", Exit, OutOfRange);
  }
}

void *
tmap_map_driver_core_thread_worker(void *arg)
{
  tmap_map_driver_thread_data_t *thread_data = (tmap_map_driver_thread_data_t*)arg;

  tmap_map_driver_core_worker(thread_data->num_ends, thread_data->seq_buffer, thread_data->sams, thread_data->seq_buffer_length, 
                           thread_data->index,
                           thread_data->func_thread_init, 
                           thread_data->func_thread_map, 
                           thread_data->func_mapq, 
                           thread_data->func_thread_cleanup,
                           thread_data->stat,
                           thread_data->rand,
                           thread_data->tid, thread_data->opt);

  return arg;
}

void 
tmap_map_driver_core(tmap_map_driver_func_init func_init,
                  tmap_map_driver_func_thread_init func_thread_init, 
                  tmap_map_driver_func_thread_map func_thread_map, 
                  tmap_map_driver_func_mapq func_mapq,
                  tmap_map_driver_func_thread_cleanup func_thread_cleanup,
                  tmap_map_opt_t *opt)
{
  uint32_t i, j, n_reads_processed=0;
  int32_t seq_buffer_length=0;
  tmap_seq_io_t **seqio=NULL;
  tmap_seq_t ***seq_buffer = NULL;
  tmap_map_sams_t ***sams = NULL;
  tmap_index_t *index = NULL;
  tmap_map_stats_t *stat = NULL;
#ifdef HAVE_LIBPTHREAD
  tmap_rand_t **rand = NULL;
  tmap_map_stats_t **stats = NULL;
#else
  tmap_rand_t *rand = NULL;
#endif
  int32_t seq_type, reads_queue_size, num_ends;

  // check input functions
  __tmap_map_driver_check_func(func_init, func_thread_init, func_thread_map, func_thread_cleanup, func_mapq, opt);
  
  if(NULL == opt->fn_reads) {
      tmap_progress_set_verbosity(0); 
  }
  
  // open the reads file for reading
  seq_type = tmap_reads_format_to_seq_type(opt->reads_format); 
  num_ends = (0 == opt->fn_reads_num) ? 0 : opt->fn_reads_num;
  seqio = tmap_malloc(sizeof(tmap_seq_io_t*)*num_ends, "seqio");
  for(i=0;i<num_ends;i++) {
      seqio[i] = tmap_seq_io_init(opt->fn_reads[i], seq_type, 0, opt->input_compr);
  }

  // get the index
  index = tmap_index_init(opt->fn_fasta, opt->shm_key);

  // initialize the options and print any relevant information
  if(0 != func_init(index->refseq, opt)) {
      tmap_error("the main function could not be initialized", Exit, OutOfRange);
  }

  // allocate the buffer
  if(-1 == opt->reads_queue_size) {
      reads_queue_size = 1;
  }
  else {
      reads_queue_size = opt->reads_queue_size;
  }
  seq_buffer = tmap_malloc(sizeof(tmap_seq_t**)*num_ends, "seq_buffer");
  sams = tmap_malloc(sizeof(tmap_map_sams_t**)*num_ends, "sams");
  for(i=0;i<num_ends;i++) {
      seq_buffer[i] = tmap_malloc(sizeof(tmap_seq_t*)*reads_queue_size, "seq_buffer[i]");
      sams[i] = tmap_malloc(sizeof(tmap_map_sams_t*)*reads_queue_size, "sams[i]");
      for(j=0;j<reads_queue_size;j++) { // initialize the buffer
          seq_buffer[i][j] = tmap_seq_init(seq_type);
      }
  }

  stat = tmap_map_stats_init();
#ifdef HAVE_LIBPTHREAD
  stats = tmap_malloc(opt->num_threads * sizeof(tmap_map_stats_t*), "stats");
  rand = tmap_malloc(opt->num_threads * sizeof(tmap_rand_t*), "rand");
  for(i=0;i<opt->num_threads;i++) {
      stats[i] = tmap_map_stats_init();
      rand[i] = tmap_rand_init(i);
  }
#else
  rand = tmap_rand_init(13);
#endif

  // Note: 'tmap_file_stdout' should not have been previously modified
  if(NULL == opt->fn_sam) {
      tmap_file_stdout = tmap_file_fdopen(fileno(stdout), "wb", opt->output_compr);
  }
  else {
      tmap_file_stdout = tmap_file_fopen(opt->fn_sam, "wb", opt->output_compr);
  }

  // SAM header
  tmap_sam_print_header(tmap_file_stdout, index->refseq, (1 == num_ends) ? seqio[0] : NULL, 
                        opt->sam_rg, opt->flow_order, opt->key_seq, opt->sam_sff_tags, opt->argc, opt->argv);

  tmap_progress_print("processing reads");
  while(1) {
      // get the reads
      seq_buffer_length = tmap_seq_io_read_buffer(seqio[0], seq_buffer[0], reads_queue_size);
      for(i=1;i<num_ends;i++) {
          if(seq_buffer_length != tmap_seq_io_read_buffer(seqio[i], seq_buffer[i], reads_queue_size)) {
              tmap_error("the input read files were of differing length", Exit, OutOfRange);
          }
      }
      if(0 == seq_buffer_length) {
          break;
      }

      // do alignment
#ifdef HAVE_LIBPTHREAD
      if(1 == opt->num_threads) {
          tmap_map_driver_core_worker(num_ends, seq_buffer, sams, seq_buffer_length, index,
                                   func_thread_init, func_thread_map, func_mapq, func_thread_cleanup, 
                                   stat,
                                   rand[0],
                                   0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          tmap_map_driver_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = tmap_calloc(opt->num_threads, sizeof(pthread_t), "threads");
          thread_data = tmap_calloc(opt->num_threads, sizeof(tmap_map_driver_thread_data_t), "thread_data");

          for(i=0;i<opt->num_threads;i++) {
              thread_data[i].num_ends = num_ends;
              thread_data[i].seq_buffer = seq_buffer;
              thread_data[i].sams = sams;
              thread_data[i].seq_buffer_length = seq_buffer_length;
              thread_data[i].index = index;
              thread_data[i].func_thread_init = func_thread_init;
              thread_data[i].func_thread_map = func_thread_map;
              thread_data[i].func_mapq = func_mapq;
              thread_data[i].func_thread_cleanup = func_thread_cleanup;
              thread_data[i].tid = i;
              thread_data[i].stat = stats[i];
              thread_data[i].rand = rand[i];
              thread_data[i].opt = opt; 
              if(0 != pthread_create(&threads[i], &attr, tmap_map_driver_core_thread_worker, &thread_data[i])) {
                  tmap_error("error creating threads", Exit, ThreadError);
              }
          }
          for(i=0;i<opt->num_threads;i++) {
              if(0 != pthread_join(threads[i], NULL)) {
                  tmap_error("error joining threads", Exit, ThreadError);
              }
              tmap_map_stats_add(stat, stats[i]);
          }

          free(threads);
          free(thread_data);
      }
#else 
      tmap_map_driver_core_worker(num_ends, seq_buffer, sams, seq_buffer_length, index,
                                  func_thread_init, func_thread_map, func_mapq, func_thread_cleanup, 
                                  stat,
                                  rand,
                                  0, opt);
#endif

      if(-1 != opt->reads_queue_size) {
          tmap_progress_print("writing alignments");
      }
      for(i=0;i<seq_buffer_length;i++) {
          // write
          if(1 == num_ends) {
              tmap_map_sams_print(seq_buffer[0][i], index->refseq, sams[0][i], 
                                  0, NULL, opt->sam_sff_tags);
          }
          else {
              for(j=0;j<num_ends;j++) {
                  tmap_map_sams_print(seq_buffer[j][i], index->refseq, sams[j][i], 
                                    (0 == j) ? 1 : ((num_ends-1 == j) ? 2 : 0),
                                    sams[(j+1) % num_ends][i], opt->sam_sff_tags);
              }
          }

          // free alignments
          for(j=0;j<num_ends;j++) {
              tmap_map_sams_destroy(sams[j][i]);
              sams[j][i] = NULL;
          }
      }
      if(-1 == opt->reads_queue_size) {
          tmap_file_fflush(tmap_file_stdout, 1);
      }
      else {
          tmap_file_fflush(tmap_file_stdout, 0); // flush
      }

      n_reads_processed += seq_buffer_length;
      if(-1 != opt->reads_queue_size) {
          tmap_progress_print2("processed %d reads", n_reads_processed);
          tmap_progress_print2("stats [%.2lf,%.2lf,%.2lf,%.2lf]",
                               stat->num_after_seeding/(double)stat->num_reads,
                               stat->num_after_scoring/(double)stat->num_reads,
                               stat->num_after_rmdup/(double)stat->num_reads,
                               stat->num_after_filter/(double)stat->num_reads);
      }
  }
  if(-1 == opt->reads_queue_size) {
      tmap_progress_print2("processed %d reads", n_reads_processed);
      tmap_progress_print2("stats [%.2lf,%.2lf,%.2lf,%.2lf]",
                           stat->num_after_seeding/(double)stat->num_reads,
                           stat->num_after_scoring/(double)stat->num_reads,
                           stat->num_after_rmdup/(double)stat->num_reads,
                           stat->num_after_filter/(double)stat->num_reads);
  }

  // close the input/output
  tmap_file_fclose(tmap_file_stdout);

  // free memory
  tmap_index_destroy(index);
  for(i=0;i<num_ends;i++) {
      tmap_seq_io_destroy(seqio[i]);
      for(j=0;j<reads_queue_size;j++) {
          tmap_seq_destroy(seq_buffer[i][j]);
      }
      free(seq_buffer[i]);
      free(sams[i]);
  }
  free(seqio);
  free(seq_buffer);
  free(sams);
  tmap_map_stats_destroy(stat);
#ifdef HAVE_LIBPTHREAD
  for(i=0;i<opt->num_threads;i++) {
      tmap_map_stats_destroy(stats[i]);
      tmap_rand_destroy(rand[i]);
  }
  free(stats);
  free(rand);
#else
  tmap_rand_destroy(rand);
#endif
}
