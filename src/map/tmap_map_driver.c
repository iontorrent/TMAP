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
#include "util/tmap_map_stats.h"
#include "util/tmap_map_util.h"
#include "tmap_map_driver.h"

#define __tmap_map_sam_sort_score_lt(a, b) ((a).score > (b).score)
TMAP_SORT_INIT(tmap_map_sam_sort_score, tmap_map_sam_t, __tmap_map_sam_sort_score_lt)

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
tmap_map_driver_core_worker(int32_t num_ends, 
                            tmap_seq_t ***seq_buffer, 
                            tmap_map_record_t **records, 
                            int32_t seq_buffer_length,
                            tmap_index_t *index,
                            tmap_map_driver_t *driver,
                            tmap_map_stats_t *stat,
                            tmap_rand_t *rand,
                            int32_t tid)
{
  int32_t i, j, k, low = 0;
  int32_t flow_order_len = 0, key_seq_len = 0;
  uint8_t *flow_order = NULL, *key_seq = NULL;
  int32_t found;
  tmap_fsw_flowseq_t *fs = NULL;
  tmap_seq_t ***seqs = NULL;
  void **data = NULL; // one per algorithm

  // initialize thread data
  data = tmap_malloc(sizeof(void*) * driver->num_algorithms, "data");
  for(i=0;i<driver->num_algorithms;i++) {
      if(NULL != driver->algorithms[i]->func_thread_init && 0 != driver->algorithms[i]->func_thread_init(&data[i], driver->opt)) {
          tmap_error("the thread function could not be initialized", Exit, OutOfRange);
      }
  }
  // initialize flow space info
  if(0 < seq_buffer_length) {
      fs = tmap_map_driver_get_flow_info(seq_buffer[0][0], driver->opt, &flow_order, &flow_order_len, &key_seq, &key_seq_len);
  }
  // init memory
  seqs = tmap_malloc(sizeof(tmap_seq_t**)*num_ends, "seqs");
  for(i=0;i<num_ends;i++) {
      seqs[i] = tmap_malloc(sizeof(tmap_seq_t*)*4, "seqs[i]");
  }

  // Go through the buffer
  while(low < seq_buffer_length) {
      if(tid == (low % driver->opt->num_threads)) {
          tmap_map_stats_t *curstat = NULL;
          tmap_map_record_t *record_prev = NULL;
              
          // remove key sequences
          for(i=0;i<num_ends;i++) {
              tmap_seq_t *seq = NULL;
              seq = seq_buffer[i][low];
              if(0 == tmap_seq_remove_key_sequence(seq, driver->opt->remove_sff_clipping, key_seq, key_seq_len)) {
                  // key sequence did not match
                  continue;
              }
          }

          // init
          for(i=0;i<num_ends;i++) {
              tmap_seq_t *seq = NULL;
              seq = seq_buffer[i][low];
              // init seqs
              for(j=0;j<4;j++) {
                  // TODO: only if necessary
                  seqs[i][j]= tmap_seq_clone(seq); // clone the sequence 
                  switch(j) {
                    case 0: // forward
                      break;
                    case 1: // reverse compliment
                      tmap_seq_reverse_compliment(seqs[i][j]); break;
                    case 2: // reverse
                      tmap_seq_reverse(seqs[i][j]); break;
                    case 3: // compliment
                      tmap_seq_compliment(seqs[i][j]); break;
                  }
                  tmap_seq_to_int(seqs[i][j]); // convert to integers
              }
          }

          // init records
          records[low] = tmap_map_record_init(num_ends);

          // go through each stage
          for(i=0;i<driver->opt->num_stages;i++) { // for each stage
              curstat = tmap_calloc(1, sizeof(tmap_map_stats_t), "curstat");
              // seed
              found = 0;
              for(j=0;j<num_ends;j++) { // for each end
                  curstat->num_reads++;
                  for(k=0;k<driver->num_algorithms;k++) { // for each algorithm
                      tmap_map_sams_t *sams = NULL;
                      if(i+1 != driver->algorithms[k]->opt->algo_stage) continue; // next algorithm
                      // map
                      sams = driver->algorithms[k]->func_thread_map(&data[k], seqs[i], index, curstat, rand, driver->algorithms[k]->opt);
                      if(NULL == sams) {
                          tmap_error("the thread function did not return a mapping", Exit, OutOfRange);
                      }
                      // append
                      tmap_map_sams_merge(records[low]->sams[j], sams);
                      // destroy
                      tmap_map_sams_destroy(sams);
                  }
                  curstat->num_after_seeding += records[low]->sams[j]->n;
              }

              // keep mappings for subsequent stages or restore mappings from
              // previous stages
              if(1 == driver->opt->mapall_keep_all) {
                  // merge from the previous stage
                  if(0 < i) {
                      tmap_map_record_merge(records[low], record_prev);
                  }

                  // keep for the next stage
                  if(i < driver->opt->num_stages-1) { // more stages left
                      record_prev = tmap_map_record_clone(records[low]);
                  }
              }

              // generate scores with smith waterman
              for(j=0;j<num_ends;j++) { // for each end
                  records[low]->sams[j] = tmap_map_util_sw_gen_score(index->refseq, records[low]->sams[j], seqs[j], rand, driver->opt);
                  curstat->num_after_scoring += records[low]->sams[j]->n;
              }

              // remove duplicates
              for(j=0;j<num_ends;j++) { // for each end
                  tmap_map_util_remove_duplicates(records[low]->sams[j], driver->opt->dup_window, rand);
                  curstat->num_after_rmdup += records[low]->sams[j]->n;
              }
              
              // (single-end) mapping quality
              for(j=0;j<num_ends;j++) { // for each end
                  driver->func_mapq(records[low]->sams[j], tmap_seq_get_bases_length(seqs[j][0]), driver->opt);
              }

              // filter if we have more stages
              if(i < driver->opt->num_stages-1) {
                  for(j=0;j<num_ends;j++) { // for each end
                      tmap_map_sams_filter2(records[low]->sams[j], driver->opt->mapall_score_thr, driver->opt->mapall_mapq_thr);
                  }
              }

              // choose alignments
              for(j=0;j<num_ends;j++) { // for each end
                  tmap_map_sams_filter1(records[low]->sams[j], driver->opt->aln_output_mode, TMAP_MAP_ALGO_NONE, rand);
                  curstat->num_after_filter += records[low]->sams[j]->n;
              }

              // TODO: PAIRING HERE
              //
              // TODO: mapping quality
              //
              // TODO: filtering 
              //
              // TODO: if we have one end for a pair, do we go onto the second
              // stage?

              // generate the cigars
              found = 0;
              for(j=0;j<num_ends;j++) { // for each end
                  records[low]->sams[j] = tmap_map_util_sw_gen_cigar(index->refseq, records[low]->sams[j], seqs[j], driver->opt);
                  if(0 < records[low]->sams[j]->n) {
                      curstat->num_with_mapping++;
                      found = 1;
                  }
              }

              // did we find any mappings?
              // TODO:
              if(1 == found) { // yes
                  tmap_map_stats_add(stat, curstat);
                  tmap_map_stats_destroy(curstat);
                  break;
              }
              else { // no
                  tmap_map_record_destroy(records[low]);
                  // re-init
                  records[low] = tmap_map_record_init(num_ends);
              }
              tmap_map_stats_destroy(curstat);
          }

          // flowspace re-align and sorting
          for(i=0;i<num_ends;i++) {
              if(0 < records[low]->sams[i]->n) {
                  // re-align the alignments in flow-space
                  if(NULL != fs) {
                      // TODO: if this is run, we do not need to run
                      // tmap_sw_global_banded_core...
                      // NB: seq_buffer should have its key sequence if 0 <
                      // key_seq_len
                      tmap_map_util_fsw(fs, seqs[i][0],
                                        flow_order, flow_order_len,
                                        key_seq, key_seq_len,
                                        records[low]->sams[i], index->refseq, 
                                        driver->opt->bw, driver->opt->softclip_type, driver->opt->score_thr,
                                        driver->opt->score_match, driver->opt->pen_mm, driver->opt->pen_gapo,
                                        driver->opt->pen_gape, driver->opt->fscore, 1-driver->opt->ignore_flowgram);
                  }

                  // sort by alignment score
                  if(1 < records[low]->sams[i]->n) {
                      tmap_sort_introsort(tmap_map_sam_sort_score,
                                          records[low]->sams[i]->n, 
                                          records[low]->sams[i]->sams);
                  }
              }
              if(NULL == seq_buffer[i][low]) {
                  tmap_error("bug encoutereed", Exit, OutOfRange);
              }
          }

          // free seqs
          for(i=0;i<num_ends;i++) {
              for(j=0;j<4;j++) {
                  tmap_seq_destroy(seqs[i][j]);
                  seqs[i][j] = NULL;
              }
          }
          // TODO
          // free other things?
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
  for(i=0;i<num_ends;i++) {
      free(seqs[i]);
  }
  free(seqs);

  // cleanup
  for(i=0;i<driver->num_algorithms;i++) {
      if(NULL != driver->algorithms[i]->func_thread_cleanup && 0 != driver->algorithms[i]->func_thread_cleanup(&data[i], driver->opt)) {
          tmap_error("the thread function could not be initialized", Exit, OutOfRange);
      }
  }
  free(data);
}

void *
tmap_map_driver_core_thread_worker(void *arg)
{
  tmap_map_driver_thread_data_t *thread_data = (tmap_map_driver_thread_data_t*)arg;

  tmap_map_driver_core_worker(thread_data->num_ends, thread_data->seq_buffer, thread_data->records, thread_data->seq_buffer_length, 
                           thread_data->index, thread_data->driver, thread_data->stat, thread_data->rand, thread_data->tid);

  return arg;
}

void 
tmap_map_driver_core(tmap_map_driver_t *driver)
{
  uint32_t i, j, n_reads_processed=0;
  int32_t seq_buffer_length=0;
  tmap_seq_io_t **seqio=NULL;
  tmap_seq_t ***seq_buffer = NULL;
  tmap_map_record_t **records=NULL;
  tmap_index_t *index = NULL;
  tmap_map_stats_t *stat = NULL;
#ifdef HAVE_LIBPTHREAD
  tmap_rand_t **rand = NULL;
  tmap_map_stats_t **stats = NULL;
#else
  tmap_rand_t *rand = NULL;
#endif
  int32_t seq_type, reads_queue_size, num_ends;

  // HERE
  //tmap_map_opt_print(driver->opt);

  if(NULL == driver->opt->fn_reads) {
      tmap_progress_set_verbosity(0); 
  }
  
  // open the reads file for reading
  seq_type = tmap_reads_format_to_seq_type(driver->opt->reads_format); 
  num_ends = (0 == driver->opt->fn_reads_num) ? 0 : driver->opt->fn_reads_num;
  seqio = tmap_malloc(sizeof(tmap_seq_io_t*)*num_ends, "seqio");
  for(i=0;i<num_ends;i++) {
      seqio[i] = tmap_seq_io_init(driver->opt->fn_reads[i], seq_type, 0, driver->opt->input_compr);
  }

  // get the index
  index = tmap_index_init(driver->opt->fn_fasta, driver->opt->shm_key);

  // initialize the driver->options and print any relevant information
  for(i=0;i<driver->num_algorithms;i++) {
      if(NULL != driver->algorithms[i]->func_init && 0 != driver->algorithms[i]->func_init(&driver->data[i], index->refseq, driver->opt)) {
          tmap_error("the main function could not be initialized", Exit, OutOfRange);
      }
  }

  // allocate the buffer
  if(-1 == driver->opt->reads_queue_size) {
      reads_queue_size = 1;
  }
  else {
      reads_queue_size = driver->opt->reads_queue_size;
  }
  seq_buffer = tmap_malloc(sizeof(tmap_seq_t**)*num_ends, "seq_buffer");
  for(i=0;i<num_ends;i++) {
      seq_buffer[i] = tmap_malloc(sizeof(tmap_seq_t*)*reads_queue_size, "seq_buffer[i]");
      for(j=0;j<reads_queue_size;j++) { // initialize the buffer
          seq_buffer[i][j] = tmap_seq_init(seq_type);
      }
  }
  records = tmap_malloc(sizeof(tmap_map_record_t*)*reads_queue_size, "records");

  stat = tmap_map_stats_init();
#ifdef HAVE_LIBPTHREAD
  stats = tmap_malloc(driver->opt->num_threads * sizeof(tmap_map_stats_t*), "stats");
  rand = tmap_malloc(driver->opt->num_threads * sizeof(tmap_rand_t*), "rand");
  for(i=0;i<driver->opt->num_threads;i++) {
      stats[i] = tmap_map_stats_init();
      rand[i] = tmap_rand_init(i);
  }
#else
  rand = tmap_rand_init(13);
#endif

  // Note: 'tmap_file_stdout' should not have been previously modified
  if(NULL == driver->opt->fn_sam) {
      tmap_file_stdout = tmap_file_fdopen(fileno(stdout), "wb", driver->opt->output_compr);
  }
  else {
      tmap_file_stdout = tmap_file_fopen(driver->opt->fn_sam, "wb", driver->opt->output_compr);
  }

  // SAM header
  tmap_sam_print_header(tmap_file_stdout, index->refseq, (1 == num_ends) ? seqio[0] : NULL, 
                        driver->opt->sam_rg, driver->opt->flow_order, driver->opt->key_seq, driver->opt->sam_sff_tags, driver->opt->argc, driver->opt->argv);

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
      if(1 == driver->opt->num_threads) {
          tmap_map_driver_core_worker(num_ends, seq_buffer, records, seq_buffer_length, index,
                                      driver, stat, rand[0], 0);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          tmap_map_driver_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = tmap_calloc(driver->opt->num_threads, sizeof(pthread_t), "threads");
          thread_data = tmap_calloc(driver->opt->num_threads, sizeof(tmap_map_driver_thread_data_t), "thread_data");

          // create threads
          for(i=0;i<driver->opt->num_threads;i++) {
              thread_data[i].num_ends = num_ends;
              thread_data[i].seq_buffer = seq_buffer;
              thread_data[i].seq_buffer_length = seq_buffer_length;
              thread_data[i].records = records;
              thread_data[i].index = index;
              thread_data[i].driver = driver;
              thread_data[i].stat = stats[i];
              thread_data[i].rand = rand[i];
              thread_data[i].tid = i;
              if(0 != pthread_create(&threads[i], &attr, tmap_map_driver_core_thread_worker, &thread_data[i])) {
                  tmap_error("error creating threads", Exit, ThreadError);
              }
          }

          // join threads
          for(i=0;i<driver->opt->num_threads;i++) {
              if(0 != pthread_join(threads[i], NULL)) {
                  tmap_error("error joining threads", Exit, ThreadError);
              }
              // add the stats
              tmap_map_stats_add(stat, stats[i]);
          }

          free(threads);
          free(thread_data);
      }
#else 
      tmap_map_driver_core_worker(num_ends, seq_buffer, records, seq_buffer_length, index,
                                  driver, stat, rand, 0);
#endif

      if(-1 != driver->opt->reads_queue_size) {
          tmap_progress_print("writing alignments");
      }
      for(i=0;i<seq_buffer_length;i++) {
          // write
          if(1 == num_ends) {
              tmap_map_sams_print(seq_buffer[0][i], index->refseq, records[i]->sams[0], 
                                  0, NULL, driver->opt->sam_sff_tags);
          }
          else {
              for(j=0;j<num_ends;j++) {
                  tmap_map_sams_print(seq_buffer[j][i], index->refseq, records[i]->sams[j],
                                    (0 == j) ? 1 : ((num_ends-1 == j) ? 2 : 0),
                                    records[i]->sams[(j+1) % num_ends], 
                                    driver->opt->sam_sff_tags);
              }
          }

          // free alignments
          tmap_map_record_destroy(records[i]); // TODO: needs implementing
          records[i] = NULL;
      }
      if(-1 == driver->opt->reads_queue_size) {
          tmap_file_fflush(tmap_file_stdout, 1);
      }
      else {
          tmap_file_fflush(tmap_file_stdout, 0); // flush
      }

      n_reads_processed += seq_buffer_length;
      if(-1 != driver->opt->reads_queue_size) {
          tmap_progress_print2("processed %d reads", n_reads_processed);
          tmap_progress_print2("stats [%.2lf,%.2lf,%.2lf,%.2lf,%.2lf]",
                               stat->num_with_mapping * 100.0 / (double)stat->num_reads,
                               stat->num_after_seeding/(double)stat->num_with_mapping,
                               stat->num_after_scoring/(double)stat->num_with_mapping,
                               stat->num_after_rmdup/(double)stat->num_with_mapping,
                               stat->num_after_filter/(double)stat->num_with_mapping);
      }
  }
  if(-1 == driver->opt->reads_queue_size) {
      tmap_progress_print2("processed %d reads", n_reads_processed);
      tmap_progress_print2("stats [%.2lf,%.2lf,%.2lf,%.2lf,%.2lf]",
                           stat->num_with_mapping * 100.0 / (double)stat->num_reads,
                           stat->num_after_seeding/(double)stat->num_with_mapping,
                           stat->num_after_scoring/(double)stat->num_with_mapping,
                           stat->num_after_rmdup/(double)stat->num_with_mapping,
                           stat->num_after_filter/(double)stat->num_with_mapping);
  }

  // cleanup the algorithm persistent data
  for(i=0;i<driver->num_algorithms;i++) {
      if(NULL != driver->algorithms[i]->func_cleanup && 0 != driver->algorithms[i]->func_cleanup(&driver->data[i])) {
          tmap_error("the main function could not be cleaned up", Exit, OutOfRange);
      }
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
  }
  free(seqio);
  free(seq_buffer);
  free(records);
  tmap_map_stats_destroy(stat);
#ifdef HAVE_LIBPTHREAD
  for(i=0;i<driver->opt->num_threads;i++) {
      tmap_map_stats_destroy(stats[i]);
      tmap_rand_destroy(rand[i]);
  }
  free(stats);
  free(rand);
#else
  tmap_rand_destroy(rand);
#endif
}

/* MAIN API */

tmap_map_driver_algorithm_t*
tmap_map_driver_algorithm_init(tmap_map_driver_func_init func_init,
                    tmap_map_driver_func_thread_init func_thread_init,
                    tmap_map_driver_func_thread_map func_thread_map,
                    tmap_map_driver_func_thread_cleanup func_thread_cleanup,
                    tmap_map_driver_func_cleanup func_cleanup,
                    tmap_map_opt_t *opt)
{
  tmap_map_driver_algorithm_t *algorithm = NULL;
  algorithm = tmap_calloc(1, sizeof(tmap_map_driver_algorithm_t), "algorithm");
  algorithm->func_init = func_init;
  algorithm->func_thread_init = func_thread_init;
  algorithm->func_thread_map = func_thread_map;
  algorithm->func_thread_cleanup = func_thread_cleanup;
  algorithm->func_cleanup = func_cleanup;
  algorithm->opt = opt;
  return algorithm;
}

void
tmap_map_driver_algorithm_destroy(tmap_map_driver_algorithm_t *algorithm)
{
  free(algorithm);
}

tmap_map_driver_t*
tmap_map_driver_init(int32_t algo_id, tmap_map_driver_func_mapq func_mapq)
{
  tmap_map_driver_t *driver = NULL;
  driver = tmap_calloc(1, sizeof(tmap_map_driver_t), "driver");
  driver->opt = tmap_map_opt_init(algo_id);
  driver->func_mapq = func_mapq;
  return driver;
}

void
tmap_map_driver_add(tmap_map_driver_t *driver,
                    tmap_map_driver_func_init func_init,
                    tmap_map_driver_func_thread_init func_thread_init,
                    tmap_map_driver_func_thread_map func_thread_map,
                    tmap_map_driver_func_thread_cleanup func_thread_cleanup,
                    tmap_map_driver_func_cleanup func_cleanup,
                    tmap_map_opt_t *opt)
{
  // TODO: check input functions
  // TODO: check that the options match the global options
  // add a new algorithm
  driver->algorithms = tmap_realloc(driver->algorithms, sizeof(tmap_map_driver_algorithm_t*) * (1+driver->num_algorithms), "driver->algorithms");
  driver->algorithms[driver->num_algorithms] = tmap_map_driver_algorithm_init(func_init, func_thread_init, func_thread_map,
                                                                              func_thread_cleanup, func_cleanup, opt);
  driver->data = tmap_realloc(driver->data, sizeof(void*) * driver->num_algorithms, "driver->data");
  driver->num_algorithms++;
  if(driver->opt->num_stages < opt->algo_stage) {
      driver->opt->num_stages = opt->algo_stage;
  }
}

void
tmap_map_driver_run(tmap_map_driver_t *driver)
{
  tmap_map_driver_core(driver);
}

void
tmap_map_driver_destroy(tmap_map_driver_t *driver)
{
  int32_t i;
  for(i=0;i<driver->num_algorithms;i++) {
      tmap_map_driver_algorithm_destroy(driver->algorithms[i]);
  }
  tmap_map_opt_destroy(driver->opt);
  free(driver->data);
  free(driver->algorithms);
  free(driver);
}
