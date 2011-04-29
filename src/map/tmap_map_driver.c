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
#include "../util/tmap_sam_print.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt_gen.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_bwt_match.h"
#include "../index/tmap_sa.h"
#include "../io/tmap_seq_io.h"
#include "../server/tmap_shm.h"
#include "../sw/tmap_sw.h"
#include "tmap_map_util.h"
#include "tmap_map1.h"
#include "tmap_map1_aux.h"
#include "tmap_map2.h"
#include "tmap_map2_aux.h"
#include "tmap_map3.h"
#include "tmap_map3_aux.h"
#include "tmap_map_driver.h"

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

void
tmap_map_driver_core_worker(tmap_seq_t **seq_buffer, tmap_map_sams_t **sams, int32_t seq_buffer_length, 
                         tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2],
                         tmap_driver_func_thread_init func_thread_init, 
                         tmap_driver_func_thread_map func_thread_map, 
                         tmap_driver_func_mapq func_mapq,
                         tmap_driver_func_thread_cleanup func_thread_cleanup,
                         int32_t tid, tmap_map_opt_t *opt)
{
  int32_t i, low = 0;
  void *data = NULL;
  int32_t flow_order_len = 0;
  uint8_t *flow_order = NULL;

  // initialize thread data
  if(0 != func_thread_init(&data, opt)) {
      tmap_error("the thread function could not be initialized", Exit, OutOfRange);
  }
              
  if(1 == opt->flow_order_use_sff) { // initialize the flow order from the SFF header
      flow_order_len = seq_buffer[low]->data.sff->gheader->flow->l;
      flow_order = tmap_malloc(sizeof(uint8_t) * flow_order_len, "flow_order");
      for(i=0;i<flow_order_len;i++) {
          flow_order[i] = tmap_nt_char_to_int[(int)seq_buffer[low]->data.sff->gheader->flow->s[i]];
      }
  }
  else if(NULL != opt->flow_order) {
      flow_order_len = strlen(opt->flow_order);
      flow_order = tmap_malloc(sizeof(uint8_t) * flow_order_len, "flow_order");
      for(i=0;i<flow_order_len;i++) {
          flow_order[i] = tmap_nt_char_to_int[(int)opt->flow_order[i]];
      }
  }

  while(low < seq_buffer_length) {
      if(tid == (low % opt->num_threads)) {
          tmap_seq_t *seq = NULL;

          // Note: for fsw re-alignment, we need the key bases, so do not get
          // rid of them just yet
          if(0 < flow_order_len) {
              seq = tmap_seq_clone(seq_buffer[low]);
          }
          else {
              seq = seq_buffer[low];
          }

          // remove key sequence for seeding
          tmap_seq_remove_key_sequence(seq);

          // map thread data,
          sams[low] = func_thread_map(&data, seq, refseq, bwt, sa, opt);
          if(sams[low] == NULL) {
              tmap_error("the thread function did not return a mapping", Exit, OutOfRange);
          }

          if(0 < sams[low]->n) {
              // mapall should have already done this!
              if(TMAP_MAP_ALGO_MAPALL != opt->algo_id) {
                  // TODO: we could intelligently remove duplicates before smith
                  // waterman, or hash the smith waterman results!
                  
                  // smith waterman
                  sams[low] = tmap_map_util_sw(refseq, sams[low], seq_buffer[low], opt);

                  // remove duplicates
                  tmap_map_util_remove_duplicates(sams[low], opt->dup_window);

                  // mapping quality
                  func_mapq(sams[low], tmap_seq_get_bases(seq)->l, opt);

                  // filter alignments
                  tmap_map_sams_filter(sams[low], opt->aln_output_mode);
              }

              // re-align the alignments in flow-space
              if(0 < flow_order_len) {
                  // Note: seq_buffer should have its key sequence
                  tmap_map_util_fsw(seq_buffer[low],
                                    flow_order, flow_order_len,
                                    sams[low], refseq, 
                                    opt->bw, opt->softclip_type, opt->score_thr,
                                    opt->score_match, opt->pen_mm, opt->pen_gapo,
                                    opt->pen_gape, opt->fscore);
                  // remove key sequence, do not output the key sequence part
                  tmap_seq_remove_key_sequence(seq_buffer[low]);
              }
          }
      }
      // next
      low++;
  }
                  
  // free thread variables
  free(flow_order);

  // cleanup
  if(0 != func_thread_cleanup(&data, opt)) {
      tmap_error("the thread function could not be initialized", Exit, OutOfRange);
  }
}

void *
tmap_map_driver_core_thread_worker(void *arg)
{
  tmap_map_driver_thread_data_t *thread_data = (tmap_map_driver_thread_data_t*)arg;

  tmap_map_driver_core_worker(thread_data->seq_buffer, thread_data->sams, thread_data->seq_buffer_length, 
                           thread_data->refseq, thread_data->bwt, thread_data->sa, 
                           thread_data->func_thread_init, 
                           thread_data->func_thread_map, 
                           thread_data->func_mapq, 
                           thread_data->func_thread_cleanup,
                           thread_data->tid, thread_data->opt);

  return arg;
}

void 
tmap_map_driver_core(tmap_driver_func_init func_init,
                  tmap_driver_func_thread_init func_thread_init, 
                  tmap_driver_func_thread_map func_thread_map, 
                  tmap_driver_func_mapq func_mapq,
                  tmap_driver_func_thread_cleanup func_thread_cleanup,
                  tmap_map_opt_t *opt)
{
  uint32_t i, n_reads_processed=0;
  int32_t seq_buffer_length;
  tmap_refseq_t *refseq=NULL;
  tmap_bwt_t *bwt[2]={NULL, NULL};
  tmap_sa_t *sa[2]={NULL, NULL};
  tmap_seq_io_t *seqio = NULL;
  tmap_seq_t **seq_buffer = NULL;
  tmap_map_sams_t **sams = NULL;
  tmap_shm_t *shm = NULL;
  int32_t seq_type, reads_queue_size;

  // check input functions
  __tmap_map_driver_check_func(func_init, func_thread_init, func_thread_map, func_thread_cleanup, func_mapq, opt);
  
  if(NULL == opt->fn_reads) {
      tmap_progress_set_verbosity(0); 
  }
  
  // open the reads file for reading
  seq_type = tmap_reads_format_to_seq_type(opt->reads_format); 
  seqio = tmap_seq_io_init(opt->fn_reads, seq_type, 0, opt->input_compr);

  // get the reference information
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

  // initialize the options and print any relevant information
  if(0 != func_init(refseq, opt)) {
      tmap_error("the main function could not be initialized", Exit, OutOfRange);
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
  for(i=0;i<reads_queue_size;i++) { // initialize the buffer
      seq_buffer[i] = tmap_seq_init(seq_type);
  }

  // Note: 'tmap_file_stdout' should not have been previously modified
  tmap_file_stdout = tmap_file_fdopen(fileno(stdout), "wb", opt->output_compr);

  // SAM header
  tmap_sam_print_header(tmap_file_stdout, refseq, seqio, opt->sam_rg, opt->flow_order, opt->sam_sff_tags, opt->argc, opt->argv);

  tmap_progress_print("processing reads");
  while(0 < (seq_buffer_length = tmap_seq_io_read_buffer(seqio, seq_buffer, reads_queue_size))) {

      // do alignment
#ifdef HAVE_LIBPTHREAD
      if(1 == opt->num_threads) {
          tmap_map_driver_core_worker(seq_buffer, sams, seq_buffer_length, refseq, bwt, sa, 
                                   func_thread_init, func_thread_map, func_mapq, func_thread_cleanup, 
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
              thread_data[i].seq_buffer = seq_buffer;
              thread_data[i].seq_buffer_length = seq_buffer_length;
              thread_data[i].sams = sams;
              thread_data[i].refseq = refseq;
              thread_data[i].bwt[0] = bwt[0];
              thread_data[i].bwt[1] = bwt[1];
              thread_data[i].sa[0] = sa[0];
              thread_data[i].sa[1] = sa[1];
              thread_data[i].func_thread_init = func_thread_init;
              thread_data[i].func_thread_map = func_thread_map;
              thread_data[i].func_mapq = func_mapq;
              thread_data[i].func_thread_cleanup = func_thread_cleanup;
              thread_data[i].tid = i;
              thread_data[i].opt = opt; 
              if(0 != pthread_create(&threads[i], &attr, tmap_map_driver_core_thread_worker, &thread_data[i])) {
                  tmap_error("error creating threads", Exit, ThreadError);
              }
          }
          for(i=0;i<opt->num_threads;i++) {
              if(0 != pthread_join(threads[i], NULL)) {
                  tmap_error("error joining threads", Exit, ThreadError);
              }
          }

          free(threads);
          free(thread_data);
      }
#else 
      tmap_map_driver_core_worker(seq_buffer, sams, seq_buffer_length, refseq, bwt, sa, 
                                  func_thread_init, func_thread_map, func_mapq, func_thread_cleanup, 
                                  0, opt);
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
      else {
          tmap_file_fflush(tmap_file_stdout, 0); // flush
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
