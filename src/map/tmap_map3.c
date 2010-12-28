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
#include "tmap_map3_aux.h"
#include "tmap_map3.h"

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t tmap_map3_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t tmap_map3_read_lock_low = 0;
#define TMAP_MAP3_THREAD_BLOCK_SIZE 512
#endif

int32_t
tmap_map3_get_seed_length(uint64_t ref_len)
{
  int32_t k = 0;
  while(0 < ref_len) {
      ref_len >>= 2; // divide by two
      k++;
  }
  return k;
}

static inline void
tmap_map3_mapq(tmap_map_sams_t *sams, int32_t score_thr, int32_t score_match, int32_t aln_output_mode)
{
  int32_t i;
  int32_t n_best = 0;
  int32_t best_score, cur_score, best_subo;
  int32_t n_seeds = 0, tot_seeds = 0;
  int32_t mapq;

  // estimate mapping quality TODO: this needs to be refined
  best_score = INT32_MIN;
  best_subo = INT32_MIN;
  n_best = 0;
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      tot_seeds += sams->sams[i].aux.map3_aux->n_seeds;
      if(best_score < cur_score) {
          best_subo = best_score;
          best_score = cur_score;
          n_best = 1;
          n_seeds = sams->sams[i].aux.map3_aux->n_seeds;
      }
      else if(cur_score == best_score) { // qual
          n_best++;
      }
      else {
          if(best_subo < cur_score) {
              best_subo = cur_score;
          }
          cur_score = sams->sams[i].score_subo;
          if(best_subo < cur_score) {
              best_subo = cur_score;
          } 
      }
  }
  if(1 < n_best) {
      mapq = 0;
  }
  else {
      double c = 0.0;
      c = n_seeds / (double)tot_seeds;
      if(best_subo < score_thr) best_subo = score_thr;
      mapq = (int32_t)(c * (best_score - best_subo) * (250.0 / best_score + 0.03 / score_match) + .499);
      if(mapq > 250) mapq = 250;
  }
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      if(cur_score == best_score) { // qual
          sams->sams[i].mapq = mapq;
      }
      else {
          sams->sams[i].mapq = 0;
      }
  }
}

static void
tmap_map3_core_worker(tmap_seq_t **seq_buffer, tmap_map_sams_t **sams, int32_t seq_buffer_length, 
                      tmap_refseq_t *refseq, tmap_bwt_t *bwt, tmap_sa_t *sa,
                      int32_t tid, tmap_map_opt_t *opt)
{
  int32_t i, low = 0, high;
  uint8_t *flow[2] = {NULL, NULL};

  // set up the flow order
  if(0 < seq_buffer_length && 0 < opt->hp_diff) {
      if(TMAP_SEQ_TYPE_SFF == seq_buffer[0]->type) {
          flow[0] = tmap_malloc(sizeof(uint8_t) * 4, "flow[0]");
          flow[1] = tmap_malloc(sizeof(uint8_t) * 4, "flow[1]");
          for(i=0;i<4;i++) {
              flow[0][i] = tmap_nt_char_to_int[(int)seq_buffer[0]->data.sff->gheader->flow->s[i]];
              flow[1][i] = 3 - tmap_nt_char_to_int[(int)seq_buffer[0]->data.sff->gheader->flow->s[i]];
          }
      }
      else {
          tmap_error("bug encountered", Exit, OutOfRange);
      }
  }

  while(low < seq_buffer_length) {
#ifdef HAVE_LIBPTHREAD
      if(1 < opt->num_threads) {
          pthread_mutex_lock(&tmap_map3_read_lock);

          // update bounds
          low = tmap_map3_read_lock_low;
          tmap_map3_read_lock_low += TMAP_MAP3_THREAD_BLOCK_SIZE;
          high = low + TMAP_MAP3_THREAD_BLOCK_SIZE;
          if(seq_buffer_length < high) {
              high = seq_buffer_length; 
          }

          pthread_mutex_unlock(&tmap_map3_read_lock);
      }
      else {
          high = seq_buffer_length; // process all
      }
#else 
      high = seq_buffer_length; // process all
#endif
      while(low<high) {
          tmap_seq_t *seq[2]={NULL, NULL}, *orig_seq=NULL;
          orig_seq = seq_buffer[low];

          // clone the sequence 
          seq[0] = tmap_seq_clone(seq_buffer[low]);
          seq[1] = tmap_seq_clone(seq_buffer[low]);
          
          // Adjust for SFF
          tmap_seq_remove_key_sequence(seq[0]);
          tmap_seq_remove_key_sequence(seq[1]);

          // Adjust for SFF
          tmap_seq_remove_key_sequence(seq[0]);
          tmap_seq_remove_key_sequence(seq[1]);

          // reverse compliment
          tmap_seq_reverse_compliment(seq[1]);

          // convert to integers
          tmap_seq_to_int(seq[0]);
          tmap_seq_to_int(seq[1]);

          // align
          sams[low] = tmap_map3_aux_core(seq, flow, refseq, bwt, sa, opt);

          // mapping quality
          tmap_map3_mapq(sams[low], opt->score_thr, opt->score_match, opt->aln_output_mode);

          // filter the alignments
          tmap_map_sams_filter(sams[low], opt->aln_output_mode);

          // re-align the alignments in flow-space
          if(TMAP_SEQ_TYPE_SFF == seq_buffer[low]->type) {
              tmap_map_util_fsw(seq_buffer[low]->data.sff, 
                                sams[low], refseq, 
                                opt->bw, opt->aln_global, opt->score_thr,
                                opt->score_match, opt->pen_mm, opt->pen_gapo,
                                opt->pen_gape, opt->fscore);
          }

          // destroy
          tmap_seq_destroy(seq[0]);
          tmap_seq_destroy(seq[1]);

          // next
          low++;
      }
  }

  // free
  free(flow[0]);
  free(flow[1]);
}

static void *
tmap_map3_core_thread_worker(void *arg)
{
  tmap_map3_thread_data_t *thread_data = (tmap_map3_thread_data_t*)arg;

  tmap_map3_core_worker(thread_data->seq_buffer, thread_data->sams, thread_data->seq_buffer_length, 
                        thread_data->refseq, thread_data->bwt, thread_data->sa, 
                        thread_data->tid, thread_data->opt);

  return arg;
}

static void 
tmap_map3_core(tmap_map_opt_t *opt)
{
  uint32_t i, n_reads_processed=0;
  int32_t seq_buffer_length;
  tmap_refseq_t *refseq=NULL;
  tmap_bwt_t *bwt=NULL;
  tmap_sa_t *sa=NULL;
  tmap_file_t *fp_reads=NULL;
  tmap_seq_io_t *seqio = NULL;
  tmap_seq_t **seq_buffer = NULL;
  tmap_map_sams_t **sams= NULL;
  tmap_shm_t *shm = NULL;
  int32_t reads_queue_size;
  
  if(NULL == opt->fn_reads) {
      tmap_progress_set_verbosity(0); 
  }

  // adjust opt for opt->score_match
  opt->score_thr *= opt->score_match;

  // For suffix search we need the reverse bwt/sa and forward refseq
  if(0 == opt->shm_key) {
      tmap_progress_print("reading in reference data");
      refseq = tmap_refseq_read(opt->fn_fasta, 0);
      bwt = tmap_bwt_read(opt->fn_fasta, 1);
      sa = tmap_sa_read(opt->fn_fasta, 1);
      tmap_progress_print2("reference data read in");
  }
  else {
      tmap_progress_print("retrieving reference data from shared memory");
      shm = tmap_shm_init(opt->shm_key, 0, 0);
      if(NULL == (refseq = tmap_refseq_shm_unpack(tmap_shm_get_buffer(shm, TMAP_SHM_LISTING_REFSEQ)))) {
          tmap_error("the packed reference sequence was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (bwt = tmap_bwt_shm_unpack(tmap_shm_get_buffer(shm, TMAP_SHM_LISTING_REV_BWT)))) {
          tmap_error("the reverse BWT string was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (sa = tmap_sa_shm_unpack(tmap_shm_get_buffer(shm, TMAP_SHM_LISTING_REV_SA)))) {
          tmap_error("the reverse SA was not found in shared memory", Exit, SharedMemoryListing);
      }
      tmap_progress_print2("reference data retrieved from shared memory");
  }

  // Set the seed length
  if(-1 == opt->seed_length) {
      opt->seed_length = tmap_map3_get_seed_length(refseq->len);
      tmap_progress_print("setting the seed length to %d", opt->seed_length);
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
      if(seq_buffer_length < num_threads * TMAP_MAP3_THREAD_BLOCK_SIZE) {
          num_threads = 1 + (seq_buffer_length / TMAP_MAP3_THREAD_BLOCK_SIZE);
      }
      tmap_map3_read_lock_low = 0; // ALWAYS set before running threads 
      if(1 == num_threads) {
          tmap_map3_core_worker(seq_buffer, sams, seq_buffer_length, refseq, bwt, sa, 0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          tmap_map3_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = tmap_calloc(num_threads, sizeof(pthread_t), "threads");
          thread_data = tmap_calloc(num_threads, sizeof(tmap_map3_thread_data_t), "thread_data");

          for(i=0;i<num_threads;i++) {
              thread_data[i].seq_buffer = seq_buffer;
              thread_data[i].seq_buffer_length = seq_buffer_length;
              thread_data[i].sams = sams;
              thread_data[i].refseq = refseq;
              thread_data[i].bwt = bwt;
              thread_data[i].sa = sa;;
              thread_data[i].tid = i;
              thread_data[i].opt = opt; 
              if(0 != pthread_create(&threads[i], &attr, tmap_map3_core_thread_worker, &thread_data[i])) {
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
      tmap_map3_core_worker(seq_buffer, sams, seq_buffer_length, refseq, bwt, sa, 0, opt);
#endif

      if(-1 != opt->reads_queue_size) {
          tmap_progress_print("writing alignments");
      }
      for(i=0;i<seq_buffer_length;i++) {
          // print
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
  tmap_bwt_destroy(bwt);
  tmap_sa_destroy(sa);
  tmap_seq_io_destroy(seqio);
  if(0 < opt->shm_key) {
      tmap_shm_destroy(shm, 0);
  }
}

int 
tmap_map3_main(int argc, char *argv[])
{
  tmap_map_opt_t *opt = NULL;

  // random seed
  srand48(0); 

  // init opt
  opt = tmap_map_opt_init(TMAP_MAP_ALGO_MAP3);

  // get options
  if(1 != tmap_map_opt_parse(argc, argv, opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return tmap_map_opt_usage(opt);
  }
  else { 
      // check command line arguments
      tmap_map_opt_check(opt);
  }

  // run map3
  tmap_map3_core(opt);

  // destroy opt
  tmap_map_opt_destroy(opt);

  tmap_progress_print2("terminating successfully");

  return 0;
}
