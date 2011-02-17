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
#include "../map/tmap_map_util.h"
#include "../map/tmap_map1.h"
#include "../map/tmap_map1_aux.h"
#include "../map/tmap_map2.h"
#include "../map/tmap_map2_aux.h"
#include "../map/tmap_map3.h"
#include "../map/tmap_map3_aux.h"
#include "../map/tmap_map_all.h"
#include "tmap_mappability.h"

/* Notes:
   We could avoid some computation give we know which algorithms to 
   run. This includes stacks, memory pools, as well as not loading 
   in all reference data.
   */

#ifdef HAVE_LIBPTHREAD
extern pthread_mutex_t tmap_map_all_read_lock;
extern int32_t tmap_map_all_read_lock_low;
#endif

static void 
tmap_mappability_core(tmap_map_opt_t *opt)
{
  uint32_t i, j, n_reads_processed=0;
  int32_t seq_buffer_length;
  tmap_refseq_t *refseq=NULL;
  tmap_bwt_t *bwt[2]={NULL, NULL};
  tmap_sa_t *sa[2]={NULL, NULL};
  tmap_seq_t **seq_buffer = NULL;
  tmap_map_sams_t **sams = NULL;
  tmap_shm_t *shm = NULL;
  int32_t reads_queue_size;
  uint32_t tid, pos, strand;
  uint32_t tid_start, tid_end, pos_start, pos_end;
  int32_t read_length;
  
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

  if(NULL == opt->region) {
      tid_start = 0;
      tid_end = refseq->num_annos-1;
      pos_start = 0;
      pos_end = refseq->annos[tid_end].len-1;
  }
  else {
      i=0;
      while(i < strlen(opt->region)) {
          if(':' == opt->region[i]) {
              break;
          }
          i++;
      }
      for(j=0;j<refseq->num_annos;j++) {
          if(i == (int32_t)refseq->annos[j].name->l &&
             0 == strncmp(opt->region, (char*)refseq->annos[j].name->s, i)) {
              break;
          }
      }
      if(refseq->num_annos <= j) {
          tmap_file_fprintf(tmap_file_stderr, "region=[%s]\n", opt->region);
          tmap_error("Could not find the given contig in the reference", Exit, OutOfRange);
      }
      tid_start = tid_end = j;
      if(i == strlen(opt->region)) {
          pos_start = 0;
          pos_end = refseq->annos[tid_start].len-1;
      }
      else {
          i++; // skip over the colon
          pos_start = atoi(opt->region + i) - 1;
          j = i;
          while(j < strlen(opt->region)) {
              if('-' == opt->region[j]) {
                  break;
              }
              j++;
          }
          if(j == strlen(opt->region)) {
              pos_end = refseq->annos[tid_end].len-1;
          }
          else {
              j++; // skip over the dash
              pos_end = atoi(opt->region + j) - 1;
          }
      }
  }
  if(refseq->num_annos <= tid_start) {
      tmap_file_fprintf(tmap_file_stderr, "region=[%s]\n", opt->region);
      tmap_error("region start contig out of range", Exit, OutOfRange);
  }
  if(refseq->num_annos <= tid_end) {
      tmap_file_fprintf(tmap_file_stderr, "region=[%s]\n", opt->region);
      tmap_error("region end contig out of range", Exit, OutOfRange);
  }
  if(refseq->annos[tid_start].len <= pos_start) {
      tmap_file_fprintf(tmap_file_stderr, "region=[%s]\n", opt->region);
      tmap_error("region start position out of range", Exit, OutOfRange);
  }
  if(refseq->annos[tid_end].len <= pos_end) {
      tmap_file_fprintf(tmap_file_stderr, "region=[%s]\n", opt->region);
      tmap_error("region end position out of range", Exit, OutOfRange);
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
  read_length = opt->read_length;
  // initialize the buffer
  for(i=0;i<reads_queue_size;i++) { 
      tmap_fq_t *fq = NULL;
      seq_buffer[i] = tmap_seq_init(TMAP_SEQ_TYPE_FQ);
      // initialize as fastqs
      fq = seq_buffer[i]->data.fq;
      tmap_string_destroy(fq->seq);
      fq->seq = tmap_string_init(read_length+1);
  }

  // Note: 'tmap_file_stdout' should not have been previously modified
  tmap_file_stdout = tmap_file_fdopen(fileno(stdout), "wb", opt->output_compr);

  // SAM header
  tmap_sam_print_header(tmap_file_stdout, refseq, NULL, opt->sam_rg, opt->sam_sff_tags, opt->argc, opt->argv);

  tmap_progress_print("processing reads");

  tid = tid_start;
  pos = pos_start;
  strand = 0;
  while(1) {
      uint32_t tmp_tid_start, tmp_pos_start, tmp_strand_start;
      uint32_t tmp_tid_end, tmp_pos_end, tmp_strand_end;
      tmp_tid_start = tid;
      tmp_pos_start = pos;
      tmp_strand_start = strand;
      tmp_tid_end = tid;
      tmp_pos_end = pos;
      tmp_strand_end = strand;
      if(-1 != opt->reads_queue_size) {
          tmap_progress_print("simulating reads");
      }
      seq_buffer_length = 0;
      while(tid < tid_end || (tid == tid_end && pos <= pos_end && pos + read_length <= refseq->annos[tid_end].len)) {
          // TODO: simulate the reads
          tmap_fq_t *fq = NULL;

          tmp_tid_end = tid;
          tmp_pos_end = pos;
          tmp_strand_end = strand;

          fq = seq_buffer[seq_buffer_length]->data.fq;
          tmap_string_lsprintf(fq->name, 0, "%s:%c:%d-%d", 
                               (char*)refseq->annos[tid].name->s, "+-"[strand], pos+1, pos+read_length);
          if(0 == strand) {
              for(i=0;i<read_length;i++) {
                  fq->seq->s[i] = "ACGTN"[tmap_refseq_seq_i(refseq, (pos+i+refseq->annos[tid].offset))];
              }
          }
          else {
              for(i=0;i<read_length;i++) {
                  fq->seq->s[read_length-i-1] = "TGCAN"[tmap_refseq_seq_i(refseq, (pos+i+refseq->annos[tid].offset))];
              }
          }
          fq->seq->s[read_length] = '\0';
          fq->seq->l = read_length;
          seq_buffer_length++;

          if(1 == strand) {
              strand = 0;
              pos++;
              if(pos == refseq->annos[tid].len - read_length) {
                  pos = 0;
                  tid++;
              }
          }
          else {
              strand++;
          }

          if(reads_queue_size <= seq_buffer_length) {
              break;
          }
      }
      if(-1 != opt->reads_queue_size) {
          tmap_progress_print2("simulated reads from %s:%d:%c to %s:%d:%c",
                               (char*)refseq->annos[tmp_tid_start].name->s, tmp_pos_start, "+-"[tmp_strand_start],
                               (char*)refseq->annos[tmp_tid_end].name->s, tmp_pos_end, "+-"[tmp_strand_end]);
      }
      
      if(0 == seq_buffer_length) {
          break;
      }
      
      if(-1 != opt->reads_queue_size) {
          tmap_progress_print("aligning reads");
      }

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

      if(tid_end < tid ||
         (tid_end == tid && 
          (pos_end + read_length < pos || 
           refseq->annos[tid_end].len < pos + read_length))) {
          break;
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
  if(0 < opt->shm_key) {
      tmap_shm_destroy(shm, 0);
  }
}

int 
tmap_mappability_main(int argc, char *argv[])
{
  tmap_map_opt_t *opt = NULL;

  // random seed
  srand48(0); 

  // init opt
  opt = tmap_map_opt_init(TMAP_MAP_ALGO_MAPPABILTY);

  opt->reads_format = TMAP_READS_FORMAT_FASTA;
      
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
  tmap_mappability_core(opt);

  // destroy opt
  tmap_map_opt_destroy(opt);

  tmap_progress_print2("terminating successfully");

  return 0;
}
