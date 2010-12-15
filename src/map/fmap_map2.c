#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
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
#include "../util/tmap_sam.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_sa.h"
#include "../io/tmap_seq_io.h"
#include "../server/tmap_shm.h"
#include "../sw/tmap_fsw.h"
#include "tmap_map_util.h"
#include "tmap_map2_mempool.h"
#include "tmap_map2_aux.h"
#include "tmap_map2_core.h"
#include "tmap_map2.h"

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t tmap_map2_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t tmap_map2_read_lock_low = 0;
#define TMAP_MAP2_THREAD_BLOCK_SIZE 512
#endif

static void
tmap_map2_core_worker(tmap_seq_t **seq_buffer, int32_t seq_buffer_length, tmap_map_sams_t **sams,
                      tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2],
                      int32_t tid, tmap_map_opt_t * opt)
{
  int32_t low, high;
  tmap_map2_global_mempool_t *pool = NULL;

  pool = tmap_map2_global_mempool_init();

  low = 0;
  while(low < seq_buffer_length) {
#ifdef HAVE_LIBPTHREAD
      if(1 < opt->num_threads) {
          pthread_mutex_lock(&tmap_map2_read_lock);

          // update bounds
          low = tmap_map2_read_lock_low;
          tmap_map2_read_lock_low += TMAP_MAP2_THREAD_BLOCK_SIZE;
          high = low + TMAP_MAP2_THREAD_BLOCK_SIZE;
          if(seq_buffer_length < high) {
              high = seq_buffer_length;
          }

          pthread_mutex_unlock(&tmap_map2_read_lock);
      }
      else {
          high = seq_buffer_length; // process all
      }
#else
      high = seq_buffer_length; // process all
#endif
      while(low<high) {

          tmap_seq_t *seq = NULL;

          // Clone the buffer
          seq = tmap_seq_clone(seq_buffer[low]);
          // Adjust for SFF
          tmap_seq_remove_key_sequence(seq);

          // process
          sams[low] = tmap_map2_aux_core(opt, seq, refseq, bwt, sa, pool);

          // filter
          tmap_map_sams_filter(sams[low], opt->aln_output_mode);

          // re-align the alignments in flow-space
          /*
          if(TMAP_SEQ_TYPE_SFF == seq_buffer[low]->type) {
              tmap_map_util_fsw(seq_buffer[low]->data.sff, 
                                sams[low], refseq, 
                                opt->bw, opt->aln_global, opt->score_thr,
                                opt->score_match, opt->pen_mm, opt->pen_gapo,
                                opt->pen_gape, opt->fscore);
          }
          */

          // destroy
          tmap_seq_destroy(seq);

          // next
          low++;
      }
  }

  tmap_map2_global_mempool_destroy(pool);
}

static void *
tmap_map2_core_thread_worker(void *arg)
{
  tmap_map2_thread_data_t *data = (tmap_map2_thread_data_t*)arg;

  tmap_map2_core_worker(data->seq_buffer, data->seq_buffer_length, data->sams,
                        data->refseq, data->bwt, data->sa,
                        data->tid, data->opt);

  return arg;
}

static void
tmap_map2_core(tmap_map_opt_t *opt)
{
  uint32_t i, n_reads_processed=0;
  int32_t seq_buffer_length;
  double scalar;
  tmap_refseq_t *refseq = NULL;
  tmap_bwt_t *bwt[2] = {NULL, NULL};
  tmap_sa_t *sa[2] = {NULL, NULL};
  tmap_file_t *fp_reads=NULL;
  tmap_seq_io_t *seqio = NULL;
  tmap_shm_t *shm = NULL;
  tmap_seq_t **seq_buffer = NULL;
  tmap_map_sams_t **sams = NULL;
  int32_t reads_queue_size;

  if(NULL == opt->fn_reads) {
      tmap_progress_set_verbosity(0); 
  }

  scalar = opt->score_match / log(opt->yita);
  /*
     tmap_progress_print( "mismatch: %lf, gap_open: %lf, gap_ext: %lf",
     exp(-opt->pen_mm / scalar) / opt->yita,
     exp(-opt->pen_gapo / scalar),
     exp(-opt->pen_gape / scalar));
     */

  // adjust opt for opt->score_match
  opt->score_thr *= opt->score_match;
  opt->length_coef *= opt->score_match;

  // get the data
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

  // allocate the buffer
  if(-1 == opt->reads_queue_size) {
      reads_queue_size = 1;
  }
  else {
      reads_queue_size = opt->reads_queue_size;
  }
  seq_buffer = tmap_malloc(sizeof(tmap_seq_t*)*reads_queue_size, "seq_buffer");
  sams = tmap_malloc(sizeof(tmap_map_sams_t*)*reads_queue_size, "alnseqs");

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
      if(seq_buffer_length < num_threads * TMAP_MAP2_THREAD_BLOCK_SIZE) {
          num_threads = 1 + (seq_buffer_length / TMAP_MAP2_THREAD_BLOCK_SIZE);
      }
      tmap_map2_read_lock_low = 0; // ALWAYS set before running threads 
      if(1 == num_threads) {
          tmap_map2_core_worker(seq_buffer, seq_buffer_length, sams,
                                refseq, bwt, sa, 0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          tmap_map2_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = tmap_calloc(num_threads, sizeof(pthread_t), "threads");
          thread_data = tmap_calloc(num_threads, sizeof(tmap_map2_thread_data_t), "thread_data");

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
              if(0 != pthread_create(&threads[i], &attr, tmap_map2_core_thread_worker, &thread_data[i])) {
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
      tmap_map2_core_worker(seq_buffer, seq_buffer_length, sams,
                            refseq, bwt, sa, 0, opt);
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
  tmap_file_fclose(fp_reads);
  tmap_file_fclose(tmap_file_stdout);

  // free memory
  free(sams);
  for(i=0;i<reads_queue_size;i++) {
      tmap_seq_destroy(seq_buffer[i]);
  }
  free(seq_buffer);
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

int 
tmap_map2_main(int argc, char *argv[])
{
  tmap_map_opt_t *opt = NULL;

  // random seed
  srand48(0);

  // init opt
  opt = tmap_map_opt_init(TMAP_MAP_ALGO_MAP2);

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

  // run map2
  tmap_map2_core(opt);

  // destroy opt
  tmap_map_opt_destroy(opt);

  tmap_progress_print2("terminating successfully");

  return 0;
}
