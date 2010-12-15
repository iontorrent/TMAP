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

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "../util/fmap_progress.h"
#include "../util/fmap_sam.h"
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_sa.h"
#include "../io/fmap_seq_io.h"
#include "../server/fmap_shm.h"
#include "../sw/fmap_fsw.h"
#include "fmap_map_util.h"
#include "fmap_map2_mempool.h"
#include "fmap_map2_aux.h"
#include "fmap_map2_core.h"
#include "fmap_map2.h"

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t fmap_map2_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t fmap_map2_read_lock_low = 0;
#define FMAP_MAP2_THREAD_BLOCK_SIZE 512
#endif

static void
fmap_map2_core_worker(fmap_seq_t **seq_buffer, int32_t seq_buffer_length, fmap_map_sams_t **sams,
                      fmap_refseq_t *refseq, fmap_bwt_t *bwt[2], fmap_sa_t *sa[2],
                      int32_t tid, fmap_map_opt_t * opt)
{
  int32_t low, high;
  fmap_map2_global_mempool_t *pool = NULL;

  pool = fmap_map2_global_mempool_init();

  low = 0;
  while(low < seq_buffer_length) {
#ifdef HAVE_LIBPTHREAD
      if(1 < opt->num_threads) {
          pthread_mutex_lock(&fmap_map2_read_lock);

          // update bounds
          low = fmap_map2_read_lock_low;
          fmap_map2_read_lock_low += FMAP_MAP2_THREAD_BLOCK_SIZE;
          high = low + FMAP_MAP2_THREAD_BLOCK_SIZE;
          if(seq_buffer_length < high) {
              high = seq_buffer_length;
          }

          pthread_mutex_unlock(&fmap_map2_read_lock);
      }
      else {
          high = seq_buffer_length; // process all
      }
#else
      high = seq_buffer_length; // process all
#endif
      while(low<high) {

          fmap_seq_t *seq = NULL;

          // Clone the buffer
          seq = fmap_seq_clone(seq_buffer[low]);
          // Adjust for SFF
          fmap_seq_remove_key_sequence(seq);

          // process
          sams[low] = fmap_map2_aux_core(opt, seq, refseq, bwt, sa, pool);

          // filter
          fmap_map_sams_filter(sams[low], opt->aln_output_mode);

          // re-align the alignments in flow-space
          /*
          if(FMAP_SEQ_TYPE_SFF == seq_buffer[low]->type) {
              fmap_map_util_fsw(seq_buffer[low]->data.sff, 
                                sams[low], refseq, 
                                opt->bw, opt->aln_global, opt->score_thr,
                                opt->score_match, opt->pen_mm, opt->pen_gapo,
                                opt->pen_gape, opt->fscore);
          }
          */

          // destroy
          fmap_seq_destroy(seq);

          // next
          low++;
      }
  }

  fmap_map2_global_mempool_destroy(pool);
}

static void *
fmap_map2_core_thread_worker(void *arg)
{
  fmap_map2_thread_data_t *data = (fmap_map2_thread_data_t*)arg;

  fmap_map2_core_worker(data->seq_buffer, data->seq_buffer_length, data->sams,
                        data->refseq, data->bwt, data->sa,
                        data->tid, data->opt);

  return arg;
}

static void
fmap_map2_core(fmap_map_opt_t *opt)
{
  uint32_t i, n_reads_processed=0;
  int32_t seq_buffer_length;
  double scalar;
  fmap_refseq_t *refseq = NULL;
  fmap_bwt_t *bwt[2] = {NULL, NULL};
  fmap_sa_t *sa[2] = {NULL, NULL};
  fmap_file_t *fp_reads=NULL;
  fmap_seq_io_t *seqio = NULL;
  fmap_shm_t *shm = NULL;
  fmap_seq_t **seq_buffer = NULL;
  fmap_map_sams_t **sams = NULL;
  int32_t reads_queue_size;

  if(NULL == opt->fn_reads) {
      fmap_progress_set_verbosity(0); 
  }

  scalar = opt->score_match / log(opt->yita);
  /*
     fmap_progress_print( "mismatch: %lf, gap_open: %lf, gap_ext: %lf",
     exp(-opt->pen_mm / scalar) / opt->yita,
     exp(-opt->pen_gapo / scalar),
     exp(-opt->pen_gape / scalar));
     */

  // adjust opt for opt->score_match
  opt->score_thr *= opt->score_match;
  opt->length_coef *= opt->score_match;

  // get the data
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

  // allocate the buffer
  if(-1 == opt->reads_queue_size) {
      reads_queue_size = 1;
  }
  else {
      reads_queue_size = opt->reads_queue_size;
  }
  seq_buffer = fmap_malloc(sizeof(fmap_seq_t*)*reads_queue_size, "seq_buffer");
  sams = fmap_malloc(sizeof(fmap_map_sams_t*)*reads_queue_size, "alnseqs");

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
      if(seq_buffer_length < num_threads * FMAP_MAP2_THREAD_BLOCK_SIZE) {
          num_threads = 1 + (seq_buffer_length / FMAP_MAP2_THREAD_BLOCK_SIZE);
      }
      fmap_map2_read_lock_low = 0; // ALWAYS set before running threads 
      if(1 == num_threads) {
          fmap_map2_core_worker(seq_buffer, seq_buffer_length, sams,
                                refseq, bwt, sa, 0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          fmap_map2_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = fmap_calloc(num_threads, sizeof(pthread_t), "threads");
          thread_data = fmap_calloc(num_threads, sizeof(fmap_map2_thread_data_t), "thread_data");

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
              if(0 != pthread_create(&threads[i], &attr, fmap_map2_core_thread_worker, &thread_data[i])) {
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
      fmap_map2_core_worker(seq_buffer, seq_buffer_length, sams,
                            refseq, bwt, sa, 0, opt);
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
  fmap_file_fclose(fp_reads);
  fmap_file_fclose(fmap_file_stdout);

  // free memory
  free(sams);
  for(i=0;i<reads_queue_size;i++) {
      fmap_seq_destroy(seq_buffer[i]);
  }
  free(seq_buffer);
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
fmap_map2_main(int argc, char *argv[])
{
  fmap_map_opt_t *opt = NULL;

  // random seed
  srand48(0);

  // init opt
  opt = fmap_map_opt_init(FMAP_MAP_ALGO_MAP2);

  // get options
  if(1 != fmap_map_opt_parse(argc, argv, opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return fmap_map_opt_usage(opt);
  }
  else { 
      // check command line arguments
      fmap_map_opt_check(opt);
  }

  // run map2
  fmap_map2_core(opt);

  // destroy opt
  fmap_map_opt_destroy(opt);

  fmap_progress_print2("terminating successfully");

  return 0;
}
