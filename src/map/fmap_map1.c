#include <stdlib.h>
#include <math.h>
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
#include "../index/tmap_bwt_gen.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_bwt_match.h"
#include "../index/tmap_sa.h"
#include "../io/tmap_seq_io.h"
#include "../server/tmap_shm.h"
#include "../sw/tmap_sw.h"
#include "../sw/tmap_fsw.h"
#include "tmap_map_util.h"
#include "tmap_map1_aux.h"
#include "tmap_map1.h"

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t tmap_map1_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t tmap_map1_read_lock_low = 0;
#define TMAP_MAP1_THREAD_BLOCK_SIZE 512
#endif

static int32_t g_log_n[256];

static inline void
tmap_map1_set_g_log_n()
{
  int32_t i;
  for(i=1;i<256;i++) {
      g_log_n[i] = (int32_t)(4.343 * log(i) + 0.5);
  }
}

static inline uint8_t
tmap_map1_sam_mapq(int32_t num_best_sa, int32_t num_all_sa, int32_t max_mm, int32_t num_mm)
{
  int32_t n;

  if(1 < num_best_sa) {
      return 0; // multiple best hits
  }
  else if(max_mm == num_mm) {
      return 25; // maximum possible score
  }

  n = (num_all_sa - num_best_sa);
  if(0 == n) {
      return 37; // no second best hits
  }
  else if(255 < n) {
      n = 255;
  }

  // use MAQ-like mapping qualities
  return (23 < g_log_n[n])? 0 : 23 - g_log_n[n];
}

static inline void
tmap_map1_sams_mapq(tmap_map_sams_t *sams, tmap_map_opt_t *opt)
{
  int32_t i;
  int32_t num_best_sa, num_best, num_all_sa;

  if(0 == sams->n) {
      return;
  }

  //Note: assumes that the alignments are sorted by increasing score
  num_best = num_best_sa = num_all_sa = 0;
  for(i=0;i<sams->n;i++) {
      if(sams->sams[0].score < sams->sams[i].score) {
          break;
      }
      num_best++;
      num_best_sa++;
  }
  for(i=0;i<sams->n;i++) {
      num_all_sa++;
  }
  for(i=0;i<num_best;i++) {
      sams->sams[i].mapq = tmap_map1_sam_mapq(num_best_sa, num_all_sa, opt->max_mm, sams->sams[i].aux.map1_aux->n_mm);
  }
  for(i=num_best;i<sams->n;i++) {
      sams->sams[i].mapq = 0;
  }
}

static void
tmap_map1_core_worker(tmap_seq_t **seq_buffer, int32_t seq_buffer_length, tmap_map_sams_t **sams,
                      tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa, 
                      int32_t tid, tmap_map_opt_t *opt)
{
  int32_t low = 0, high;
  tmap_bwt_match_width_t *width[2]={NULL,NULL}, *seed_width[2]={NULL,NULL};
  int32_t width_length = 0;
  tmap_map1_aux_stack_t *stack = NULL;

  // for calculating mapping qualities
  tmap_map1_set_g_log_n();

  seed_width[0] = tmap_calloc(opt->seed_length, sizeof(tmap_bwt_match_width_t), "seed_width[0]");
  seed_width[1] = tmap_calloc(opt->seed_length, sizeof(tmap_bwt_match_width_t), "seed_width[1]");

  stack = tmap_map1_aux_stack_init();

  while(low < seq_buffer_length) {
#ifdef HAVE_LIBPTHREAD
      if(1 < opt->num_threads) {
          pthread_mutex_lock(&tmap_map1_read_lock);

          // update bounds
          low = tmap_map1_read_lock_low;
          tmap_map1_read_lock_low += TMAP_MAP1_THREAD_BLOCK_SIZE;
          high = low + TMAP_MAP1_THREAD_BLOCK_SIZE;
          if(seq_buffer_length < high) {
              high = seq_buffer_length; 
          }

          pthread_mutex_unlock(&tmap_map1_read_lock);
      }
      else {
          high = seq_buffer_length; // process all
      }
#else 
      high = seq_buffer_length; // process all
#endif
      while(low<high) {
          tmap_map_opt_t opt_local = (*opt); // copy over values
          tmap_seq_t *seq[2]={NULL, NULL}, *orig_seq=NULL;
          orig_seq = seq_buffer[low];
          tmap_string_t *bases[2]={NULL, NULL};

          // clone the sequence 
          seq[0] = tmap_seq_clone(orig_seq);
          seq[1] = tmap_seq_clone(orig_seq);

          // Adjust for SFF
          tmap_seq_remove_key_sequence(seq[0]);
          tmap_seq_remove_key_sequence(seq[1]);

          // reverse compliment
          tmap_seq_reverse_compliment(seq[1]);

          // convert to integers
          tmap_seq_to_int(seq[0]);
          tmap_seq_to_int(seq[1]);

          // get bases
          bases[0] = tmap_seq_get_bases(seq[0]);
          bases[1] = tmap_seq_get_bases(seq[1]);

          // remember to round up
          opt_local.max_mm = (opt->max_mm < 0) ? (int)(0.99 + opt->max_mm_frac * bases[0]->l) : opt->max_mm; 
          opt_local.max_gape = (opt->max_gape < 0) ? (int)(0.99 + opt->max_gape_frac * bases[0]->l) : opt->max_gape; 
          opt_local.max_gapo = (opt->max_gapo < 0) ? (int)(0.99 + opt->max_gapo_frac * bases[0]->l) : opt->max_gapo; 
          if(width_length < bases[0]->l) {
              free(width[0]); free(width[1]);
              width_length = bases[0]->l;
              width[0] = tmap_calloc(width_length, sizeof(tmap_bwt_match_width_t), "width[0]");
              width[1] = tmap_calloc(width_length, sizeof(tmap_bwt_match_width_t), "width[1]");
          }
          tmap_bwt_match_cal_width(bwt[0], bases[0]->l, bases[0]->s, width[0]);
          tmap_bwt_match_cal_width(bwt[0], bases[1]->l, bases[1]->s, width[1]);

          if(bases[0]->l < opt->seed_length) {
              opt_local.seed_length = -1;
          }
          else {
              tmap_bwt_match_cal_width(bwt[0], opt->seed_length, bases[0]->s, seed_width[0]);
              tmap_bwt_match_cal_width(bwt[0], opt->seed_length, bases[1]->s, seed_width[1]);
          }

          sams[low] = tmap_map1_aux_core(seq, refseq, bwt[1], sa, width, (0 < opt_local.seed_length) ? seed_width : NULL, &opt_local, stack);
                  
          // adjust map1 scoring, since it does not consider opt->score_match
          tmap_map_util_map1_adjust_score(sams[low], opt->score_match, opt->pen_mm, opt->pen_gapo, opt->pen_gape);

          // mapping quality
          tmap_map1_sams_mapq(sams[low], opt);

          // filter alignments
          tmap_map_sams_filter(sams[low], opt->aln_output_mode);

          // re-align the alignments in flow-space
          /*
          if(TMAP_SEQ_TYPE_SFF == seq_buffer[low]->type) {
              tmap_map_util_fsw(seq_buffer[low]->data.sff, 
                                sams[low], refseq, 
                                TMAP_MAP1_FSW_BW, 1, INT32_MIN,
                                0, opt->pen_mm, opt->pen_gapo,
                                opt->pen_gape, opt->fscore);
          }
          */

          // destroy
          tmap_seq_destroy(seq[0]);
          tmap_seq_destroy(seq[1]);

          // next
          low++;
      }
  }

  tmap_map1_aux_stack_destroy(stack);
  free(seed_width[0]);
  free(seed_width[1]);
  free(width[0]);
  free(width[1]);
}

static void *
tmap_map1_core_thread_worker(void *arg)
{
  tmap_map1_thread_data_t *thread_data = (tmap_map1_thread_data_t*)arg;

  tmap_map1_core_worker(thread_data->seq_buffer, thread_data->seq_buffer_length, thread_data->sams,
                        thread_data->refseq, thread_data->bwt, thread_data->sa, 
                        thread_data->tid, thread_data->opt);

  return arg;
}

static void 
tmap_map1_core(tmap_map_opt_t *opt)
{
  uint32_t i, n_reads_processed=0;
  int32_t seq_buffer_length;
  tmap_refseq_t *refseq=NULL;
  tmap_bwt_t *bwt[2]={NULL,NULL};
  tmap_sa_t *sa=NULL;
  tmap_file_t *fp_reads=NULL;
  tmap_seq_io_t *seqio = NULL;
  tmap_seq_t **seq_buffer = NULL;
  tmap_map_sams_t **sams = NULL;
  tmap_shm_t *shm = NULL;
  int32_t reads_queue_size;

  if(NULL == opt->fn_reads) {
      tmap_progress_set_verbosity(0); 
  }

  // For suffix search we need the reverse bwt/sa and forward refseq
  if(0 == opt->shm_key) {
      tmap_progress_print("reading in reference data");
      refseq = tmap_refseq_read(opt->fn_fasta, 0);
      bwt[0] = tmap_bwt_read(opt->fn_fasta, 0);
      bwt[1] = tmap_bwt_read(opt->fn_fasta, 1);
      sa = tmap_sa_read(opt->fn_fasta, 1);
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
      if(NULL == (sa = tmap_sa_shm_unpack(tmap_shm_get_buffer(shm, TMAP_SHM_LISTING_REV_SA)))) {
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
  sams = tmap_malloc(sizeof(tmap_map_sams_t*)*reads_queue_size, "sams");

  if(NULL == opt->fn_reads) {
      fp_reads = tmap_file_fdopen(fileno(stdin), "rb", opt->input_compr);
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
      tmap_map1_read_lock_low = 0; // ALWAYS set before running threads 
      if(seq_buffer_length < num_threads * TMAP_MAP1_THREAD_BLOCK_SIZE) {
          num_threads = 1 + (seq_buffer_length / TMAP_MAP1_THREAD_BLOCK_SIZE);
      }
      if(1 == num_threads) {
          tmap_map1_core_worker(seq_buffer, seq_buffer_length, sams, refseq, bwt, sa, 0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          tmap_map1_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = tmap_calloc(num_threads, sizeof(pthread_t), "threads");
          thread_data = tmap_calloc(num_threads, sizeof(tmap_map1_thread_data_t), "thread_data");

          for(i=0;i<num_threads;i++) {
              thread_data[i].seq_buffer = seq_buffer;
              thread_data[i].seq_buffer_length = seq_buffer_length;
              thread_data[i].sams = sams;
              thread_data[i].refseq = refseq;
              thread_data[i].bwt[0] = bwt[0];
              thread_data[i].bwt[1] = bwt[1];
              thread_data[i].sa = sa;
              thread_data[i].tid = i;
              thread_data[i].opt = opt; 
              if(0 != pthread_create(&threads[i], &attr, tmap_map1_core_thread_worker, &thread_data[i])) {
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
      tmap_map1_core_worker(seq_buffer, seq_buffer_length, sams, refseq, bwt, sa, 0, opt);
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

  // close the output
  tmap_file_fclose(tmap_file_stdout);

  // free memory
  for(i=0;i<reads_queue_size;i++) {
      tmap_seq_destroy(seq_buffer[i]);
  }
  free(seq_buffer);
  free(sams);
  tmap_file_fclose(fp_reads);
  tmap_refseq_destroy(refseq);
  tmap_bwt_destroy(bwt[0]);
  tmap_bwt_destroy(bwt[1]);
  tmap_sa_destroy(sa);
  tmap_seq_io_destroy(seqio);
  if(0 < opt->shm_key) {
      tmap_shm_destroy(shm, 0);
  }
}

int 
tmap_map1_main(int argc, char *argv[])
{
  tmap_map_opt_t *opt = NULL;

  // random seed
  srand48(0); 

  // init opt
  opt = tmap_map_opt_init(TMAP_MAP_ALGO_MAP1);

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

  // run map1
  tmap_map1_core(opt);

  // destroy opt
  tmap_map_opt_destroy(opt);

  tmap_progress_print2("terminating successfully");

  return 0;
}
