#include <stdlib.h>
#include <math.h>
#include <config.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif
#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "../util/fmap_progress.h"
#include "../util/fmap_sam.h"
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt_gen.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_bwt_match.h"
#include "../index/fmap_sa.h"
#include "../io/fmap_seq_io.h"
#include "../server/fmap_shm.h"
#include "fmap_map_util.h"
#include "fmap_map1_aux.h"
#include "fmap_map1.h"

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t fmap_map1_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t fmap_map1_read_lock_low = 0;
#define FMAP_MAP1_THREAD_BLOCK_SIZE 512
#endif

static int32_t g_log_n[256];

static inline void
fmap_map1_set_g_log_n()
{
  int32_t i;
  for(i=1;i<256;i++) {
      g_log_n[i] = (int32_t)(4.343 * log(i) + 0.5);
  }
}

static inline uint8_t
fmap_map1_aln_mapq(int32_t num_best_sa, int32_t num_all_sa, int32_t max_mm, int32_t num_mm)
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
fmap_map1_hit_overwrite(fmap_map1_hit_t *dest, fmap_map1_hit_t *src)
{
  // cigar
  free(dest->cigar);
  dest->cigar = src->cigar;
  dest->n_cigar = src->n_cigar;
  src->cigar = NULL;
  src->n_cigar = 0;

  // rest
  dest->score = src->score;
  dest->n_mm = src->n_mm;
  dest->n_gapo = src->n_gapo;
  dest->n_gape = src->n_gape;
  dest->mapq = src->mapq;
  dest->strand = src->strand;
  dest->k = src->k;
  dest->l = src->l;
}

inline fmap_map1_aln_t *
fmap_map1_aln_init()
{
  return fmap_calloc(1, sizeof(fmap_map1_aln_t), "return");
}

inline void
fmap_map1_aln_destroy(fmap_map1_aln_t *aln)
{
  int32_t i;
  for(i=0;i<aln->n;i++) {
      free(aln->hits[i].cigar);
  }
  free(aln->hits);
  free(aln);
}

static inline void
fmap_map1_aln_filter(fmap_map1_opt_t *opt, fmap_map1_aln_t *alns, int32_t max_mm)
{
  int32_t i, was_rand = 0;
  int32_t num_best_sa, num_best, num_all_sa;

  if(0 == alns->n) {
      return;
  }

  //Note: assumes that the alignments are sorted by increasing score
  num_best = num_best_sa = num_all_sa = 0;
  for(i=0;i<alns->n;i++) {
      if(alns->hits[0].score < alns->hits[i].score) {
          break;
      }
      num_best++;
      num_best_sa += alns->hits[i].l - alns->hits[i].k + 1;
  }
  for(i=0;i<alns->n;i++) {
      num_all_sa += alns->hits[i].l - alns->hits[i].k + 1;
  }

  if(FMAP_MAP1_ALN_OUTPUT_MODE_ALL == opt->aln_output_mode) { // all hits
      for(i=0;i<alns->n;i++) {
          alns->hits[i].mapq = 0;
      }
      alns->hits[0].mapq = fmap_map1_aln_mapq(num_best_sa, num_all_sa, max_mm, alns->hits[0].n_mm);
  }
  else if(FMAP_MAP1_ALN_OUTPUT_MODE_BEST == opt->aln_output_mode  // unique best hit 
          || FMAP_MAP1_ALN_OUTPUT_MODE_BEST_RAND == opt->aln_output_mode) { // random best hit

      if(1 < num_best_sa && FMAP_MAP1_ALN_OUTPUT_MODE_BEST_RAND == opt->aln_output_mode) { // pick a random one
          int32_t best_index = 0, rand;
          rand = drand48() * num_best_sa; 
          for(i=0;i<num_best;i++) { // get the "rand"th best score
              if(rand < alns->hits[i].l - alns->hits[i].k + 1) {
                  best_index = i;
                  break;
              }
              rand -= alns->hits[i].l - alns->hits[i].k + 1;
          }
          // copy to the front
          if(0 < best_index) {
              fmap_map1_hit_overwrite(&alns->hits[0], &alns->hits[best_index]);
          }
          alns->hits[0].k += rand;
          if(alns->hits[0].l < alns->hits[0].k) {
              fmap_error("control reach an unexpected point", Exit, OutOfRange);
          }
          alns->hits[0].l = alns->hits[0].k;
          num_best_sa = 1;
          num_best = 1; 
          // WARNING: make sure mapping quality is zero after this
          was_rand = 1;
      }

      if(1 < num_best_sa) {
          // free them all
          for(i=0;i<alns->n;i++) {
              free(alns->hits[i].cigar);
          }
          free(alns->hits);
          alns->n = 0;
      }
      else {
          if(FMAP_MAP1_ALN_OUTPUT_MODE_BEST_RAND == opt->aln_output_mode && 1 == was_rand) {
              alns->hits[0].mapq = 0;
          }
          else {
              alns->hits[0].mapq = fmap_map1_aln_mapq(num_best_sa, num_all_sa, max_mm, alns->hits[0].n_mm);
          }

          // free the rest
          for(i=1;i<alns->n;i++) { 
              free(alns->hits[i].cigar);
          }
          // save only the first
          alns->hits = fmap_realloc(alns->hits, sizeof(fmap_map1_hit_t*), "alns->hits");
          alns->n = 1;
      }
  }
  else if(FMAP_MAP1_ALN_OUTPUT_MODE_BEST_ALL == opt->aln_output_mode) { // all best hits
      for(i=0;i<num_best;i++) {
          alns->hits[i].mapq = fmap_map1_aln_mapq(num_best_sa, num_all_sa, max_mm, alns->hits[i].n_mm);
      }
      if(num_best < alns->n) {
          for(i=num_best;i<alns->n;i++) { // free the rest
              free(alns->hits[i].cigar);
          }
          alns->hits = fmap_realloc(alns->hits, sizeof(fmap_map1_hit_t*)*num_best, "alns->hits");
          alns->n = num_best;
      }
  }
}

static inline int 
fmap_map1_print_sam(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_bwt_t *bwt, fmap_sa_t *sa, fmap_map1_hit_t *h)
{
  uint32_t i, j, n = 0;
      int32_t seq_len = 0;

      seq_len = fmap_seq_get_bases(seq)->l;

      if(FMAP_SEQ_TYPE_SFF == seq->type) {
          seq_len -= seq->data.sff->gheader->key_length; // soft clip the key sequence
      }

  for(i=h->k;i<=h->l;i++) {
      uint32_t pos = 0, seqid = 0, pacpos = 0; 
      int32_t aln_ref_l = 0;

      // get the number of non-inserted bases 
      for(j=0;j<h->n_cigar;j++) {
          switch((h->cigar[j] & 0xf)) {
            case BAM_CMATCH:
            case BAM_CDEL:
              aln_ref_l += (h->cigar[j] >> 4); break;
            default:
              break;
          }
      }

      // SA position to packed refseq position
      pacpos = bwt->seq_len - fmap_sa_pac_pos(sa, bwt, i) - aln_ref_l + 1;

      if(0 < fmap_refseq_pac2real(refseq, pacpos, seq_len, &seqid, &pos)) {

          fmap_sam_print_mapped(fmap_file_stdout, seq, refseq,
                                h->strand, seqid, pos-1, 
                                h->mapq, h->cigar, h->n_cigar, 
                                "\tAS:i:%d\tNM:i:%d\tXM:i:%d\tXO:i:%d\tXG:i:%d",
                                h->score,
                                (h->n_mm + h->n_gapo + h->n_gape),
                                h->n_mm, h->n_gapo, h->n_gape);

          n++;
      }
  }

  return n;
}

static void
fmap_map1_core_worker(fmap_seq_t **seq_buffer, int32_t seq_buffer_length, fmap_map1_aln_t **alns,
                      fmap_bwt_t *bwt[2], int32_t tid, fmap_map1_opt_t *opt)
{
  int32_t low = 0, high;
  fmap_bwt_match_width_t *width[2]={NULL,NULL}, *seed_width[2]={NULL,NULL};
  int32_t width_length = 0;
  fmap_map1_aux_stack_t *stack = NULL;

  // for calculating mapping qualities
  fmap_map1_set_g_log_n();

  seed_width[0] = fmap_calloc(opt->seed_length, sizeof(fmap_bwt_match_width_t), "seed_width[0]");
  seed_width[1] = fmap_calloc(opt->seed_length, sizeof(fmap_bwt_match_width_t), "seed_width[1]");

  stack = fmap_map1_aux_stack_init();

  while(low < seq_buffer_length) {
#ifdef HAVE_LIBPTHREAD
      if(1 < opt->num_threads) {
          pthread_mutex_lock(&fmap_map1_read_lock);

          // update bounds
          low = fmap_map1_read_lock_low;
          fmap_map1_read_lock_low += FMAP_MAP1_THREAD_BLOCK_SIZE;
          high = low + FMAP_MAP1_THREAD_BLOCK_SIZE;
          if(seq_buffer_length < high) {
              high = seq_buffer_length; 
          }

          pthread_mutex_unlock(&fmap_map1_read_lock);
      }
      else {
          high = seq_buffer_length; // process all
      }
#else 
      high = seq_buffer_length; // process all
#endif
      while(low<high) {
          fmap_map1_opt_t opt_local = (*opt); // copy over values
          fmap_seq_t *seq[2]={NULL, NULL}, *orig_seq=NULL;
          orig_seq = seq_buffer[low];
          fmap_string_t *bases[2]={NULL, NULL};

          // clone the sequence 
          seq[0] = fmap_seq_clone(orig_seq);
          seq[1] = fmap_seq_clone(orig_seq);

          // Adjust for SFF
          fmap_seq_remove_key_sequence(seq[0]);
          fmap_seq_remove_key_sequence(seq[1]);

          // reverse compliment
          fmap_seq_reverse_compliment(seq[1]);

          // convert to integers
          fmap_seq_to_int(seq[0]);
          fmap_seq_to_int(seq[1]);

          // get bases
          bases[0] = fmap_seq_get_bases(seq[0]);
          bases[1] = fmap_seq_get_bases(seq[1]);

          // remember to round up
          opt_local.max_mm = (opt->max_mm < 0) ? (int)(0.99 + opt->max_mm_frac * bases[0]->l) : opt->max_mm; 
          opt_local.max_gape = (opt->max_gape < 0) ? (int)(0.99 + opt->max_gape_frac * bases[0]->l) : opt->max_gape; 
          opt_local.max_gapo = (opt->max_gapo < 0) ? (int)(0.99 + opt->max_gapo_frac * bases[0]->l) : opt->max_gapo; 
          if(width_length < bases[0]->l) {
              free(width[0]); free(width[1]);
              width_length = bases[0]->l;
              width[0] = fmap_calloc(width_length, sizeof(fmap_bwt_match_width_t), "width[0]");
              width[1] = fmap_calloc(width_length, sizeof(fmap_bwt_match_width_t), "width[1]");
          }
          fmap_bwt_match_cal_width(bwt[0], bases[0]->l, bases[0]->s, width[0]);
          fmap_bwt_match_cal_width(bwt[0], bases[1]->l, bases[1]->s, width[1]);

          if(bases[0]->l < opt->seed_length) {
              opt_local.seed_length = -1;
          }
          else {
              fmap_bwt_match_cal_width(bwt[0], opt->seed_length, bases[0]->s, seed_width[0]);
              fmap_bwt_match_cal_width(bwt[0], opt->seed_length, bases[1]->s, seed_width[1]);
          }

          alns[low] = fmap_map1_aux_core(seq, bwt[1], width, (0 < opt_local.seed_length) ? seed_width : NULL, &opt_local, stack);

          // filter alignments
          fmap_map1_aln_filter(&opt_local, alns[low], opt->max_mm);

          // destroy
          fmap_seq_destroy(seq[0]);
          fmap_seq_destroy(seq[1]);

          // next
          low++;
      }
  }

  fmap_map1_aux_stack_destroy(stack);
  free(seed_width[0]);
  free(seed_width[1]);
  free(width[0]);
  free(width[1]);
}

static void *
fmap_map1_core_thread_worker(void *arg)
{
  fmap_map1_thread_data_t *thread_data = (fmap_map1_thread_data_t*)arg;

  fmap_map1_core_worker(thread_data->seq_buffer, thread_data->seq_buffer_length, thread_data->alns,
                        thread_data->bwt, thread_data->tid, thread_data->opt);

  return arg;
}

static void 
fmap_map1_core(fmap_map1_opt_t *opt)
{
  uint32_t i, j, n_reads_processed=0;
  int32_t seq_buffer_length;
  fmap_refseq_t *refseq=NULL;
  fmap_bwt_t *bwt[2]={NULL,NULL};
  fmap_sa_t *sa=NULL;
  fmap_file_t *fp_reads=NULL;
  fmap_seq_io_t *seqio = NULL;
  fmap_seq_t **seq_buffer = NULL;
  fmap_map1_aln_t **alns = NULL;
  fmap_shm_t *shm = NULL;
  int32_t reads_queue_size;

  if(NULL == opt->fn_reads) {
      fmap_progress_set_verbosity(0); 
  }

  // For suffix search we need the reverse bwt/sa and forward refseq
  if(0 == opt->shm_key) {
      fmap_progress_print("reading in reference data");
      refseq = fmap_refseq_read(opt->fn_fasta, 0);
      bwt[0] = fmap_bwt_read(opt->fn_fasta, 0);
      bwt[1] = fmap_bwt_read(opt->fn_fasta, 1);
      sa = fmap_sa_read(opt->fn_fasta, 1);
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
      if(NULL == (sa = fmap_sa_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_REV_SA)))) {
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
  alns = fmap_malloc(sizeof(fmap_map1_aln_t*)*reads_queue_size, "alns");

  if(NULL == opt->fn_reads) {
      fp_reads = fmap_file_fdopen(fileno(stdin), "rb", opt->input_compr);
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
  fmap_sam_print_header(fmap_file_stdout, refseq, seqio, opt->sam_rg, opt->argc, opt->argv);

  fmap_progress_print("processing reads");
  while(0 < (seq_buffer_length = fmap_seq_io_read_buffer(seqio, seq_buffer, reads_queue_size))) {

      // do alignment
#ifdef HAVE_LIBPTHREAD
      int32_t num_threads = opt->num_threads;
      fmap_map1_read_lock_low = 0; // ALWAYS set before running threads 
      if(seq_buffer_length < num_threads * FMAP_MAP1_THREAD_BLOCK_SIZE) {
          num_threads = 1 + (seq_buffer_length / FMAP_MAP1_THREAD_BLOCK_SIZE);
      }
      if(1 == num_threads) {
          fmap_map1_core_worker(seq_buffer, seq_buffer_length, alns, bwt, 0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          fmap_map1_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = fmap_calloc(num_threads, sizeof(pthread_t), "threads");
          thread_data = fmap_calloc(num_threads, sizeof(fmap_map1_thread_data_t), "thread_data");

          for(i=0;i<num_threads;i++) {
              thread_data[i].seq_buffer = seq_buffer;
              thread_data[i].seq_buffer_length = seq_buffer_length;
              thread_data[i].alns = alns;
              thread_data[i].bwt[0] = bwt[0];
              thread_data[i].bwt[1] = bwt[1];
              thread_data[i].tid = i;
              thread_data[i].opt = opt; 
              if(0 != pthread_create(&threads[i], &attr, fmap_map1_core_thread_worker, &thread_data[i])) {
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
      fmap_map1_core_worker(seq_buffer, seq_buffer_length, alns, bwt, 0, opt);
#endif

      if(-1 != opt->reads_queue_size) {
          fmap_progress_print("writing alignments");
      }
      for(i=0;i<seq_buffer_length;i++) {
          int32_t n_mapped = 0;
          fmap_map1_aln_t *a = alns[i];

          // print alignments
          if(0 < a->n) {
              for(j=0;j<a->n;j++) { 
                  n_mapped += fmap_map1_print_sam(seq_buffer[i], refseq, bwt[1], sa, &a->hits[j]);
              }
          }
          if(0 == n_mapped) {
              fmap_sam_print_unmapped(fmap_file_stdout, seq_buffer[i]);
          }

          // free alignments
          fmap_map1_aln_destroy(alns[i]);
          alns[i] = NULL;
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

  // close the output
  fmap_file_fclose(fmap_file_stdout);

  // free memory
  for(i=0;i<reads_queue_size;i++) {
      fmap_seq_destroy(seq_buffer[i]);
  }
  free(seq_buffer);
  free(alns);
  fmap_file_fclose(fp_reads);
  fmap_refseq_destroy(refseq);
  fmap_bwt_destroy(bwt[0]);
  fmap_bwt_destroy(bwt[1]);
  fmap_sa_destroy(sa);
  fmap_seq_io_destroy(seqio);
  if(0 < opt->shm_key) {
      fmap_shm_destroy(shm, 0);
  }
}

int 
fmap_map1_usage(fmap_map1_opt_t *opt)
{
  char *reads_format = fmap_get_reads_file_format_string(opt->reads_format);
  // Future options:
  // - adapter trimming ?
  // - homopolymer enumeration ?
  // - add option to try various seed offsets
  // - add option for "how many edits away" to search
  // - add an option to only output all alignments
  // - add an option to randomize best scoring alignments

  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Usage: %s map1 [options]", PACKAGE);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (required):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -f FILE     the FASTA reference file name [%s]\n", opt->fn_fasta);
  fmap_file_fprintf(fmap_file_stderr, "         -r FILE     the reads file name [%s]\n", (NULL == opt->fn_reads) ? "stdin" : opt->fn_reads);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (optional):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff) [%s]\n", reads_format);
  fmap_file_fprintf(fmap_file_stderr, "         -l INT      the k-mer length to seed CALs (-1 to disable) [%d]\n", opt->seed_length);
  fmap_file_fprintf(fmap_file_stderr, "         -k INT      maximum number of mismatches in the seed [%d]\n", opt->seed_max_mm);

  fmap_file_fprintf(fmap_file_stderr, "         -m NUM      maximum number of or (read length) fraction of mismatches");
  if(opt->max_mm < 0) fmap_file_fprintf(fmap_file_stderr, " [fraction: %lf]\n", opt->max_mm_frac);
  else fmap_file_fprintf(fmap_file_stderr, " [number: %d]\n", opt->max_mm); 

  fmap_file_fprintf(fmap_file_stderr, "         -o NUM      maximum number of or (read length) fraction of indel starts");
  if(opt->max_gapo < 0) fmap_file_fprintf(fmap_file_stderr, " [fraction: %lf]\n", opt->max_gapo_frac);
  else fmap_file_fprintf(fmap_file_stderr, " [number: %d]\n", opt->max_gapo); 

  fmap_file_fprintf(fmap_file_stderr, "         -e NUM      maximum number of or (read length) fraction of indel extensions");
  if(opt->max_gape < 0) fmap_file_fprintf(fmap_file_stderr, " [fraction: %lf]\n", opt->max_gape_frac);
  else fmap_file_fprintf(fmap_file_stderr, " [number: %d]\n", opt->max_gape); 

  fmap_file_fprintf(fmap_file_stderr, "         -M INT      the mismatch penalty [%d]\n", opt->pen_mm); 
  fmap_file_fprintf(fmap_file_stderr, "         -O INT      the indel start penalty [%d]\n", opt->pen_gapo); 
  fmap_file_fprintf(fmap_file_stderr, "         -E INT      the indel extend penalty [%d]\n", opt->pen_gape); 
  fmap_file_fprintf(fmap_file_stderr, "         -d INT      the maximum number of CALs to extend a deletion [%d]\n", opt->max_cals_del); 
  fmap_file_fprintf(fmap_file_stderr, "         -i INT      indels are not allowed within INT number of bps from the end of the read [%d]\n", opt->indel_ends_bound);
  fmap_file_fprintf(fmap_file_stderr, "         -b INT      stop searching when INT optimal CALs have been found [%d]\n", opt->max_best_cals);
  fmap_file_fprintf(fmap_file_stderr, "         -Q INT      maximum number of alignment nodes [%d]\n", opt->max_entries);
  fmap_file_fprintf(fmap_file_stderr, "         -q INT      the queue size for the reads (-1 disables) [%d]\n", opt->reads_queue_size);
  fmap_file_fprintf(fmap_file_stderr, "         -n INT      the number of threads [%d]\n", opt->num_threads);
  fmap_file_fprintf(fmap_file_stderr, "         -a INT      output filter [%d]\n", opt->aln_output_mode);
  fmap_file_fprintf(fmap_file_stderr, "                             0 - unique best hits\n");
  fmap_file_fprintf(fmap_file_stderr, "                             1 - random best hit\n");
  fmap_file_fprintf(fmap_file_stderr, "                             2 - all best hits\n");
  fmap_file_fprintf(fmap_file_stderr, "                             3 - all alignments\n");
  fmap_file_fprintf(fmap_file_stderr, "         -R STRING   the RG line in the SAM header [%s]\n", opt->sam_rg);
  fmap_file_fprintf(fmap_file_stderr, "         -j          the input is bz2 compressed (bzip2) [%s]\n",
                    (FMAP_FILE_BZ2_COMPRESSION == opt->input_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -z          the input is gz compressed (gzip) [%s]\n",
                    (FMAP_FILE_GZ_COMPRESSION == opt->input_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -J          the output is bz2 compressed (bzip2) [%s]\n",
                    (FMAP_FILE_BZ2_COMPRESSION == opt->output_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -Z          the output is gz compressed (gzip) [%s]\n",
                    (FMAP_FILE_GZ_COMPRESSION == opt->output_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -s INT      use shared memory with the following key [%d]\n", opt->shm_key);
  fmap_file_fprintf(fmap_file_stderr, "         -v          print verbose progress information\n");
  fmap_file_fprintf(fmap_file_stderr, "         -h          print this message\n");
  fmap_file_fprintf(fmap_file_stderr, "\n");

  free(reads_format);

  return 1;
}

fmap_map1_opt_t *
fmap_map1_opt_init()
{
  fmap_map1_opt_t *opt = NULL;

  opt = fmap_calloc(1, sizeof(fmap_map1_opt_t), "opt");

  // program defaults
  opt->argv = NULL;
  opt->argc = -1;
  opt->fn_fasta = opt->fn_reads = NULL;
  opt->reads_format = FMAP_READS_FORMAT_UNKNOWN;
  opt->seed_length = 32; // move this to a define block
  opt->seed_max_mm = 3; // move this to a define block
  opt->max_mm = -1; opt->max_mm_frac = 0.02; // TODO: move this to a define block 
  opt->max_gapo = -1; opt->max_gapo_frac = 0.01; // TODO: move this to a define block
  opt->max_gape = -1; opt->max_gape_frac = 0.025; // TODO: move this to a define block
  opt->pen_mm = FMAP_MAP_UTIL_PEN_MM; 
  opt->pen_gapo = FMAP_MAP_UTIL_PEN_GAPO;
  opt->pen_gape = FMAP_MAP_UTIL_PEN_GAPE;
  opt->max_cals_del = 10; // TODO: move this to a define block
  opt->indel_ends_bound = 5; // TODO: move this to a define block
  opt->max_best_cals = 32; // TODO: move this to a define block
  opt->reads_queue_size = 65536; // TODO: move this to a define block
  opt->max_entries= 2000000; // TODO: move this to a define block
  opt->num_threads = 1;
  opt->aln_output_mode = FMAP_MAP1_ALN_OUTPUT_MODE_BEST_RAND; // TODO: move this to a define block
  opt->sam_rg = NULL;
  opt->input_compr = FMAP_FILE_NO_COMPRESSION;
  opt->output_compr = FMAP_FILE_NO_COMPRESSION;
  opt->shm_key = 0;

  return opt;
}

void
fmap_map1_opt_destroy(fmap_map1_opt_t *opt)
{
  free(opt->fn_fasta);
  free(opt->fn_reads);
  free(opt->sam_rg);
  free(opt);
}

int32_t 
fmap_map1_opt_parse(int argc, char *argv[], fmap_map1_opt_t *opt)
{
  int c;

  opt->argc = argc; opt->argv = argv;

  while((c = getopt(argc, argv, "f:r:F:l:k:m:o:e:M:O:E:d:i:b:Q:q:n:a:R:jzJZs:vh")) >= 0) {
      switch(c) {
        case 'f':
          opt->fn_fasta = fmap_strdup(optarg); break;
        case 'r':
          opt->fn_reads = fmap_strdup(optarg); 
          fmap_get_reads_file_format_from_fn_int(opt->fn_reads, &opt->reads_format, &opt->input_compr);
          break;
        case 'F':
          opt->reads_format = fmap_get_reads_file_format_int(optarg); break;
        case 'l':
          opt->seed_length = atoi(optarg); break;
        case 'k':
          opt->seed_max_mm = atoi(optarg); break;
        case 'm':
          if(NULL != strstr(optarg, ".")) opt->max_mm = -1, opt->max_mm_frac = atof(optarg);
          else opt->max_mm = atoi(optarg), opt->max_mm_frac = -1.0;
          break;
        case 'o':
          if(NULL != strstr(optarg, ".")) opt->max_gapo = -1, opt->max_gapo_frac = atof(optarg);
          else opt->max_gapo = atoi(optarg), opt->max_gapo_frac = -1.0;
          break;
        case 'e':
          if(NULL != strstr(optarg, ".")) opt->max_gape = -1, opt->max_gape_frac = atof(optarg);
          else opt->max_gape = atoi(optarg), opt->max_gape_frac = -1.0;
          break;
        case 'M':
          opt->pen_mm = atoi(optarg); break;
        case 'O':
          opt->pen_gapo = atoi(optarg); break;
        case 'E':
          opt->pen_gape = atoi(optarg); break;
        case 'd':
          opt->max_cals_del = atoi(optarg); break;
        case 'i':
          opt->indel_ends_bound = atoi(optarg); break;
        case 'b':
          opt->max_best_cals = atoi(optarg); break;
        case 'Q': 
          opt->max_entries = atoi(optarg); break;
        case 'q': 
          opt->reads_queue_size = atoi(optarg); break;
        case 'n':
          opt->num_threads = atoi(optarg); break;
        case 'a':
          opt->aln_output_mode = atoi(optarg); break;
        case 'R':
          opt->sam_rg = fmap_strdup(optarg); break;
        case 'j':
          opt->input_compr = FMAP_FILE_BZ2_COMPRESSION; 
          fmap_get_reads_file_format_from_fn_int(opt->fn_reads, &opt->reads_format, &opt->input_compr);
          break;
        case 'z':
          opt->input_compr = FMAP_FILE_GZ_COMPRESSION; 
          fmap_get_reads_file_format_from_fn_int(opt->fn_reads, &opt->reads_format, &opt->input_compr);
          break;
        case 'J':
          opt->output_compr = FMAP_FILE_BZ2_COMPRESSION; break;
        case 'Z':
          opt->output_compr = FMAP_FILE_GZ_COMPRESSION; break;
        case 's':
          opt->shm_key = atoi(optarg); break;
        case 'v':
          fmap_progress_set_verbosity(1); break;
        case 'h':
        default:
          return 0;
      }
  }
  return 1;
}

void 
fmap_map1_opt_check(fmap_map1_opt_t *opt)
{
  if(NULL == opt->fn_fasta && 0 == opt->shm_key) {
      fmap_error("option -f or option -s must be specified", Exit, CommandLineArgument);
  }
  else if(NULL != opt->fn_fasta && 0 < opt->shm_key) {
      fmap_error("option -f and option -s may not be specified together", Exit, CommandLineArgument);
  }
  if(NULL == opt->fn_reads && FMAP_READS_FORMAT_UNKNOWN == opt->reads_format) {
      fmap_error("option -F or option -r must be specified", Exit, CommandLineArgument);
  }
  if(FMAP_READS_FORMAT_UNKNOWN == opt->reads_format) {
      fmap_error("the reads format (-r) was unrecognized", Exit, CommandLineArgument);
  }
  if(-1 != opt->seed_length) fmap_error_cmd_check_int(opt->seed_length, 1, INT32_MAX, "-l");

  // this will take care of the case where they are both < 0
  fmap_error_cmd_check_int((opt->max_mm_frac < 0) ? opt->max_mm : (int32_t)opt->max_mm_frac, 0, INT32_MAX, "-m"); 
  // this will take care of the case where they are both < 0
  fmap_error_cmd_check_int((opt->max_gapo_frac < 0) ? opt->max_gapo : (int32_t)opt->max_gapo_frac, 0, INT32_MAX, "-m"); 
  // this will take care of the case where they are both < 0
  fmap_error_cmd_check_int((opt->max_gape_frac < 0) ? opt->max_gape : (int32_t)opt->max_gape_frac, 0, INT32_MAX, "-m"); 

  fmap_error_cmd_check_int(opt->pen_mm, 0, INT32_MAX, "-M");
  fmap_error_cmd_check_int(opt->pen_gapo, 0, INT32_MAX, "-O");
  fmap_error_cmd_check_int(opt->pen_gape, 0, INT32_MAX, "-E");
  fmap_error_cmd_check_int(opt->max_cals_del, 1, INT32_MAX, "-d");
  fmap_error_cmd_check_int(opt->indel_ends_bound, 0, INT32_MAX, "-i");
  fmap_error_cmd_check_int(opt->max_best_cals, 0, INT32_MAX, "-b");
  fmap_error_cmd_check_int(opt->max_entries, 1, INT32_MAX, "-Q");
  if(-1 != opt->reads_queue_size) fmap_error_cmd_check_int(opt->reads_queue_size, 1, INT32_MAX, "-q");
  fmap_error_cmd_check_int(opt->num_threads, 1, INT32_MAX, "-n");
  fmap_error_cmd_check_int(opt->aln_output_mode, 0, 3, "-a");

  if(FMAP_FILE_BZ2_COMPRESSION == opt->output_compr 
     && -1 == opt->reads_queue_size) {
      fmap_error("cannot buffer reads with bzip2 output (options \"-q 1 -J\")", Exit, OutOfRange);
  }
}

int 
fmap_map1_main(int argc, char *argv[])
{
  fmap_map1_opt_t *opt = NULL;

  // random seed
  srand48(0); 

  // init opt
  opt = fmap_map1_opt_init();

  // get options
  if(1 != fmap_map1_opt_parse(argc, argv, opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return fmap_map1_usage(opt);
  }
  else { 
      // check command line arguments
      fmap_map1_opt_check(opt);
  }

  // run map1
  fmap_map1_core(opt);

  // destroy opt
  fmap_map1_opt_destroy(opt);

  fmap_progress_print2("terminating successfully");

  return 0;
}
