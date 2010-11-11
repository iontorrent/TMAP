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
#include "fmap_map_util.h"
#include "fmap_map2_mempool.h"
#include "fmap_map2_aux.h"
#include "fmap_map2_core.h"
#include "fmap_map2.h"

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t fmap_map2_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t fmap_map2_read_lock_low = 0;
#define FMAP_MAP1_THREAD_BLOCK_SIZE 1024
#endif

static void
fmap_map2_filter_sam(fmap_seq_t *seq, fmap_map2_sam_t *sam, int32_t aln_output_mode)
{
  int32_t i, j;
  int32_t n_best = 0;
  int32_t best_score, cur_score;
  if(FMAP_MAP_UTIL_ALN_MODE_ALL == aln_output_mode || sam->num_entries <= 1) {
      return;
  }

  best_score = INT32_MIN;
  n_best = 0;
  for(i=0;i<sam->num_entries;i++) {
      cur_score = sam->entries[i].AS;
      if(best_score < cur_score) {
          best_score = cur_score;
          n_best = 1;
      }
      else if(!(cur_score < best_score)) { // equal
          n_best++;
      }
  }

  // copy to the front
  if(n_best < sam->num_entries) {
      for(i=j=0;i<sam->num_entries;i++) {
          cur_score = sam->entries[i].AS;
          if(cur_score < best_score) { // not the best
              free(sam->entries[i].cigar);
              sam->entries[i].cigar = NULL;
          }
          else if(j < i) { // copy if we are not on the same index
              sam->entries[j] = sam->entries[i];
              sam->entries[i].cigar = NULL;
              j++;
          }
      }
      // reallocate
      fmap_map2_sam_realloc(sam, n_best);
  }

  if(FMAP_MAP_UTIL_ALN_MODE_UNIQ_BEST == aln_output_mode) {
      if(1 < n_best) { // there can only be one
          fmap_map2_sam_realloc(sam, 0);
      }
  } 
  else if(FMAP_MAP_UTIL_ALN_MODE_RAND_BEST == aln_output_mode) {
      // get a random
      i = drand48() * n_best;
      if(0 != i ) {
          free(sam->entries[0].cigar);
          sam->entries[0] = sam->entries[i];
          sam->entries[i].cigar = NULL;
      }
      // reallocate
      fmap_map2_sam_realloc(sam, 1);
  }
  else if(FMAP_MAP_UTIL_ALN_MODE_ALL_BEST == aln_output_mode) {
      // do nothing
  }
  else {
      fmap_error("bug encountered", Exit, OutOfRange);
  }
}

static inline void
fmap_map2_print_sam(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_map2_sam_entry_t *sam)
{

  if(0 < sam->XI) {
      fmap_sam_print_mapped(fmap_file_stdout, seq, refseq, 
                            sam->strand, sam->seqid, sam->pos, sam->mapq,
                            sam->cigar, sam->n_cigar, 
                            "\tAS:i:%d\tXS:i:%d\tXF:i:%d\tXE:i:%d\tXI:i:%d",
                            sam->AS, sam->XS, sam->XF, sam->XE, sam->XI);
  }
  else {
      fmap_sam_print_mapped(fmap_file_stdout, seq, refseq, 
                            sam->strand, sam->seqid, sam->pos, sam->mapq,
                            sam->cigar, sam->n_cigar, 
                            "\tAS:i:%d\tXS:i:%d\tXF:i:%d\tXE:i:%d",
                            sam->AS, sam->XS, sam->XF, sam->XE);
  }
}

static void
fmap_map2_core_worker(fmap_seq_t **seq_buffer, int32_t seq_buffer_length, fmap_map2_sam_t **sams,
                      fmap_refseq_t *refseq, fmap_bwt_t *bwt[2], fmap_sa_t *sa[2],
                      int32_t tid, fmap_map2_opt_t * opt)
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
          fmap_map2_read_lock_low += FMAP_MAP1_THREAD_BLOCK_SIZE;
          high = low + FMAP_MAP1_THREAD_BLOCK_SIZE;
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
          if(NULL != sams[low]) {
              fmap_map2_filter_sam(seq_buffer[low], sams[low], opt->aln_output_mode);
          }

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
fmap_map2_core(fmap_map2_opt_t *opt)
{
  uint32_t i, j, n_reads_processed=0;
  int32_t seq_buffer_length;
  double scalar;
  fmap_refseq_t *refseq = NULL;
  fmap_bwt_t *bwt[2] = {NULL, NULL};
  fmap_sa_t *sa[2] = {NULL, NULL};
  fmap_file_t *fp_reads=NULL;
  fmap_seq_io_t *seqio = NULL;
  fmap_shm_t *shm = NULL;
  fmap_seq_t **seq_buffer = NULL;
  fmap_map2_sam_t **sams = NULL;
  int32_t reads_queue_size;

  scalar = opt->score_match / log(opt->yita);
  fmap_progress_print( "mismatch: %lf, gap_open: %lf, gap_ext: %lf",
                      exp(-opt->pen_mm / scalar) / opt->yita,
                      exp(-opt->pen_gapo / scalar),
                      exp(-opt->pen_gape / scalar));

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
          fmap_error("the reverse SA was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (sa[1] = fmap_sa_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_REV_SA)))) {
          fmap_error("the reverse SA was not found in shared memory", Exit, SharedMemoryListing);
      }
      fmap_progress_print2("reference data retrieved from shared memory");

  }

  // Note: 'fmap_file_stdout' should not have been previously modified
  fmap_file_stdout = fmap_file_fdopen(fileno(stdout), "wb", opt->output_compr);

  // SAM header
  fmap_sam_print_header(fmap_file_stdout, refseq, opt->argc, opt->argv);

  // allocate the buffer
  if(-1 == opt->reads_queue_size) {
      reads_queue_size = 1;
  }
  else {
      reads_queue_size = opt->reads_queue_size;
  }
  seq_buffer = fmap_malloc(sizeof(fmap_seq_t*)*reads_queue_size, "seq_buffer");
  sams = fmap_malloc(sizeof(fmap_map2_sam_t*)*reads_queue_size, "alnseqs");

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

  fmap_progress_print("processing reads");
  while(0 < (seq_buffer_length = fmap_seq_io_read_buffer(seqio, seq_buffer, reads_queue_size))) {

      // do alignment
#ifdef HAVE_LIBPTHREAD
      if(1 == opt->num_threads) {
          fmap_map2_core_worker(seq_buffer, seq_buffer_length, sams,
                                refseq, bwt, sa, 0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          fmap_map2_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = fmap_calloc(opt->num_threads, sizeof(pthread_t), "threads");
          thread_data = fmap_calloc(opt->num_threads, sizeof(fmap_map2_thread_data_t), "thread_data");
          fmap_map2_read_lock_low = 0; // ALWAYS set before running threads 

          for(i=0;i<opt->num_threads;i++) {
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
          for(i=0;i<opt->num_threads;i++) {
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

      fmap_progress_print("writing alignments");

      for(i=0;i<seq_buffer_length;i++) {
          if(NULL != sams[i] && 0 < sams[i]->num_entries) {
              // print mapped reads
              for(j=0;j<sams[i]->num_entries;j++) {
                  fmap_map2_print_sam(seq_buffer[i], refseq, &sams[i]->entries[j]);
              }
          }
          else {
              // print unmapped reads
              fmap_sam_print_unmapped(fmap_file_stdout, seq_buffer[i]);
          }

          // free alignments
          fmap_map2_sam_destroy(sams[i]);
          sams[i] = NULL;
      }

      if(-1 == opt->reads_queue_size) {
          fmap_file_fflush(fmap_file_stdout, 1);
      }

      n_reads_processed += seq_buffer_length;
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
fmap_map2_usage(fmap_map2_opt_t *opt)
{
  char *reads_format = fmap_get_reads_file_format_string(opt->reads_format);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Usage: %s map2 [options]", PACKAGE);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (required):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -f FILE     the FASTA reference file name [%s]\n", opt->fn_fasta);
  fmap_file_fprintf(fmap_file_stderr, "         -r FILE     the reads file name [%s]\n", (NULL == opt->fn_reads) ? "stdin" : opt->fn_reads);
  fmap_file_fprintf(fmap_file_stderr, "Options (optional):\n");
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff) [%s]\n", reads_format);
  fmap_file_fprintf(fmap_file_stderr, "         -A INT      score for a match [%d]\n", opt->score_match);
  fmap_file_fprintf(fmap_file_stderr, "         -M INT      mismatch penalty [%d]\n", opt->pen_mm);
  fmap_file_fprintf(fmap_file_stderr, "         -O INT      gap open penalty [%d]\n", opt->pen_gapo);
  fmap_file_fprintf(fmap_file_stderr, "         -E INT      gap extension penalty [%d]\n", opt->pen_gape);
  //fmap_file_fprintf(fmap_file_stderr, "         -y FLOAT    error recurrence coef. (4..16) [%.1lf]\n", opt->yita);
  fmap_file_fprintf(fmap_file_stderr, "         -m FLOAT    mask level [%.2f]\n", opt->mask_level);
  fmap_file_fprintf(fmap_file_stderr, "         -c FLOAT    coefficient of length-threshold adjustment [%.1lf]\n", opt->length_coef);
  fmap_file_fprintf(fmap_file_stderr, "         -w INT      band width [%d]\n", opt->band_width);
  fmap_file_fprintf(fmap_file_stderr, "         -T INT      score threshold divided by the match score [%d]\n", opt->score_thr);
  fmap_file_fprintf(fmap_file_stderr, "         -S INT      maximum seeding interval size [%d]\n", opt->max_seed_intv);
  fmap_file_fprintf(fmap_file_stderr, "         -b INT      Z-best [%d]\n", opt->z_best);
  fmap_file_fprintf(fmap_file_stderr, "         -N INT      # seeds to trigger reverse alignment [%d]\n", opt->seeds_rev);
  fmap_file_fprintf(fmap_file_stderr, "         -g          align the full read (global alignment) [%s]\n", (0 == opt->aln_global) ? "false" : "true");
  fmap_file_fprintf(fmap_file_stderr, "         -q INT      the queue size for the reads (-1 disables) [%d]\n", opt->reads_queue_size);
  fmap_file_fprintf(fmap_file_stderr, "         -n INT      the number of threads [%d]\n", opt->num_threads);
  fmap_file_fprintf(fmap_file_stderr, "         -a INT      output filter [%d]\n", opt->aln_output_mode);
  fmap_file_fprintf(fmap_file_stderr, "                             0 - unique best hits\n");
  fmap_file_fprintf(fmap_file_stderr, "                             1 - random best hit\n");
  fmap_file_fprintf(fmap_file_stderr, "                             2 - all best hits\n");
  fmap_file_fprintf(fmap_file_stderr, "                             3 - all alignments\n");
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

  return 1;
}

fmap_map2_opt_t *
fmap_map2_opt_init()
{
  fmap_map2_opt_t *opt = NULL;

  opt = fmap_calloc(1, sizeof(fmap_map2_opt_t), "opt");

  opt->fn_fasta = opt->fn_reads = NULL;
  opt->reads_format = FMAP_READS_FORMAT_UNKNOWN;
  opt->score_match = 1; opt->pen_mm = 3; opt->pen_gapo = 5; opt->pen_gape = 2;
  opt->yita = 5.5f; opt->mask_level = 0.50; opt->length_coef = 5.5f;
  opt->band_width = 50; opt->score_thr = 30;
  opt->max_seed_intv = 3; opt->z_best = 5; opt->seeds_rev = 5;
  opt->aln_global = 0;
  opt->reads_queue_size = 65536;
  opt->num_threads = 1;
  opt->aln_output_mode = FMAP_MAP_UTIL_ALN_MODE_RAND_BEST; 
  opt->input_compr = FMAP_FILE_NO_COMPRESSION;
  opt->output_compr = FMAP_FILE_NO_COMPRESSION;
  opt->shm_key = 0;

  return opt;
}

void
fmap_map2_opt_destroy(fmap_map2_opt_t *opt)
{
  free(opt->fn_fasta);
  free(opt->fn_reads);
  free(opt);
}

int 
fmap_map2_opt_parse(int argc, char *argv[], fmap_map2_opt_t *opt)
{
  int c;

  opt->argc = argc; opt->argv = argv;

  while((c = getopt(argc, argv, "f:r:F:A:M:O:E:y:m:c:w:T:S:b:N:gq:n:a:jzJZs:vh")) >= 0) {
      switch (c) {
        case 'f':
          opt->fn_fasta = fmap_strdup(optarg); break;
        case 'r':
          opt->fn_reads = fmap_strdup(optarg);
          fmap_get_reads_file_format_from_fn_int(opt->fn_reads, &opt->reads_format, &opt->input_compr);
          break;
        case 'F':
          opt->reads_format = fmap_get_reads_file_format_int(optarg); break;
        case 'A':
          opt->score_match = atoi(optarg); break;
        case 'M':
          opt->pen_mm = atoi(optarg); break;
        case 'O':
          opt->pen_gapo = atoi(optarg); break;
        case 'E':
          opt->pen_gape = atoi(optarg); break;
          /*
             case 'y': 
             opt->yita = atof(optarg); break;
             */
        case 'm': 
          opt->mask_level = atof(optarg); break;
        case 'c': 
          opt->length_coef = atof(optarg); break;
        case 'w': 
          opt->band_width = atoi(optarg); break;
        case 'T': 
          opt->score_thr = atoi(optarg); break;
        case 'S':
          opt->max_seed_intv = atoi(optarg); break;
        case 'b': 
          opt->z_best= atoi(optarg); break;
        case 'N':
          opt->seeds_rev = atoi(optarg); break;
        case 'g':
          opt->aln_global = 1; break;
        case 'q':
          opt->reads_queue_size = atoi(optarg); break;
        case 'n':
          opt->num_threads = atoi(optarg); break;
        case 'a':
          opt->aln_output_mode = atoi(optarg); break;
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
fmap_map2_opt_check(fmap_map2_opt_t *opt)
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

  fmap_error_cmd_check_int(opt->score_match, 0, INT32_MAX, "-A");
  fmap_error_cmd_check_int(opt->pen_mm, 0, INT32_MAX, "-M");
  fmap_error_cmd_check_int(opt->pen_gapo, 0, INT32_MAX, "-O");
  fmap_error_cmd_check_int(opt->pen_gape, 0, INT32_MAX, "-E");
  //fmap_error_cmd_check_int(opt->yita, 0, 1, "-y");
  fmap_error_cmd_check_int(opt->mask_level, 0, 1, "-m");
  fmap_error_cmd_check_int(opt->length_coef, 0, INT32_MAX, "-c");
  fmap_error_cmd_check_int(opt->band_width, 0, INT32_MAX, "-w");
  fmap_error_cmd_check_int(opt->score_thr, 0, INT32_MAX, "-T");
  fmap_error_cmd_check_int(opt->max_seed_intv, 0, INT32_MAX, "-S");
  fmap_error_cmd_check_int(opt->z_best, 1, INT32_MAX, "-Z");
  fmap_error_cmd_check_int(opt->seeds_rev, 0, INT32_MAX, "-N");
  if(-1 != opt->reads_queue_size) fmap_error_cmd_check_int(opt->reads_queue_size, 1, INT32_MAX, "-q");
  fmap_error_cmd_check_int(opt->num_threads, 1, INT32_MAX, "-n");

  if(FMAP_FILE_BZ2_COMPRESSION == opt->output_compr
     && -1 == opt->reads_queue_size) {
      fmap_error("cannot buffer reads with bzip2 output (options \"-q 1 -J\")", Exit, OutOfRange);
  }   
}

int 
fmap_map2_main(int argc, char *argv[])
{
  fmap_map2_opt_t *opt = NULL;

  // random seed
  srand48(0);

  // init opt
  opt = fmap_map2_opt_init();

  // get options
  if(1 != fmap_map2_opt_parse(argc, argv, opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return fmap_map2_usage(opt);
  }
  else { 
      // check command line arguments
      fmap_map2_opt_check(opt);
  }

  // run map2
  fmap_map2_core(opt);

  // destroy opt
  fmap_map2_opt_destroy(opt);

  fmap_progress_print2("terminating successfully");

  return 0;
}
