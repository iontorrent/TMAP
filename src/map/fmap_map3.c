#include <stdlib.h>
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
#include "../index/fmap_bwt_gen.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_bwt_match.h"
#include "../index/fmap_sa.h"
#include "../io/fmap_seq_io.h"
#include "../server/fmap_shm.h"
#include "fmap_map_util.h"
#include "fmap_map3_aux.h"
#include "fmap_map3.h"

#ifdef HAVE_LIBPTHREAD
static pthread_mutex_t fmap_map3_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t fmap_map3_read_lock_low = 0;
#define FMAP_MAP3_THREAD_BLOCK_SIZE 1024
#endif

static int32_t
fmap_map3_get_seed_length(uint64_t ref_len)
{
  int32_t k = 0;
  while(0 < ref_len) {
      ref_len >>= 2; // divide by four
      k++;
  }
  // add two just to be sure
  return k+2;
}

static inline void
fmap_map3_aln_filter(fmap_seq_t *seq, fmap_map3_aln_t *aln, 
                     int32_t score_thr, int32_t score_match, int32_t aln_output_mode)
{
  int32_t i, j;
  int32_t n_best = 0;
  int32_t best_score, cur_score, best_subo;
  int32_t n_seeds = 0, tot_seeds = 0;
  int32_t mapq;

  // estimate mapping quality TODO: this needs to be refined
  best_score = INT32_MIN;
  best_subo = INT32_MIN;
  n_best = 0;
  for(i=0;i<aln->n;i++) {
      cur_score = aln->hits[i].score;
      tot_seeds += aln->hits[i].n_seeds;
      if(best_score < cur_score) {
          best_score = cur_score;
          n_best = 1;
          n_seeds = aln->hits[i].n_seeds;
      }
      else if(!(cur_score < best_score)) { // qual
          n_best++;
      }
      else {
          if(best_subo < cur_score) {
              best_subo = cur_score;
          }
          cur_score = aln->hits[i].score_subo;
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
  for(i=0;i<aln->n;i++) {
      aln->hits[i].mapq = mapq;
  }

  if(FMAP_MAP_UTIL_ALN_MODE_ALL == aln_output_mode
     || aln->n <= 1) {
      return;
  }

  best_score = INT32_MIN;
  n_best = 0;
  for(i=0;i<aln->n;i++) {
      cur_score = aln->hits[i].score;
      if(best_score < cur_score) {
          best_score = cur_score;
          n_best = 1;
      }
      else if(!(cur_score < best_score)) { // equal
          n_best++;
      }
  }

  // copy to the front
  if(n_best < aln->n) {
      for(i=j=0;i<aln->n;i++) {
          cur_score = aln->hits[i].score;
          if(cur_score < best_score) { // not the best
              free(aln->hits[i].cigar);
              aln->hits[i].cigar = NULL;
              aln->hits[i].n_cigar = 0;
          }
          else {
              if(j < i) { // copy if we are not on the same index
                  aln->hits[j] = aln->hits[i];
                  aln->hits[i].cigar = NULL;
              }
              j++;
          }
      }
      // reallocate
      fmap_map3_aln_realloc(aln, n_best);
  }

  if(FMAP_MAP_UTIL_ALN_MODE_UNIQ_BEST == aln_output_mode) {
      if(1 < n_best) { // there can only be one
          fmap_map3_aln_realloc(aln, 0);
      }
  }
  else if(FMAP_MAP_UTIL_ALN_MODE_RAND_BEST == aln_output_mode) { // get a random
      i = (int32_t)(drand48() * aln->n);
      if(0 != i) {
          free(aln->hits[0].cigar);
          aln->hits[0] = aln->hits[i];
          aln->hits[i].cigar = NULL;
      }
      // reallocate
      fmap_map3_aln_realloc(aln, 1);
  }
  else if(FMAP_MAP_UTIL_ALN_MODE_ALL_BEST == aln_output_mode) {
      // do nothing
  }
  else {
      fmap_error("bug encountered", Exit, OutOfRange);
  }
}

static inline void
fmap_map3_print_sam(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_map3_hit_t *hit)
{
  fmap_sam_print_mapped(fmap_file_stdout, seq, refseq,
                        hit->strand, hit->seqid, hit->pos,
                        hit->mapq, hit->cigar, hit->n_cigar,
                        "\tAS:i:%d\tXS:i:%d", 
                        hit->score, hit->n_seeds);
}


static void
fmap_map3_core_worker(fmap_seq_t **seq_buffer, fmap_map3_aln_t **alns, int32_t seq_buffer_length, 
                      fmap_refseq_t *refseq, fmap_bwt_t *bwt, fmap_sa_t *sa,
                      int32_t tid, fmap_map3_opt_t *opt)
{
  int32_t low = 0, high;

  while(low < seq_buffer_length) {
#ifdef HAVE_LIBPTHREAD
      if(1 < opt->num_threads) {
          pthread_mutex_lock(&fmap_map3_read_lock);

          // update bounds
          low = fmap_map3_read_lock_low;
          fmap_map3_read_lock_low += FMAP_MAP3_THREAD_BLOCK_SIZE;
          high = low + FMAP_MAP3_THREAD_BLOCK_SIZE;
          if(seq_buffer_length < high) {
              high = seq_buffer_length; 
          }

          pthread_mutex_unlock(&fmap_map3_read_lock);
      }
      else {
          high = seq_buffer_length; // process all
      }
#else 
      high = seq_buffer_length; // process all
#endif
      while(low<high) {
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

          alns[low] = fmap_map3_aux_core(seq, refseq, bwt, sa, opt);

          // filter the alignments
          fmap_map3_aln_filter(seq_buffer[low], alns[low],
                               opt->score_thr, opt->score_match, opt->aln_output_mode);

          // destroy
          fmap_seq_destroy(seq[0]);
          fmap_seq_destroy(seq[1]);

          // next
          low++;
      }
  }
}

static void *
fmap_map3_core_thread_worker(void *arg)
{
  fmap_map3_thread_data_t *thread_data = (fmap_map3_thread_data_t*)arg;

  fmap_map3_core_worker(thread_data->seq_buffer, thread_data->alns, thread_data->seq_buffer_length, 
                        thread_data->refseq, thread_data->bwt, thread_data->sa, 
                        thread_data->tid, thread_data->opt);

  return arg;
}

static void 
fmap_map3_core(fmap_map3_opt_t *opt)
{
  uint32_t i, j, n_reads_processed=0;
  int32_t seq_buffer_length;
  fmap_refseq_t *refseq=NULL;
  fmap_bwt_t *bwt=NULL;
  fmap_sa_t *sa=NULL;
  fmap_file_t *fp_reads=NULL;
  fmap_seq_io_t *seqio = NULL;
  fmap_seq_t **seq_buffer = NULL;
  fmap_map3_aln_t **alns = NULL;
  fmap_shm_t *shm = NULL;
  int32_t reads_queue_size;

  // adjust opt for opt->score_match
  opt->score_thr *= opt->score_match;

  // For suffix search we need the reverse bwt/sa and forward refseq
  if(0 == opt->shm_key) {
      fmap_progress_print("reading in reference data");
      refseq = fmap_refseq_read(opt->fn_fasta, 0);
      bwt = fmap_bwt_read(opt->fn_fasta, 1);
      sa = fmap_sa_read(opt->fn_fasta, 1);
      fmap_progress_print2("reference data read in");
  }
  else {
      fmap_progress_print("retrieving reference data from shared memory");
      shm = fmap_shm_init(opt->shm_key, 0, 0);
      if(NULL == (refseq = fmap_refseq_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_REFSEQ)))) {
          fmap_error("the packed reference sequence was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (bwt = fmap_bwt_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_REV_BWT)))) {
          fmap_error("the reverse BWT string was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (sa = fmap_sa_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_REV_SA)))) {
          fmap_error("the reverse SA was not found in shared memory", Exit, SharedMemoryListing);
      }
      fmap_progress_print2("reference data retrieved from shared memory");
  }

  // Set the seed length
  if(-1 == opt->seed_length) {
      opt->seed_length = fmap_map3_get_seed_length(refseq->len);
      fmap_progress_print("setting the seed length to %d", opt->seed_length);
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
  alns = fmap_malloc(sizeof(fmap_map3_aln_t*)*reads_queue_size, "alns");

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
          fmap_map3_core_worker(seq_buffer, alns, seq_buffer_length, refseq, bwt, sa, 0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          fmap_map3_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = fmap_calloc(opt->num_threads, sizeof(pthread_t), "threads");
          thread_data = fmap_calloc(opt->num_threads, sizeof(fmap_map3_thread_data_t), "thread_data");
          fmap_map3_read_lock_low = 0; // ALWAYS set before running threads 

          for(i=0;i<opt->num_threads;i++) {
              thread_data[i].seq_buffer = seq_buffer;
              thread_data[i].seq_buffer_length = seq_buffer_length;
              thread_data[i].alns = alns;
              thread_data[i].refseq = refseq;
              thread_data[i].bwt = bwt;
              thread_data[i].sa = sa;;
              thread_data[i].tid = i;
              thread_data[i].opt = opt; 
              if(0 != pthread_create(&threads[i], &attr, fmap_map3_core_thread_worker, &thread_data[i])) {
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
      fmap_map3_core_worker(seq_buffer, alns, seq_buffer_length, refseq, bwt, sa, 0, opt);
#endif

      fmap_progress_print("writing alignments");
      for(i=0;i<seq_buffer_length;i++) {

          if(0 < alns[i]->n) {
              for(j=0;j<alns[i]->n;j++) {
                  fmap_map3_print_sam(seq_buffer[i], refseq, &alns[i]->hits[j]);
              }
          }
          else {
              fmap_sam_print_unmapped(fmap_file_stdout, seq_buffer[i]);
          }

          // free alignments
          fmap_map3_aln_destroy(alns[i]);
          alns[i] = NULL;
      }

      if(-1 == opt->reads_queue_size) {
          fmap_file_fflush(fmap_file_stdout, 1);
      }

      n_reads_processed += seq_buffer_length;
      fmap_progress_print2("processed %d reads", n_reads_processed);
  }

  // close the input/output
  fmap_file_fclose(fmap_file_stdout);
  fmap_file_fclose(fp_reads);

  // free memory
  for(i=0;i<reads_queue_size;i++) {
      fmap_seq_destroy(seq_buffer[i]);
  }
  free(seq_buffer);
  free(alns);
  fmap_refseq_destroy(refseq);
  fmap_bwt_destroy(bwt);
  fmap_sa_destroy(sa);
  fmap_seq_io_destroy(seqio);
  if(0 < opt->shm_key) {
      fmap_shm_destroy(shm, 0);
  }
}

static int 
usage(fmap_map3_opt_t *opt)
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
  fmap_file_fprintf(fmap_file_stderr, "Usage: %s map3 [options]", PACKAGE);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (optional):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -f FILE     the FASTA reference file name [%s]\n", opt->fn_fasta);
  fmap_file_fprintf(fmap_file_stderr, "         -r FILE     the reads file name [%s]\n", (NULL == opt->fn_reads) ? "stdin" : opt->fn_reads);
  fmap_file_fprintf(fmap_file_stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff) [%s]\n", reads_format);
  fmap_file_fprintf(fmap_file_stderr, "         -l INT      the k-mer length to seed CALs (-1 tunes to the genome size) [%d]\n", opt->seed_length);
  fmap_file_fprintf(fmap_file_stderr, "         -S INT      the maximum number of hits returned by a seed [%d]\n", opt->max_seed_hits);
  fmap_file_fprintf(fmap_file_stderr, "         -b INT      the band width to group seeds [%d]\n", opt->max_seed_band);
  fmap_file_fprintf(fmap_file_stderr, "         -w INT      the extra bases to add before and after the target during Smith-Waterman [%d]\n", opt->sw_offset);
  fmap_file_fprintf(fmap_file_stderr, "         -A INT      the match score [%d]\n", opt->score_match); 
  fmap_file_fprintf(fmap_file_stderr, "         -M INT      the mismatch penalty [%d]\n", opt->pen_mm); 
  fmap_file_fprintf(fmap_file_stderr, "         -O INT      the indel start penalty [%d]\n", opt->pen_gapo); 
  fmap_file_fprintf(fmap_file_stderr, "         -E INT      the indel extend penalty [%d]\n", opt->pen_gape); 
  fmap_file_fprintf(fmap_file_stderr, "         -T INT      score threshold divided by the match score [%d]\n", opt->score_thr);
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

  free(reads_format);

  return 1;
}

static fmap_map3_opt_t *
fmap_map3_opt_init()
{
  fmap_map3_opt_t *opt = NULL;

  opt = fmap_calloc(1, sizeof(fmap_map3_opt_t), "opt");

  // program defaults
  opt->argv = NULL;
  opt->argc = -1;
  opt->fn_fasta = opt->fn_reads = NULL;
  opt->reads_format = FMAP_READS_FORMAT_UNKNOWN;
  opt->seed_length = -1; // move this to a define block
  opt->max_seed_hits = 8; // move this to a define block
  opt->max_seed_band = 50; // move this to a define block
  opt->sw_offset = 10; // move this to a define block
  opt->score_match = 1;
  opt->pen_mm = 3; opt->pen_gapo = 5; opt->pen_gape = 2; // TODO: move this to a define block
  opt->score_thr = 20;
  opt->aln_global = 0;
  opt->hp_diff = 0;
  opt->reads_queue_size = 65536; // TODO: move this to a define block
  opt->num_threads = 1;
  opt->aln_output_mode = FMAP_MAP_UTIL_ALN_MODE_RAND_BEST;
  opt->input_compr = FMAP_FILE_NO_COMPRESSION;
  opt->output_compr = FMAP_FILE_NO_COMPRESSION;
  opt->shm_key = 0;

  return opt;
}

static void
fmap_map3_opt_destroy(fmap_map3_opt_t *opt)
{
  free(opt->fn_fasta);
  free(opt->fn_reads);
  free(opt);
}

int 
fmap_map3_main(int argc, char *argv[])
{
  int c;
  fmap_map3_opt_t *opt = NULL;

  srand48(0); // random seed
  opt = fmap_map3_opt_init(argv, argc);
  opt->argc = argc; opt->argv = argv;

  while((c = getopt(argc, argv, "f:r:F:l:S:b:w:A:M:O:E:T:gH:q:n:a:jzJZs:vh")) >= 0) {
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
        case 'S':
          opt->max_seed_hits = atoi(optarg); break;
        case 'b':
          opt->max_seed_band = atoi(optarg); break;
        case 'w':
          opt->sw_offset = atoi(optarg); break;
        case 'A':
          opt->score_match = atoi(optarg); break;
        case 'M':
          opt->pen_mm = atoi(optarg); break;
        case 'O':
          opt->pen_gapo = atoi(optarg); break;
        case 'E':
          opt->pen_gape = atoi(optarg); break;
        case 'T':
          opt->score_thr = atoi(optarg); break;
        case 'g':
          opt->aln_global = 1; break;
        case 'q': 
          opt->reads_queue_size = atoi(optarg); break;
        case 'H':
          opt->hp_diff = atoi(optarg); break;
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
          return usage(opt);
      }
  }

  if(argc != optind || 1 == argc) {
      return usage(opt);
  }
  else { // check command line arguments
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
      fmap_error_cmd_check_int(opt->max_seed_hits, 1, INT32_MAX, "-S");
      fmap_error_cmd_check_int(opt->max_seed_hits, 1, INT32_MAX, "-b");
      fmap_error_cmd_check_int(opt->sw_offset, 1, INT32_MAX, "-w");

      fmap_error_cmd_check_int(opt->score_match, 0, INT32_MAX, "-A");
      fmap_error_cmd_check_int(opt->pen_mm, 0, INT32_MAX, "-M");
      fmap_error_cmd_check_int(opt->pen_gapo, 0, INT32_MAX, "-O");
      fmap_error_cmd_check_int(opt->pen_gape, 0, INT32_MAX, "-E");
      fmap_error_cmd_check_int(opt->score_thr, 0, INT32_MAX, "-T");
      if(-1 != opt->reads_queue_size) fmap_error_cmd_check_int(opt->reads_queue_size, 1, INT32_MAX, "-q");
      fmap_error_cmd_check_int(opt->num_threads, 1, INT32_MAX, "-n");
      fmap_error_cmd_check_int(opt->aln_output_mode, 0, 3, "-a");
      fmap_error_cmd_check_int(opt->aln_output_mode, 0, INT32_MAX, "-H");

      if(FMAP_FILE_BZ2_COMPRESSION == opt->output_compr 
         && -1 == opt->reads_queue_size) {
          fmap_error("cannot buffer reads with bzip2 output (options \"-q 1 -J\")", Exit, OutOfRange);
      }
  }

  fmap_map3_core(opt);

  fmap_map3_opt_destroy(opt);

  fmap_progress_print2("terminating successfully");

  return 0;
}
