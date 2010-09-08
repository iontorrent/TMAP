#include <stdlib.h>
#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt_gen.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_bwt_match.h"
#include "../index/fmap_sa.h"
#include "../io/fmap_seq.h"
#include "fmap_map1.h"

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
static pthread_mutex_t fmap_map1_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t fmap_map1_read_lock_tid = 0;
#define FMAP_MAP1_THREAD_BLOCK_SIZE 1024
#endif

static void
fmap_map1_core_worker(fmap_seq_t **seq_buffer, int32_t seq_buffer_length,
                      fmap_refseq_t *refseq, fmap_bwt_t *bwt, fmap_sa_t *sa,
                      int32_t tid, fmap_map1_opt_t *opt)
{
  int32_t low = 0, high, i;

  while(low < seq_buffer_length) {
#ifdef HAVE_LIBPTHREAD
      if(1 < opt->num_threads) {
          pthread_mutex_lock(&fmap_map1_read_lock);

          // update bounds
          low = fmap_map1_read_lock_tid;
          fmap_map1_read_lock_tid += FMAP_MAP1_THREAD_BLOCK_SIZE;
          high = low + FMAP_MAP1_THREAD_BLOCK_SIZE;

          pthread_mutex_unlock(&fmap_map1_read_lock);
      }
      else {
          high = seq_buffer_length; // process all
      }
#else 
      high = seq_buffer_length; // process all
#endif
      for(i=low;i<high;i++) { // process each read
          // TODO
      }
  }
}
              
static void *
fmap_map1_core_thread_worker(void *arg)
{
  fmap_map1_thread_data_t *thread_data = (fmap_map1_thread_data_t*)arg;

  fmap_map1_core_worker(thread_data->seq_buffer, thread_data->seq_buffer_length,
                        thread_data->refseq, thread_data->bwt, thread_data->sa,
                        thread_data->tid, thread_data->opt);

  return arg;
}

static void 
fmap_map1_core(fmap_map1_opt_t *opt)
{
  uint32_t i;
  int32_t seq_buffer_length;
  fmap_refseq_t *refseq=NULL;
  fmap_bwt_t *bwt=NULL;
  fmap_sa_t *sa=NULL;
  fmap_file_t *fp_reads=NULL;
  fmap_seq_io_t *seqio = NULL;
  fmap_seq_t **seq_buffer = NULL;

  // SAM header
  /*
     refseq = fmap_refseq_read(opt->fn_fasta, 0);
     for(i=0;i<refseq->num_annos;i++) {
     fprintf(stdout, "@SQ\tSN:%s\tLN:%d\n",
     refseq->annos[i].name, (int)refseq->annos[i].len);
     }
     */

  // TODO modify so we use forward search
  bwt = fmap_bwt_read(opt->fn_fasta, 1);
  sa = fmap_sa_read(opt->fn_fasta, 1);

  // TODO: support FASTA/FASTQ/SFF
  fp_reads = fmap_file_fopen(opt->fn_reads, "rb", FMAP_FILE_NO_COMPRESSION);
  seqio = fmap_seq_io_init(fp_reads);

  // initialize the buffer
  seq_buffer = fmap_malloc(sizeof(fmap_seq_t*), "seq_buffer");
  for(i=0;i<opt->reads_queue_size;i++) {
      seq_buffer[i] = fmap_seq_init();
  }

  while(0 < (seq_buffer_length = fmap_seq_io_read_buffer(seqio, seq_buffer, opt->reads_queue_size))) {
#ifdef HAVE_LIBPTHREAD
      if(1 == opt->num_threads) {
          fmap_map1_core_worker(seq_buffer, seq_buffer_length, refseq, bwt, sa, 0, opt);
      }
      else {
          pthread_attr_t attr;
          pthread_t *threads = NULL;
          fmap_map1_thread_data_t *thread_data=NULL;

          pthread_attr_init(&attr);
          pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

          threads = fmap_calloc(opt->num_threads, sizeof(pthread_t), "threads");
          thread_data = fmap_calloc(opt->num_threads, sizeof(fmap_map1_thread_data_t), "thread_data");
          fmap_map1_read_lock_tid = 0; // ALWAYS set before running threads 

          for(i=0;i<opt->num_threads;i++) {
              thread_data[i].seq_buffer = seq_buffer;
              thread_data[i].seq_buffer_length = seq_buffer_length;
              thread_data[i].refseq = refseq;
              thread_data[i].bwt = bwt;
              thread_data[i].sa = sa;
              thread_data[i].tid = i;
              thread_data[i].opt = opt; 
              if(0 != pthread_create(&threads[i], &attr, fmap_map1_core_thread_worker, &thread_data[i])) {
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
      fmap_map1_core_worker(seq_buffer, seq_buffer_length, refseq, bwt, sa, 0, opt);
#endif
  }

  // free memory
  for(i=0;i<opt->reads_queue_size;i++) {
      fmap_seq_destroy(seq_buffer[i]);
  }
  free(seq_buffer);
  fmap_file_fclose(fp_reads);
  fmap_refseq_destroy(refseq);
  fmap_bwt_destroy(bwt);
  fmap_sa_destroy(sa);
}

static int 
usage(fmap_map1_opt_t *opt)
{
  char *reads_format = fmap_get_reads_file_format_string(opt->reads_format);

  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: %s map1 [options]", PACKAGE);
  fprintf(stderr, "\n");
  fprintf(stderr, "Options (required):\n");
  fprintf(stderr, "         -f FILE     the FASTA reference file name [%s]\n", opt->fn_fasta);
  fprintf(stderr, "         -r FILE     the reads file name [%s]\n", opt->fn_reads);
  fprintf(stderr, "Options (optional):\n");
  fprintf(stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff) [%s]\n", reads_format);
  
  fprintf(stderr, "         -m NUM      maximum number of or (read length) fraction of mismatches");
  if(opt->max_mm < 0) fprintf(stderr, " [fraction: %lf]\n", opt->max_mm_frac);
  else fprintf(stderr, " [number: %d]\n", opt->max_mm); 
 
  fprintf(stderr, "         -o NUM      maximum number of or (read length) fraction of indel starts");
  if(opt->max_gapo < 0) fprintf(stderr, " [fraction: %lf]\n", opt->max_gapo_frac);
  else fprintf(stderr, " [number: %d]\n", opt->max_gapo); 

  fprintf(stderr, "         -e NUM      maximum number of or (read length) fraction of indel extensions");
  if(opt->max_gape < 0) fprintf(stderr, " [fraction: %lf]\n", opt->max_gape_frac);
  else fprintf(stderr, " [number: %d]\n", opt->max_gape); 

  fprintf(stderr, "         -M INT      the mismatch penalty [%d]\n", opt->pen_mm); 
  fprintf(stderr, "         -O INT      the indel start penalty [%d]\n", opt->pen_gapo); 
  fprintf(stderr, "         -E INT      the indel extend penalty [%d]\n", opt->pen_gape); 
  fprintf(stderr, "         -d INT      the maximum number of CALs to extend a deletion [%d]\n", opt->max_cals_del); 
  fprintf(stderr, "         -i INT      indels are not allowed within INT number of bps from the end of the read [%d]\n", opt->indel_ends_bound);
  fprintf(stderr, "         -b INT      stop searching when INT optimal CALs have been found [%d]\n", opt->max_best_cals);
  // adapter trimming ?
  fprintf(stderr, "         -q INT      the queue size for the reads [%d]\n", opt->reads_queue_size);
  fprintf(stderr, "         -n INT      the number of threads [%d]\n", opt->num_threads);
  fprintf(stderr, "         -h          print this message\n");
  fprintf(stderr, "\n");

  free(reads_format);

  return 1;
}

int 
fmap_map1(int argc, char *argv[])
{
  int c;
  fmap_map1_opt_t opt;

  // program defaults
  opt.fn_fasta = opt.fn_reads = NULL;
  opt.max_mm = -1; opt.max_mm_frac = 0.02; // TODO: move this to a define block 
  opt.max_gapo = -1; opt.max_gapo_frac = 0.01; // TODO: move this to a define block
  opt.max_gape = -1; opt.max_gape_frac = 0.10; // TODO: move this to a define block
  opt.pen_mm = 3; opt.pen_gapo = 11; opt.pen_gape = 4; // TODO: move this to a define block
  opt.max_cals_del = 10; // TODO: move this to a define block
  opt.indel_ends_bound = 5; // TODO: move this to a define block
  opt.max_best_cals = 32; // TODO: move this to a define block
  opt.reads_format = FMAP_READS_FORMAT_FASTQ;
  opt.reads_queue_size = 65536; // TODO: move this to a define block
  opt.num_threads = 1;

  while((c = getopt(argc, argv, "f:r:F:q:n:h")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = fmap_strdup(optarg); break;
        case 'r':
          opt.fn_reads = fmap_strdup(optarg); break;
        case 'F':
          opt.reads_format = fmap_get_reads_file_format_int(optarg); break;
        case 'q': 
          opt.reads_queue_size = atoi(optarg); break;
        case 'n':
          opt.num_threads = atoi(optarg); break;
        case 'h':
        default:
          return usage(&opt);
      }
  }

  if(argc != optind || 1 == argc) {
      return usage(&opt);
  }
  if(NULL == opt.fn_fasta) {
      fmap_error("required option -f", Exit, CommandLineArgument);
  }
  if(NULL == opt.fn_reads) {
      fmap_error("required option -r", Exit, CommandLineArgument);
  }
  if(FMAP_READS_FORMAT_UNKNOWN == opt.reads_format) {
      fmap_error("option -F unrecognized", Exit, CommandLineArgument);
  }
  if(FMAP_READS_FORMAT_SFF == opt.reads_format) {
      fmap_error("SFF currently not supported", Exit, CommandLineArgument);
  }
  if(opt.reads_queue_size <= 0) {
      fmap_error("option -q out of range", Exit, CommandLineArgument);
  }
  if(opt.num_threads <= 0) {
      fmap_error("option -n out of range", Exit, CommandLineArgument);
  }

  fmap_map1_core(&opt);

  free(opt.fn_fasta);
  free(opt.fn_reads);

  return 0;
}
