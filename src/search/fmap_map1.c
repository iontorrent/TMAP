#include <stdlib.h>
#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "../util/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt_gen.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_bwt_match.h"
#include "../index/fmap_sa.h"
#include "../io/fmap_seq_io.h"
#include "fmap_map1_aux.h"
#include "fmap_map1.h"

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
static pthread_mutex_t fmap_map1_read_lock = PTHREAD_MUTEX_INITIALIZER;
static int32_t fmap_map1_read_lock_tid = 0;
#define FMAP_MAP1_THREAD_BLOCK_SIZE 1024
#endif

static inline int 
fmap_map1_print_sam(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_bwt_t *bwt, fmap_sa_t *sa, fmap_map1_aln_t *a)
{
  uint32_t i, j, n = 0;

  for(i=a->k;i<=a->l;i++) {
      uint16_t flag = 0x0000;
      uint32_t pos = 0, seqid = 0, pacpos = 0; 
      int32_t aln_ref_l = 0;

      // get the number of non-inserted bases 
      for(j=0;j<a->cigar_length;j++) {
          switch((a->cigar[j] & 0xf)) {
            case BAM_CMATCH:
            case BAM_CDEL:
              aln_ref_l += (a->cigar[j] >> 4); break;
            default:
              break;
          }
      }

      // SA position to packed refseq position
      pacpos = bwt->seq_len - fmap_sa_pac_pos(sa, bwt, i) - aln_ref_l + 1;

      if(0 <= fmap_refseq_pac2real(refseq, pacpos, seq->seq.l, &seqid, &pos)) {
          if(1 == a->strand) { // reverse for the output
              flag |= 0x0010;
              fmap_seq_reverse(seq, 1);
          }
          fprintf(stdout, "%s\t%u\t%s\t%u\t%u\t",
                  seq->name.s, flag, refseq->annos[seqid].name,
                  pos, 255);
          for(j=0;j<a->cigar_length;j++) {
              fprintf(stdout, "%d%c", (a->cigar[j]>>4), "MIDNSHP"[a->cigar[j] & 0xf]);
          }
          fprintf(stdout, "\t*\t0\t0\t%s\t%s\n",
                  seq->seq.s, seq->qual.s);
          if(1 == a->strand) { // reverse back
              fmap_seq_reverse(seq, 1);
          }
          n++;
      }
  }

  return n;
}

static inline void
fmap_map1_print_sam_unmapped(fmap_seq_t *seq)
{
  uint16_t flag = 0x0004;
  fprintf(stdout, "%s\t%u\t%s\t%u\t%u\t*\t*\t0\t0\t%s\t%s\n",
          seq->name.s, flag, "*",
          0, 0, seq->seq.s, seq->qual.s);
}


static int 
fmap_map1_cal_width(const fmap_bwt_t *rbwt, int len, const char *str, fmap_map1_width_t *width)
{
  // TODO: update based on new BWT ordering
  uint32_t k, l, ok, ol;
  int i, bid;
  bid = 0;
  k = 0; l = rbwt->seq_len;
  for (i = 0; i < len; ++i) {
      uint8_t c = (int)str[i];
      if (c < 4) {
          fmap_bwt_2occ(rbwt, k - 1, l, c, &ok, &ol);
          k = rbwt->L2[c] + ok + 1;
          l = rbwt->L2[c] + ol;
      }
      if (k > l || c > 3) { // then restart
          k = 0;
          l = rbwt->seq_len;
          ++bid;
      }
      width[i].w = l - k + 1;
      width[i].bid = bid;
  }
  width[len].w = 0;
  width[len].bid = ++bid;
  return bid;
}

static void
fmap_map1_core_worker(fmap_seq_t **seq_buffer, int32_t seq_buffer_length, fmap_map1_aln_t ***alns,
                      fmap_bwt_t *bwt, int32_t tid, fmap_map1_opt_t *opt)
{
  int32_t low = 0, high;
  fmap_map1_width_t *width[2]={NULL,NULL}, *seed_width[2]={NULL,NULL};
  int32_t width_length = 0;
  fmap_map1_aux_stack_t *stack;

  seed_width[0] = fmap_calloc(opt->seed_length+1, sizeof(fmap_map1_width_t), "seed_width[0]");
  seed_width[1] = fmap_calloc(opt->seed_length+1, sizeof(fmap_map1_width_t), "seed_width[1]");

  stack = fmap_map1_aux_stack_init();

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
      while(low<high) {
          fmap_map1_opt_t opt_local = (*opt); // copy over values
          fmap_seq_t *seq[2]={NULL, NULL}, *orig_seq=NULL;

          orig_seq = seq_buffer[low];

          // clone the sequence and get the reverse compliment
          seq[0] = fmap_seq_clone(orig_seq);
          seq[1] = fmap_seq_clone(orig_seq);
          fmap_seq_reverse(seq[1], 1);
          // convert to integers
          fmap_seq_to_int(seq[0]);
          fmap_seq_to_int(seq[1]);

          opt_local.max_mm = (opt->max_mm < 0) ? (int)(opt->max_mm_frac * orig_seq->seq.l) : opt->max_mm; 
          opt_local.max_gape = (opt->max_gape < 0) ? (int)(opt->max_gape_frac * orig_seq->seq.l) : opt->max_gape; 
          opt_local.max_gapo = (opt->max_gapo < 0) ? (int)(opt->max_gapo_frac * orig_seq->seq.l) : opt->max_gapo; 

          if(width_length < orig_seq->seq.l) {
              free(width[0]); free(width[1]);
              width_length = orig_seq->seq.l;
              width[0] = fmap_calloc(width_length+1, sizeof(fmap_map1_width_t), "width[0]");
              width[1] = fmap_calloc(width_length+1, sizeof(fmap_map1_width_t), "width[1]");
          }
          fmap_map1_cal_width(bwt, seq[0]->seq.l, seq[0]->seq.s, width[0]);
          fmap_map1_cal_width(bwt, seq[1]->seq.l, seq[1]->seq.s, width[1]);

          if(orig_seq->seq.l < opt->seed_length) {
              opt_local.seed_length = -1;
          }
          else {
              fmap_map1_cal_width(bwt, opt->seed_length, seq[0]->seq.s, width[0]);
              fmap_map1_cal_width(bwt, opt->seed_length, seq[1]->seq.s, width[1]);
          }

          alns[low] = fmap_map1_aux_core(seq, bwt, width, (0 < opt_local.seed_length) ? seed_width : NULL, &opt_local, stack);

          low++;

          // destroy
          fmap_seq_destroy(seq[0]);
          fmap_seq_destroy(seq[1]);
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
  uint32_t i, j;
  int32_t seq_buffer_length;
  fmap_refseq_t *refseq=NULL;
  fmap_bwt_t *bwt=NULL;
  fmap_sa_t *sa=NULL;
  fmap_file_t *fp_reads=NULL;
  fmap_seq_io_t *seqio = NULL;
  fmap_seq_t **seq_buffer = NULL;
  fmap_map1_aln_t ***alns = NULL;

  // For suffix search we need the reverse bwt/sa and forward refseq
  refseq = fmap_refseq_read(opt->fn_fasta, 0);
  bwt = fmap_bwt_read(opt->fn_fasta, 1);
  sa = fmap_sa_read(opt->fn_fasta, 1);

  // SAM header
  for(i=0;i<refseq->num_annos;i++) {
      fprintf(stdout, "@SQ\tSN:%s\tLN:%d\n",
              refseq->annos[i].name, (int)refseq->annos[i].len);
  }

  // TODO: support FASTA/FASTQ/SFF
  fp_reads = fmap_file_fopen(opt->fn_reads, "rb", FMAP_FILE_NO_COMPRESSION);
  seqio = fmap_seq_io_init(fp_reads);

  // initialize the buffer
  seq_buffer = fmap_malloc(sizeof(fmap_seq_t*)*opt->reads_queue_size, "seq_buffer");
  alns = fmap_malloc(sizeof(fmap_map1_aln_t**)*opt->reads_queue_size, "alns");
  for(i=0;i<opt->reads_queue_size;i++) {
      seq_buffer[i] = fmap_seq_init();
  }

  while(0 < (seq_buffer_length = fmap_seq_io_read_buffer(seqio, seq_buffer, opt->reads_queue_size))) {
      // do alignment
#ifdef HAVE_LIBPTHREAD
      if(1 == opt->num_threads) {
          fmap_map1_core_worker(seq_buffer, seq_buffer_length, alns, bwt, 0, opt);
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
              thread_data[i].alns = alns;
              thread_data[i].bwt = bwt;
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
      fmap_map1_core_worker(seq_buffer, seq_buffer_length, alns, bwt, 0, opt);
#endif


      for(i=0;i<seq_buffer_length;i++) {
          int32_t n_alns = 0, n_mapped = 0;
          fmap_map1_aln_t **a = alns[i];

          // get the number of alignments
          if(a != NULL) {
              while(NULL != a[n_alns]) {
                  n_alns++;
              }
          }

          // print alignments
          if(0 < n_alns) {
              for(j=0;j<n_alns;j++) { 
                  n_mapped += fmap_map1_print_sam(seq_buffer[i], refseq, bwt, sa, alns[i][j]);
              }
          }
          if(0 == n_mapped) {
              fmap_map1_print_sam_unmapped(seq_buffer[i]);
          }

          // free alignments
          for(j=0;j<n_alns;j++) { 
              free(alns[i][j]->cigar); // free cigar
              free(alns[i][j]); // free alignment
          }
          free(alns[i]);
          alns[i] = NULL;
      }
  }

  // free memory
  for(i=0;i<opt->reads_queue_size;i++) {
      fmap_seq_destroy(seq_buffer[i]);
  }
  free(seq_buffer);
  free(alns);
  fmap_file_fclose(fp_reads);
  fmap_refseq_destroy(refseq);
  fmap_bwt_destroy(bwt);
  fmap_sa_destroy(sa);
  fmap_seq_io_destroy(seqio);
}

static int 
usage(fmap_map1_opt_t *opt)
{
  char *reads_format = fmap_get_reads_file_format_string(opt->reads_format);
  // Future options:
  // - adapter trimming ?
  // - homopolymer enumeration ?
  // - add option to try various seed offsets

  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: %s map1 [options]", PACKAGE);
  fprintf(stderr, "\n");
  fprintf(stderr, "Options (required):\n");
  fprintf(stderr, "         -f FILE     the FASTA reference file name [%s]\n", opt->fn_fasta);
  fprintf(stderr, "         -r FILE     the reads file name [%s]\n", opt->fn_reads);
  fprintf(stderr, "Options (optional):\n");
  fprintf(stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff) [%s]\n", reads_format);
  fprintf(stderr, "         -l INT      the k-mer length to seed CALs (-1 to disable) [%d]\n", opt->seed_length);
  fprintf(stderr, "         -k INT      maximum number of mismatches in the seed [%d]\n", opt->seed_max_mm);

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
  opt.reads_format = FMAP_READS_FORMAT_FASTQ;
  opt.seed_length = 32; // move this to a define block
  opt.seed_max_mm = 3; // move this to a define block
  opt.max_mm = -1; opt.max_mm_frac = 0.02; // TODO: move this to a define block 
  opt.max_gapo = -1; opt.max_gapo_frac = 0.01; // TODO: move this to a define block
  opt.max_gape = -1; opt.max_gape_frac = 0.10; // TODO: move this to a define block
  opt.pen_mm = 3; opt.pen_gapo = 11; opt.pen_gape = 4; // TODO: move this to a define block
  opt.max_cals_del = 10; // TODO: move this to a define block
  opt.indel_ends_bound = 5; // TODO: move this to a define block
  opt.max_best_cals = 32; // TODO: move this to a define block
  opt.reads_queue_size = 65536; // TODO: move this to a define block
  opt.num_threads = 1;

  while((c = getopt(argc, argv, "f:r:F:l:k:m:o:e:M:O:E:d:i:b:q:n:h")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = fmap_strdup(optarg); break;
        case 'r':
          opt.fn_reads = fmap_strdup(optarg); break;
        case 'F':
          opt.reads_format = fmap_get_reads_file_format_int(optarg); break;
        case 'l':
          opt.seed_length = atoi(optarg); break;
        case 'k':
          opt.seed_max_mm = atoi(optarg); break;
        case 'm':
          if(NULL != strstr(optarg, ".")) opt.max_mm = -1, opt.max_mm_frac = atof(optarg);
          else opt.max_mm = atoi(optarg), opt.max_mm_frac = -1.0;
          break;
        case 'o':
          if(NULL != strstr(optarg, ".")) opt.max_gapo = -1, opt.max_gapo_frac = atof(optarg);
          else opt.max_gapo = atoi(optarg), opt.max_gapo_frac = -1.0;
          break;
        case 'e':
          if(NULL != strstr(optarg, ".")) opt.max_gape = -1, opt.max_gape_frac = atof(optarg);
          else opt.max_gape = atoi(optarg), opt.max_gape_frac = -1.0;
          break;
        case 'M':
          opt.pen_mm = atoi(optarg); break;
        case 'O':
          opt.pen_gapo = atoi(optarg); break;
        case 'E':
          opt.pen_gape = atoi(optarg); break;
        case 'd':
          opt.max_cals_del = atoi(optarg); break;
        case 'i':
          opt.indel_ends_bound = atoi(optarg); break;
        case 'b':
          opt.max_best_cals = atoi(optarg); break;
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
  else { // check command line arguments
      if(NULL == opt.fn_fasta) fmap_error("required option -f", Exit, CommandLineArgument);
      if(NULL == opt.fn_reads) fmap_error("required option -r", Exit, CommandLineArgument);
      if(FMAP_READS_FORMAT_UNKNOWN == opt.reads_format) fmap_error("option -F unrecognized", Exit, CommandLineArgument);
      if(FMAP_READS_FORMAT_SFF == opt.reads_format) fmap_error("SFF currently not supported", Exit, CommandLineArgument); // TODO: support SFF
      if(-1 != opt.seed_length) fmap_error_cmd_check_int(opt.seed_length, 1, INT32_MAX, "-s");

      // this will take care of the case where they are both < 0
      fmap_error_cmd_check_int((opt.max_mm_frac < 0) ? opt.max_mm : (int32_t)opt.max_mm_frac, 0, INT32_MAX, "-m"); 
      // this will take care of the case where they are both < 0
      fmap_error_cmd_check_int((opt.max_gapo_frac < 0) ? opt.max_gapo : (int32_t)opt.max_gapo_frac, 0, INT32_MAX, "-m"); 
      // this will take care of the case where they are both < 0
      fmap_error_cmd_check_int((opt.max_gape_frac < 0) ? opt.max_gape : (int32_t)opt.max_gape_frac, 0, INT32_MAX, "-m"); 

      fmap_error_cmd_check_int(opt.pen_mm, 0, INT32_MAX, "-M");
      fmap_error_cmd_check_int(opt.pen_gapo, 0, INT32_MAX, "-O");
      fmap_error_cmd_check_int(opt.pen_gape, 0, INT32_MAX, "-E");
      fmap_error_cmd_check_int(opt.max_cals_del, 1, INT32_MAX, "-d");
      fmap_error_cmd_check_int(opt.indel_ends_bound, 0, INT32_MAX, "-i");
      fmap_error_cmd_check_int(opt.max_best_cals, 0, INT32_MAX, "-b");
      fmap_error_cmd_check_int(opt.reads_queue_size, 1, INT32_MAX, "-q");
      fmap_error_cmd_check_int(opt.num_threads, 1, INT32_MAX, "-n");
  }

  fmap_map1_core(&opt);

  free(opt.fn_fasta);
  free(opt.fn_reads);

  return 0;
}
