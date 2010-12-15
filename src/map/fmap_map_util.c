#include <stdlib.h>
#include <unistd.h>
#include "../util/fmap_alloc.h"
#include "../util/fmap_error.h"
#include "../util/fmap_sam.h"
#include "../util/fmap_progress.h"
#include "../util/fmap_definitions.h"
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_sa.h"
#include "../sw/fmap_sw.h"
#include "../sw/fmap_fsw.h"
#include "fmap_map_util.h"

fmap_map_opt_t *
fmap_map_opt_init(int32_t algo_id)
{
  int32_t i;
  fmap_map_opt_t *opt = NULL;

  opt = fmap_calloc(1, sizeof(fmap_map_opt_t), "opt");

  // internal data
  opt->algo_id = algo_id;
  opt->argv = NULL;
  opt->argc = -1;

  // global options 
  opt->fn_fasta = opt->fn_reads = NULL;
  opt->reads_format = FMAP_READS_FORMAT_UNKNOWN;
  opt->score_match = FMAP_MAP_UTIL_SCORE_MATCH;
  opt->pen_mm = FMAP_MAP_UTIL_PEN_MM;
  opt->pen_gapo = FMAP_MAP_UTIL_PEN_GAPO;
  opt->pen_gape = FMAP_MAP_UTIL_PEN_GAPE;
  opt->fscore = FMAP_MAP_UTIL_FSCORE;
  opt->flow = NULL;
  opt->bw = 50; 
  opt->aln_global = 0;
  opt->reads_queue_size = 65536; // TODO: move this to a define block
  opt->num_threads = 1;
  opt->aln_output_mode = FMAP_MAP_UTIL_ALN_MODE_RAND_BEST; // TODO: move this to a define block
  opt->sam_rg = NULL;
  opt->sam_sff_tags = 0;
  opt->input_compr = FMAP_FILE_NO_COMPRESSION;
  opt->output_compr = FMAP_FILE_NO_COMPRESSION;
  opt->shm_key = 0;

  switch(algo_id) {
    case FMAP_MAP_ALGO_MAP1:
      // map1
      opt->seed_length = 32; // move this to a define block
      opt->seed_max_mm = 3; // move this to a define block
      opt->max_mm = 2; opt->max_mm_frac = -1.;
      opt->max_gapo = 1; opt->max_gapo_frac = -1.;
      opt->max_gape = 6; opt->max_gape_frac = -1.;
      opt->max_cals_del = 10; // TODO: move this to a define block
      opt->indel_ends_bound = 5; // TODO: move this to a define block
      opt->max_best_cals = 32; // TODO: move this to a define block
      opt->max_entries = 2000000; // TODO: move this to a define block
      break;
    case FMAP_MAP_ALGO_MAP2:
      // map2
      opt->yita = 5.5f;
      //opt->mask_level = 0.50; 
      opt->length_coef = 5.5f;
      opt->score_thr = 30;
      opt->max_seed_intv = 3; 
      opt->z_best = 5; 
      opt->seeds_rev = 5;
      break;
    case FMAP_MAP_ALGO_MAP3:
      // map3
      opt->seed_length = -1; // move this to a define block
      opt->seed_length_set = 0;
      opt->max_seed_hits = 8; // move this to a define block
      opt->max_seed_band = 50; // move this to a define block
      opt->score_thr = 30;
      opt->hp_diff = 0;
      break;
    case FMAP_MAP_ALGO_MAPALL:
      // mapall
      opt->dup_window = 128;
      opt->aln_output_mode_ind = 0;
      for(i=0;i<2;i++) {
          opt->algos[i] = 0;
          opt->opt_map1[i] = fmap_map_opt_init(FMAP_MAP_ALGO_MAP1);
          opt->opt_map2[i] = fmap_map_opt_init(FMAP_MAP_ALGO_MAP2);
          opt->opt_map3[i] = fmap_map_opt_init(FMAP_MAP_ALGO_MAP3);
      }
      break;
    default:
      break;
  }

  return opt;
}

void
fmap_map_opt_destroy(fmap_map_opt_t *opt)
{
  int32_t i;

  free(opt->fn_fasta);
  free(opt->fn_reads); 
  free(opt->sam_rg);
  free(opt->flow);

  switch(opt->algo_id) {
    case FMAP_MAP_ALGO_MAP1:
    case FMAP_MAP_ALGO_MAP2:
    case FMAP_MAP_ALGO_MAP3:
      break;
    case FMAP_MAP_ALGO_MAPALL:
      // mapall
      opt->dup_window = 128;
      opt->aln_output_mode_ind = 0;
      for(i=0;i<2;i++) {
          opt->algos[i] = 0;
          fmap_map_opt_destroy(opt->opt_map1[i]);
          fmap_map_opt_destroy(opt->opt_map2[i]);
          fmap_map_opt_destroy(opt->opt_map3[i]);
      }
      break;
    default:
      break;
  }

  free(opt);
}

int
fmap_map_opt_usage(fmap_map_opt_t *opt)
{
  char *reads_format = fmap_get_reads_file_format_string(opt->reads_format);

  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Usage: %s %s [options]\n", fmap_algo_id_to_name(opt->algo_id), PACKAGE);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "global options (required):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -f FILE     the FASTA reference file name [%s]\n", opt->fn_fasta);
  fmap_file_fprintf(fmap_file_stderr, "         -r FILE     the reads file name [%s]\n", (NULL == opt->fn_reads) ? "stdin" : opt->fn_reads);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "global options (optional):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff) [%s]\n", reads_format);
  fmap_file_fprintf(fmap_file_stderr, "         -A INT      score for a match [%d]\n", opt->score_match);
  fmap_file_fprintf(fmap_file_stderr, "         -M INT      the mismatch penalty [%d]\n", opt->pen_mm);
  fmap_file_fprintf(fmap_file_stderr, "         -O INT      the indel start penalty [%d]\n", opt->pen_gapo);
  fmap_file_fprintf(fmap_file_stderr, "         -E INT      the indel extend penalty [%d]\n", opt->pen_gape);
  fmap_file_fprintf(fmap_file_stderr, "         -X INT      the flow score penalty [%d]\n", opt->fscore);
  /*
  fmap_file_fprintf(fmap_file_stderr, "         -x STRING   the flow order ([ACGT]{4}) [%s]\n",
                    (NULL == opt->flow) ? "not using" : opt->flow);
                    */
  fmap_file_fprintf(fmap_file_stderr, "         -w INT      the band width [%d]\n", opt->bw);
  fmap_file_fprintf(fmap_file_stderr, "         -T INT      score threshold divided by the match score [%d]\n", opt->score_thr);
  fmap_file_fprintf(fmap_file_stderr, "         -q INT      the queue size for the reads (-1 disables) [%d]\n", opt->reads_queue_size);
  fmap_file_fprintf(fmap_file_stderr, "         -n INT      the number of threads [%d]\n", opt->num_threads);
  fmap_file_fprintf(fmap_file_stderr, "         -a INT      output filter [%d]\n", opt->aln_output_mode);
  fmap_file_fprintf(fmap_file_stderr, "                             0 - unique best hits\n");
  fmap_file_fprintf(fmap_file_stderr, "                             1 - random best hit\n");
  fmap_file_fprintf(fmap_file_stderr, "                             2 - all best hits\n");
  fmap_file_fprintf(fmap_file_stderr, "                             3 - all alignments\n");
  fmap_file_fprintf(fmap_file_stderr, "         -R STRING   the RG line in the SAM header [%s]\n", opt->sam_rg);
  fmap_file_fprintf(fmap_file_stderr, "         -Y          include SFF specific SAM tags [%s]\n",
                    (1 == opt->sam_sff_tags) ? "true" : "false");
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

  fmap_file_fprintf(fmap_file_stderr, "%s options (optional):\n", fmap_algo_id_to_name(opt->algo_id));
  switch(opt->algo_id) {
    case FMAP_MAP_ALGO_MAP1:
      // map1
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
      fmap_file_fprintf(fmap_file_stderr, "         -d INT      the maximum number of CALs to extend a deletion [%d]\n", opt->max_cals_del);
      fmap_file_fprintf(fmap_file_stderr, "         -i INT      indels are not allowed within INT number of bps from the end of the read [%d]\n", opt->indel_ends_bound);
      fmap_file_fprintf(fmap_file_stderr, "         -b INT      stop searching when INT optimal CALs have been found [%d]\n", opt->max_best_cals);
      fmap_file_fprintf(fmap_file_stderr, "         -Q INT      maximum number of alignment nodes [%d]\n", opt->max_entries);
      break;
    case FMAP_MAP_ALGO_MAP2:
      //fmap_file_fprintf(fmap_file_stderr, "         -y FLOAT    error recurrence coef. (4..16) [%.1lf]\n", opt->yita);
      //fmap_file_fprintf(fmap_file_stderr, "         -m FLOAT    mask level [%.2f]\n", opt->mask_level);
      fmap_file_fprintf(fmap_file_stderr, "         -c FLOAT    coefficient of length-threshold adjustment [%.1lf]\n", opt->length_coef);
      fmap_file_fprintf(fmap_file_stderr, "         -S INT      maximum seeding interval size [%d]\n", opt->max_seed_intv);
      fmap_file_fprintf(fmap_file_stderr, "         -b INT      Z-best [%d]\n", opt->z_best);
      fmap_file_fprintf(fmap_file_stderr, "         -N INT      # seeds to trigger reverse alignment [%d]\n", opt->seeds_rev);
      break;
    case FMAP_MAP_ALGO_MAP3:
      fmap_file_fprintf(fmap_file_stderr, "         -l INT      the k-mer length to seed CALs (-1 tunes to the genome size) [%d]\n", opt->seed_length);
      fmap_file_fprintf(fmap_file_stderr, "         -S INT      the maximum number of hits returned by a seed [%d]\n", opt->max_seed_hits);
      fmap_file_fprintf(fmap_file_stderr, "         -b INT      the window of bases in which to group seeds [%d]\n", opt->max_seed_band);
      fmap_file_fprintf(fmap_file_stderr, "         -H INT      single homopolymer error difference for enumeration [%d]\n", opt->hp_diff);
      break;
    case FMAP_MAP_ALGO_MAPALL:
      fmap_file_fprintf(fmap_file_stderr, "         -W INT      remove duplicate alignments from different algorithms within this bp window [%d]\n",
                        opt->dup_window);
      fmap_file_fprintf(fmap_file_stderr, "         -I          apply the output filter for each algorithm separately [%s]\n",
                        (1 == opt->aln_output_mode_ind) ? "true" : "false");
      break;
    default:
      break;
  }

  // free
  fmap_file_fprintf(fmap_file_stderr, "\n");
  free(reads_format);

  return 1;
}

int32_t
fmap_map_opt_parse(int argc, char *argv[], fmap_map_opt_t *opt)
{
  int c;
  int32_t found_option;
  char *getopt_format = NULL;

  opt->argc = argc; opt->argv = argv;
  switch(opt->algo_id) {
    case FMAP_MAP_ALGO_MAP1:
      getopt_format = fmap_strdup("f:r:F:A:M:O:E:X:x:gw:T:q:n:a:R:YjzJZs:vhl:k:m:o:e:d:i:b:Q:");
      break;
    case FMAP_MAP_ALGO_MAP2:
      getopt_format = fmap_strdup("f:r:F:A:M:O:E:X:x:gw:T:q:n:a:R:YjzJZs:vhc:S:b:N:");
      break;
    case FMAP_MAP_ALGO_MAP3:
      getopt_format = fmap_strdup("f:r:F:A:M:O:E:X:x:gw:T:q:n:a:R:YjzJZs:vhl:S:b:H:");
      break;
    case FMAP_MAP_ALGO_MAPALL:
      getopt_format = fmap_strdup("f:r:F:A:M:O:E:X:x:gw:T:q:n:a:R:YjzJZs:vhW:I");
      break;
    default:
      break;
  }

  // global options
  while((c = getopt(argc, argv, getopt_format)) >= 0) {
      found_option = 1;
      switch(c) { 
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
        case 'X':
          opt->fscore = atoi(optarg); break;
        case 'x':
          opt->flow = fmap_strdup(optarg); break;
        case 'g':
          opt->aln_global = 1; break;
        case 'w':
          opt->bw = atoi(optarg); break;
        case 'T':
          opt->score_thr = atoi(optarg); break;
        case 'q':
          opt->reads_queue_size = atoi(optarg); break;
        case 'n':
          opt->num_threads = atoi(optarg); break;
        case 'a':
          opt->aln_output_mode = atoi(optarg); break;
        case 'R':
          opt->sam_rg = fmap_strdup(optarg); break;
        case 'Y':
          opt->sam_sff_tags = 1; break;
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
          return 0;
          break;
        default:
          // algorithm-specific options
          switch(opt->algo_id) {
            case FMAP_MAP_ALGO_MAP1:
              switch(c) {
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
                case 'd':
                  opt->max_cals_del = atoi(optarg); break;
                case 'i':
                  opt->indel_ends_bound = atoi(optarg); break;
                case 'b':
                  opt->max_best_cals = atoi(optarg); break;
                case 'Q':
                  opt->max_entries = atoi(optarg); break;
                default:
                  return 0;
              }
              break;
            case FMAP_MAP_ALGO_MAP2:
              switch(c) {
                  /*
                     case 'y': 
                     opt->yita = atof(optarg); break;
                     */
                  /*
                     case 'm': 
                     opt->mask_level = atof(optarg); break;
                     */
                case 'c':
                  opt->length_coef = atof(optarg); break;
                case 'S':
                  opt->max_seed_intv = atoi(optarg); break;
                case 'b':
                  opt->z_best= atoi(optarg); break;
                case 'N':
                  opt->seeds_rev = atoi(optarg); break;
                  break;
                default: 
                  return 0;
              }
              break;
            case FMAP_MAP_ALGO_MAP3:
              switch(c) {
                case 'l':
                  opt->seed_length = atoi(optarg); opt->seed_length_set = 1; break;
                case 'S':
                  opt->max_seed_hits = atoi(optarg); break;
                case 'b':
                  opt->max_seed_band = atoi(optarg); break;
                case 'H':
                  opt->hp_diff = atoi(optarg); break;
                  break;
                default:
                  return 0;
              }
              break;
            case FMAP_MAP_ALGO_MAPALL:
              switch(c) {
                case 'W':
                  opt->dup_window = atoi(optarg); break;
                case 'I':
                  opt->aln_output_mode_ind = 1; break;
                default:
                  return 0;
              }
              break;
            default:
              break;
          }
          break;
      }
  }
  return 1;
}

static int32_t
fmap_map_opt_file_check_with_null(char *fn1, char *fn2)
{
    if(NULL == fn1 && NULL == fn2) {
        return 0;
    }
    else if((NULL == fn1 && NULL != fn2)
       || (NULL != fn1 && NULL == fn2)) {
        return 1;
    }
    else if(0 != strcmp(fn1, fn2)) {
        return 1;
    }
    return 0;
}

// for map1
#define __fmap_map_opt_check_common1(opt_map_all, opt_map_other) do { \
    if(0 != fmap_map_opt_file_check_with_null((opt_map_other)->fn_fasta, (opt_map_all)->fn_fasta)) { \
        fmap_error("option -f was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if(0 != fmap_map_opt_file_check_with_null((opt_map_other)->fn_reads, (opt_map_all)->fn_reads)) { \
        fmap_error("option -r was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->reads_format != (opt_map_all)->reads_format) { \
        fmap_error("option -F was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->pen_mm != (opt_map_all)->pen_mm) { \
        fmap_error("option -M was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->pen_gapo != (opt_map_all)->pen_gapo) { \
        fmap_error("option -O was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->pen_gape != (opt_map_all)->pen_gape) { \
        fmap_error("option -E was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if(0 != fmap_map_opt_file_check_with_null((opt_map_other)->flow, (opt_map_all)->flow)) { \
        fmap_error("option -x was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->reads_queue_size != (opt_map_all)->reads_queue_size) { \
        fmap_error("option -q was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->num_threads != (opt_map_all)->num_threads) { \
        fmap_error("option -n was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if(0 != fmap_map_opt_file_check_with_null((opt_map_other)->sam_rg, (opt_map_all)->sam_rg)) { \
        fmap_error("option -R was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->sam_sff_tags != (opt_map_all)->sam_sff_tags) { \
        fmap_error("option -Y was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->input_compr != (opt_map_all)->input_compr) { \
        fmap_error("option -j or -z was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->output_compr != (opt_map_all)->output_compr) { \
        fmap_error("option -J or -Z was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->shm_key != (opt_map_all)->shm_key) { \
        fmap_error("option -s was specified outside of the common options", Exit, CommandLineArgument); \
    } \
} while(0)

// for map2/map3
#define __fmap_map_opt_check_common2(opt_map_all, opt_map_other) do { \
    __fmap_map_opt_check_common1(opt_map_all, opt_map_other); \
    if((opt_map_other)->score_match != (opt_map_all)->score_match) { \
        fmap_error("option -A was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->fscore != (opt_map_all)->fscore) { \
        fmap_error("option -X was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->bw != (opt_map_all)->bw) { \
        fmap_error("option -w was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->aln_global != (opt_map_all)->aln_global) { \
        fmap_error("option -g was specified outside of the common options", Exit, CommandLineArgument); \
    } \
} while(0)


void
fmap_map_opt_check(fmap_map_opt_t *opt)
{
  int32_t i;
  // global options
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
  fmap_error_cmd_check_int(opt->score_match, 0, INT32_MAX, "-M");
  fmap_error_cmd_check_int(opt->pen_mm, 0, INT32_MAX, "-M");
  fmap_error_cmd_check_int(opt->pen_gapo, 0, INT32_MAX, "-O");
  fmap_error_cmd_check_int(opt->pen_gape, 0, INT32_MAX, "-E");
  fmap_error_cmd_check_int(opt->fscore, 0, INT32_MAX, "-X");
  if(NULL != opt->flow) fmap_error_cmd_check_int(strlen(opt->flow), 4, 4, "-x");
  fmap_error_cmd_check_int(opt->bw, 0, INT32_MAX, "-w");
  fmap_error_cmd_check_int(opt->score_thr, 0, INT32_MAX, "-T");
  if(-1 != opt->reads_queue_size) fmap_error_cmd_check_int(opt->reads_queue_size, 1, INT32_MAX, "-q");
  fmap_error_cmd_check_int(opt->num_threads, 1, INT32_MAX, "-n");
  fmap_error_cmd_check_int(opt->aln_output_mode, 0, 3, "-a");
  if(FMAP_FILE_BZ2_COMPRESSION == opt->output_compr
     && -1 == opt->reads_queue_size) {
      fmap_error("cannot buffer reads with bzip2 output (options \"-q 1 -J\")", Exit, OutOfRange);
  }

  switch(opt->algo_id) {
    case FMAP_MAP_ALGO_MAP1:
      // map1 options
      // this will take care of the case where they are both < 0
      fmap_error_cmd_check_int((opt->max_mm_frac < 0) ? opt->max_mm : (int32_t)opt->max_mm_frac, 0, INT32_MAX, "-m");
      // this will take care of the case where they are both < 0
      fmap_error_cmd_check_int((opt->max_gapo_frac < 0) ? opt->max_gapo : (int32_t)opt->max_gapo_frac, 0, INT32_MAX, "-m");
      // this will take care of the case where they are both < 0
      fmap_error_cmd_check_int((opt->max_gape_frac < 0) ? opt->max_gape : (int32_t)opt->max_gape_frac, 0, INT32_MAX, "-m");
      fmap_error_cmd_check_int(opt->max_cals_del, 1, INT32_MAX, "-d");
      fmap_error_cmd_check_int(opt->indel_ends_bound, 0, INT32_MAX, "-i");
      fmap_error_cmd_check_int(opt->max_best_cals, 0, INT32_MAX, "-b");
      fmap_error_cmd_check_int(opt->max_entries, 1, INT32_MAX, "-Q");
      if(-1 != opt->seed_length) fmap_error_cmd_check_int(opt->seed_length, 1, INT32_MAX, "-l");
      break;
    case FMAP_MAP_ALGO_MAP2:
      //fmap_error_cmd_check_int(opt->yita, 0, 1, "-y");
      //fmap_error_cmd_check_int(opt->mask_level, 0, 1, "-m");
      fmap_error_cmd_check_int(opt->length_coef, 0, INT32_MAX, "-c");
      fmap_error_cmd_check_int(opt->max_seed_intv, 0, INT32_MAX, "-S");
      fmap_error_cmd_check_int(opt->z_best, 1, INT32_MAX, "-Z");
      fmap_error_cmd_check_int(opt->seeds_rev, 0, INT32_MAX, "-N");
      break;
    case FMAP_MAP_ALGO_MAP3:
      if(-1 != opt->seed_length) fmap_error_cmd_check_int(opt->seed_length, 1, INT32_MAX, "-l");
      fmap_error_cmd_check_int(opt->max_seed_hits, 1, INT32_MAX, "-S");
      fmap_error_cmd_check_int(opt->max_seed_hits, 1, INT32_MAX, "-b");
      fmap_error_cmd_check_int(opt->hp_diff, 0, INT32_MAX, "-H");
      if(0 < opt->hp_diff && FMAP_SEQ_TYPE_SFF != opt->reads_format) fmap_error("-H option must be used with SFF only", Exit, OutOfRange); 
      break;
    case FMAP_MAP_ALGO_MAPALL:
      fmap_error_cmd_check_int(opt->dup_window, 0, INT32_MAX, "-W");
      fmap_error_cmd_check_int(opt->aln_output_mode_ind, 0, 1, "-I");
      if(0 == opt->algos[0]) {
          fmap_error("no algorithms given for stage 1", Exit, CommandLineArgument);
      }
      for(i=0;i<2;i++) {
          // check mapping algorithm specific options
          fmap_map_opt_check(opt->opt_map1[i]);
          fmap_map_opt_check(opt->opt_map2[i]);
          fmap_map_opt_check(opt->opt_map3[i]);

          // check that common values match other opt values
          __fmap_map_opt_check_common1(opt, opt->opt_map1[i]);
          __fmap_map_opt_check_common2(opt, opt->opt_map2[i]);
          __fmap_map_opt_check_common2(opt, opt->opt_map3[i]);
      }
      break;
    default:
      break;
  }
}

void
fmap_map_opt_print(fmap_map_opt_t *opt)
{
  fprintf(stderr, "algo_id=%d\n", opt->algo_id);
  fprintf(stderr, "fn_fasta=%s\n", opt->fn_fasta);
  fprintf(stderr, "fn_reads=%s\n", opt->fn_reads);
  fprintf(stderr, "reads_format=%d\n", opt->reads_format);
  fprintf(stderr, "score_match=%d\n", opt->score_match);
  fprintf(stderr, "pen_mm=%d\n", opt->pen_mm);
  fprintf(stderr, "pen_gapo=%d\n", opt->pen_gapo);
  fprintf(stderr, "pen_gape=%d\n", opt->pen_gape);
  fprintf(stderr, "fscore=%d\n", opt->fscore);
  fprintf(stderr, "flow=%s\n", opt->flow);
  fprintf(stderr, "aln_global=%d\n", opt->aln_global);
  fprintf(stderr, "reads_queue_size=%d\n", opt->reads_queue_size);
  fprintf(stderr, "num_threads=%d\n", opt->num_threads);
  fprintf(stderr, "aln_output_mode=%d\n", opt->aln_output_mode);
  fprintf(stderr, "sam_rg=%s\n", opt->sam_rg);
  fprintf(stderr, "sam_sff_tags=%d\n", opt->sam_sff_tags);
  fprintf(stderr, "input_compr=%d\n", opt->input_compr);
  fprintf(stderr, "output_compr=%d\n", opt->output_compr);
  fprintf(stderr, "shm_key=%d\n", (int)opt->shm_key);
  fprintf(stderr, "seed_length=%d\n", opt->seed_length);
  fprintf(stderr, "seed_length_set=%d\n", opt->seed_length_set);
  fprintf(stderr, "bw=%d\n", opt->bw);
  fprintf(stderr, "score_thr=%d\n", opt->score_thr);
  fprintf(stderr, "seed_max_mm=%d\n", opt->seed_max_mm);
  fprintf(stderr, "max_mm=%d\n", opt->max_mm);
  fprintf(stderr, "max_mm_frac=%lf\n", opt->max_mm_frac);
  fprintf(stderr, "max_gapo=%d\n", opt->max_gapo);
  fprintf(stderr, "max_gapo_frac=%lf\n", opt->max_gapo_frac);
  fprintf(stderr, "max_gape=%d\n", opt->max_gape);
  fprintf(stderr, "max_gape_frac=%lf\n", opt->max_gape_frac);
  fprintf(stderr, "max_cals_del=%d\n", opt->max_cals_del);
  fprintf(stderr, "indel_ends_bound=%d\n", opt->indel_ends_bound);
  fprintf(stderr, "max_best_cals=%d\n", opt->max_best_cals);
  fprintf(stderr, "max_entries=%d\n", opt->max_entries);
  fprintf(stderr, "yita=%lf\n", opt->yita);
  fprintf(stderr, "length_coef=%lf\n", opt->length_coef);
  fprintf(stderr, "max_seed_intv=%d\n", opt->max_seed_intv);
  fprintf(stderr, "z_best=%d\n", opt->z_best);
  fprintf(stderr, "seeds_rev=%d\n", opt->seeds_rev);
  fprintf(stderr, "max_seed_hits=%d\n", opt->max_seed_hits);
  fprintf(stderr, "max_seed_band=%d\n", opt->max_seed_band);
  fprintf(stderr, "hp_diff=%d\n", opt->hp_diff);
  fprintf(stderr, "dup_window=%d\n", opt->dup_window);
  fprintf(stderr, "aln_output_mode_ind=%d\n", opt->aln_output_mode_ind);
}

void
fmap_map_sam_malloc_aux(fmap_map_sam_t *s, int32_t algo_id)
{
  switch(s->algo_id) {
    case FMAP_MAP_ALGO_MAP1:
      s->aux.map1_aux = fmap_calloc(1, sizeof(fmap_map_map1_aux_t), "s->aux.map1_aux");
      break;
    case FMAP_MAP_ALGO_MAP2:
      s->aux.map2_aux = fmap_calloc(1, sizeof(fmap_map_map2_aux_t), "s->aux.map2_aux");
      break;
    case FMAP_MAP_ALGO_MAP3:
      s->aux.map3_aux = fmap_calloc(1, sizeof(fmap_map_map3_aux_t), "s->aux.map3_aux");
      break;
    default:
      break;
  }
}

inline void
fmap_map_sam_destroy_aux(fmap_map_sam_t *s)
{
  switch(s->algo_id) {
    case FMAP_MAP_ALGO_MAP1:
      free(s->aux.map1_aux);
      s->aux.map1_aux = NULL;
      break;
    case FMAP_MAP_ALGO_MAP2:
      free(s->aux.map2_aux);
      s->aux.map2_aux = NULL;
      break;
    case FMAP_MAP_ALGO_MAP3:
      free(s->aux.map3_aux);
      s->aux.map3_aux = NULL;
      break;
    default:
      break;
  }
}

void
fmap_map_sam_destroy(fmap_map_sam_t *s)
{
  fmap_map_sam_destroy_aux(s);
  free(s->cigar);
  s->cigar = NULL;
  s->n_cigar = 0;
}

fmap_map_sams_t *
fmap_map_sams_init()
{
  return fmap_calloc(1, sizeof(fmap_map_sams_t), "return");
}

void
fmap_map_sams_realloc(fmap_map_sams_t *s, int32_t n)
{
  int32_t i;
  if(n == s->n) return; 
  for(i=n;i<s->n;i++) {
      fmap_map_sam_destroy_aux(&s->sams[i]);
  }
  s->sams = fmap_realloc(s->sams, sizeof(fmap_map_sam_t) * n, "s->sams");
  s->n = n;
}

void
fmap_map_sams_destroy(fmap_map_sams_t *s)
{
  int32_t i;
  if(NULL == s) return;
  for(i=0;i<s->n;i++) {
      fmap_map_sam_destroy_aux(&s->sams[i]);
  }
  free(s->sams);
  free(s);
}

void
fmap_map_sam_copy_and_nullify(fmap_map_sam_t *dest, fmap_map_sam_t *src)
{
  (*dest) = (*src);
  src->cigar = NULL;
  switch(src->algo_id) {
    case FMAP_MAP_ALGO_MAP1:
      src->aux.map1_aux = NULL;
      break;
    case FMAP_MAP_ALGO_MAP2:
      src->aux.map2_aux = NULL;
      break;
    case FMAP_MAP_ALGO_MAP3:
      src->aux.map3_aux = NULL;
      break;
    default:
      break;
  }
}

static void
fmap_map_sam_print(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_map_sam_t *sam, int32_t sam_sff_tags)
{
  if(NULL == sam) { // unmapped
      fmap_sam_print_unmapped(fmap_file_stdout, seq, sam_sff_tags);
  }
  else {
      switch(sam->algo_id) {
        case FMAP_MAP_ALGO_MAP1:
          fmap_sam_print_mapped(fmap_file_stdout, seq, sam_sff_tags, refseq, 
                                sam->strand, sam->seqid, sam->pos,
                                sam->mapq, sam->cigar, sam->n_cigar,
                                sam->score, sam->algo_id, sam->algo_stage, 
                                "\tXM:i:%d\tXO:i:%d\tXG:i:%d",
                                sam->aux.map1_aux->n_mm,
                                sam->aux.map1_aux->n_gapo,
                                sam->aux.map1_aux->n_gape);
          break;
        case FMAP_MAP_ALGO_MAP2:
          if(0 < sam->aux.map2_aux->XI) {
              fmap_sam_print_mapped(fmap_file_stdout, seq, sam_sff_tags, refseq, 
                                    sam->strand, sam->seqid, sam->pos,
                                    sam->mapq, sam->cigar, sam->n_cigar,
                                    sam->score, sam->algo_id, sam->algo_stage, 
                                    "\tXS:i:%d\tXF:i:%d\tXE:i:%d\tXI:i:%d",
                                    sam->score_subo,
                                    sam->aux.map2_aux->XF, sam->aux.map2_aux->XE, 
                                    sam->aux.map2_aux->XI);
          }
          else {
              fmap_sam_print_mapped(fmap_file_stdout, seq, sam_sff_tags, refseq, 
                                    sam->strand, sam->seqid, sam->pos,
                                    sam->mapq, sam->cigar, sam->n_cigar,
                                    sam->score, sam->algo_id, sam->algo_stage, 
                                    "\tXS:i:%d\tXF:i:%d\tXE:i:%d",
                                    sam->score_subo,
                                    sam->aux.map2_aux->XF, sam->aux.map2_aux->XE);
          }
          break;
        case FMAP_MAP_ALGO_MAP3:
          fmap_sam_print_mapped(fmap_file_stdout, seq, sam_sff_tags, refseq, 
                                sam->strand, sam->seqid, sam->pos,
                                sam->mapq, sam->cigar, sam->n_cigar,
                                sam->score, sam->algo_id, sam->algo_stage, 
                                "\tXS:i:%d\tXE:i:%d",
                                sam->score_subo,
                                sam->aux.map3_aux->n_seeds);
          break;
      }
  }
}

void 
fmap_map_sams_print(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_map_sams_t *sams, int32_t sam_sff_tags) 
{
  int32_t i;
  if(0 < sams->n) {
      for(i=0;i<sams->n;i++) {
          fmap_map_sam_print(seq, refseq, &sams->sams[i], sam_sff_tags);
      }
  }
  else {
      fmap_map_sam_print(seq, refseq, NULL, sam_sff_tags);
  }
}

void
fmap_map_sams_filter1(fmap_map_sams_t *sams, int32_t aln_output_mode, int32_t algo_id)
{
  int32_t i, j, k;
  int32_t n_best = 0;
  int32_t best_score, cur_score;

  if(sams->n <= 1) {
      return;
  }

  best_score = INT32_MIN;
  n_best = 0;
  for(i=0;i<sams->n;i++) {
      if(FMAP_MAP_ALGO_NONE == algo_id
         || sams->sams[i].algo_id == algo_id) {
          cur_score = sams->sams[i].score;
          if(best_score < cur_score) {
              best_score = cur_score;
              n_best = 1;
          }
          else if(!(cur_score < best_score)) { // equal
              n_best++;
          }
      }
  }

  // adjust mapping qualities
  if(1 < n_best) {
      for(i=0;i<sams->n;i++) {
          if(FMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              sams->sams[i].mapq = 0;
          }
      }
  }
  else {
      for(i=0;i<sams->n;i++) {
          if(FMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              cur_score = sams->sams[i].score;
              if(cur_score < best_score) { // not the best
                  sams->sams[i].mapq = 0;
              }
          }
      }
  }

  if(FMAP_MAP_UTIL_ALN_MODE_ALL == aln_output_mode) {
      return;
  }

  // copy to the front
  if(n_best < sams->n) {
      for(i=j=0;i<sams->n;i++) {
          if(FMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              cur_score = sams->sams[i].score;
              if(cur_score < best_score) { // not the best
                  fmap_map_sam_destroy(&sams->sams[i]);
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                  }
                  j++;
              }
          }
          else {
              if(j < i) { // copy if we are not on the same index
                  fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
              }
              j++;
          }
      }
      // reallocate
      fmap_map_sams_realloc(sams, j);
  }

  if(FMAP_MAP_UTIL_ALN_MODE_UNIQ_BEST == aln_output_mode) {
      if(1 < n_best) { // there can only be one
          if(FMAP_MAP_ALGO_NONE == algo_id) {
              fmap_map_sams_realloc(sams, 0);
          }
          else {
              // get rid of all of them
              for(i=j=0;i<sams->n;i++) {
                  if(sams->sams[i].algo_id == algo_id) {
                      fmap_map_sam_destroy(&sams->sams[i]);
                  }
                  else {
                      if(j < i) { // copy if we are not on the same index
                          fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                      }
                      j++;
                  }
              }
              fmap_map_sams_realloc(sams, j);
          }
      }
  }
  else if(FMAP_MAP_UTIL_ALN_MODE_RAND_BEST == aln_output_mode) { // get a random
      int32_t r = (int32_t)(drand48() * n_best);

      // keep the rth one
      if(FMAP_MAP_ALGO_NONE == algo_id) {
          if(0 != r) {
              fmap_map_sam_copy_and_nullify(&sams->sams[0], &sams->sams[r]);
          }
          // reallocate
          fmap_map_sams_realloc(sams, 1);
      }
      else {
          // keep the rth one
          for(i=j=k=0;i<sams->n;i++) {
              if(sams->sams[i].algo_id == algo_id) {
                  if(k == r) { // keep
                      if(j < i) { // copy if we are not on the same index
                          fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                      }
                      j++;
                  }
                  else { // free
                      fmap_map_sam_destroy(&sams->sams[i]);
                  }
                  k++;
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                  }
                  j++;
              }
          }
          fmap_map_sams_realloc(sams, j);
      }
  }
  else if(FMAP_MAP_UTIL_ALN_MODE_ALL_BEST == aln_output_mode) {
      // do nothing
  }
  else {
      fmap_error("bug encountered", Exit, OutOfRange);
  }
}

void
fmap_map_sams_filter(fmap_map_sams_t *sams, int32_t aln_output_mode)
{
  fmap_map_sams_filter1(sams, aln_output_mode, FMAP_MAP_ALGO_NONE);
}

void
fmap_map_util_map1_adjust_score(fmap_map_sams_t *sams, int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape)
{
  int32_t i, j, n_match;
  for(i=0;i<sams->n;i++) {
      fmap_map_sam_t *sam = &sams->sams[i];
      n_match = 0 - sam->aux.map1_aux->n_mm;
      for(j=0;j<sam->n_cigar;j++) {
          switch(FMAP_SW_CIGAR_OP(sam->cigar[j])) {
            case BAM_CMATCH:
              n_match += FMAP_SW_CIGAR_LENGTH(sam->cigar[j]);
              break;
            default:
              break;
          }
      }

      // update the score
      sam->score = n_match * score_match;
      sam->score -= sam->aux.map1_aux->n_mm * pen_mm;
      sam->score -= sam->aux.map1_aux->n_gapo * pen_gapo;
      sam->score -= sam->aux.map1_aux->n_gape * pen_gape;
  }
}

void
fmap_map_util_fsw(fmap_sff_t *sff, 
                  fmap_map_sams_t *sams, fmap_refseq_t *refseq,
                  int32_t bw, int32_t aln_global, int32_t score_thr,
                  int32_t score_match, int32_t pen_mm, int32_t pen_gapo, 
                  int32_t pen_gape, int32_t fscore)
{
  int32_t i, j, k, l;
  uint8_t *target = NULL;
  int32_t target_mem = 0, target_len = 0;

  fmap_fsw_flowseq_t *fseq[2] = {NULL, NULL};
  fmap_fsw_path_t *path = NULL;
  int32_t path_mem = 0, path_len = 0;
  fmap_fsw_param_t param;
  int32_t matrix[25];

  // generate the alignment parameters
  param.matrix = matrix;
  param.band_width = 0;
  param.offset = FMAP_MAP_UTIL_FSW_OFFSET; // this sets the hp difference
  __fmap_fsw_gen_ap1(param, score_match, pen_mm, pen_gapo, pen_gape, fscore);

  // go through each hit
  for(i=0;i<sams->n;i++) {
      fmap_map_sam_t *s = &sams->sams[i];
      uint32_t ref_start, ref_end, pacpos;

      // get flow sequence if necessary
      if(NULL == fseq[s->strand]) {
          fseq[s->strand] = fmap_fsw_sff_to_flowseq(sff);
          if(1 == s->strand) fmap_fsw_flowseq_reverse_compliment(fseq[s->strand]);
      }

      param.band_width = 0;
      ref_start = ref_end = s->pos;
      for(j=0;j<s->n_cigar;j++) {
          int32_t op, op_len;

          op = FMAP_SW_CIGAR_OP(s->cigar[j]);
          op_len = FMAP_SW_CIGAR_LENGTH(s->cigar[j]);

          switch(op) {
            case BAM_CMATCH:
              ref_end += op_len;
              break;
            case BAM_CDEL:
              if(param.band_width < op_len) param.band_width += op_len;
              ref_end += op_len;
              break;
            case BAM_CINS:
              if(param.band_width < op_len) param.band_width += op_len;
              break;
            case BAM_CSOFT_CLIP:
              if(0 == j) {
                  if(ref_start <= op_len) {
                      ref_start = 1;
                  }
                  else {
                      ref_start = ref_start - op_len + 1;
                  }
              }
              else ref_end += op_len;
              break;
            default:
              // ignore
              break;
          }
      }

      // check bounds
      if(ref_start < 1) ref_start = 1;
      if(refseq->annos[s->seqid].len < ref_end) {
          ref_end = refseq->annos[s->seqid].len;
      }

      // get the target sequence
      target_len = ref_end - ref_start + 1;
      if(target_mem < target_len) {
          target_mem = target_len;
          fmap_roundup32(target_mem);
          target = fmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }
      for(pacpos=ref_start;pacpos<=ref_end;pacpos++) {
          // add contig offset and make zero-based
          target[pacpos-ref_start] =
            fmap_refseq_seq_i(refseq, pacpos + refseq->annos[s->seqid].offset-1);
      }

      // add to the band width
      param.band_width += 2 * bw;

      // make sure we have enough memory for the path
      while(path_mem <= target_len + fseq[s->strand]->num_flows) { // lengthen the path
          path_mem = target_len + fseq[s->strand]->num_flows;
          fmap_roundup32(path_mem);
          target = fmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }

      // re-align
      if(0 == aln_global) {
          s->score = fmap_fsw_local_core(target, target_len, fseq[s->strand], &param, path, &path_len, score_thr, &s->score_subo);
      }
      else {
          s->score = fmap_fsw_fitting_core(target, target_len, fseq[s->strand], &param, path, &path_len);
          s->score_subo = INT32_MIN;
      }


      if(path_len < 0) { // update
          s->pos = (ref_start-1) + (path[path_len-1].j-1);
          free(s->cigar);
          s->cigar = fmap_fsw_path2cigar(path, path_len, &s->n_cigar, 1);

          if(0 < path[path_len-1].i) { // skipped beginning flows
              // get the number of bases to clip
              for(j=k=0;j<path[path_len-1].i;j++) {
                  k += fseq[s->strand]->base_calls[j];
              }
              if(0 < k) { // bases should be soft-clipped
                  s->cigar = fmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
                  for(l=s->n_cigar-1;0<=l;l--) {
                      s->cigar[l+1] = s->cigar[l];
                  }
                  FMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, k);
                  s->n_cigar++;
              }
          }

          if(path[0].i < fseq[s->strand]->num_flows) { // skipped ending flows
              // get the number of bases to clip 
              for(j=path[0].i,k=0;j<fseq[s->strand]->num_flows;j++) {
                  k += fseq[s->strand]->base_calls[j];
              }
              if(0 < k) { // bases should be soft-clipped
                  s->cigar = fmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
                  s->cigar[s->n_cigar] = (k << 4) | 4;
                  s->n_cigar++;
              }
          }
      }
  }
  // free
  if(NULL != fseq[0]) fmap_fsw_flowseq_destroy(fseq[0]);
  if(NULL != fseq[1]) fmap_fsw_flowseq_destroy(fseq[1]);
  free(target);
  free(path);
}
