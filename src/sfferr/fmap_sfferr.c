#include <stdlib.h>
#include <config.h>
#include <math.h>
#ifdef HAVE_SAMTOOLS
#include <sam.h>
#endif
#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "../util/fmap_progress.h"
#include "../seq/fmap_seq.h"
#include "../io/fmap_seq_io.h"
#include "fmap_sfferr.h"

#ifdef HAVE_SAMTOOLS

static void 
fmap_sfferr_core(fmap_sfferr_opt_t *opt)
{
  samfile_t *fp_am = NULL;
}

static int 
usage(fmap_sfferr_opt_t *opt)
{
  char *reads_format = fmap_get_reads_file_format_string(opt->reads_format);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Usage: %s sfferr [options]", PACKAGE);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (optional):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -f FILE     the FASTA reference file name [%s]\n", opt->fn_fasta);
  fmap_file_fprintf(fmap_file_stderr, "         -r FILE     the reads file name [%s]\n", opt->fn_sff);
  fmap_file_fprintf(fmap_file_stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff) [%s]\n", reads_format);
  fmap_file_fprintf(fmap_file_stderr, "         -S FILE     the SAM/BAM file name [%s]\n", opt->fn_sam);
  fmap_file_fprintf(fmap_file_stderr, "         -R FLOAT    the fraction of reads to sample [%.2lf]\n", opt->rand_sample_num);
  fmap_file_fprintf(fmap_file_stderr, "         -j          the input is bz2 compressed (bzip2) [%s]\n",
                    (FMAP_FILE_BZ2_COMPRESSION == opt->input_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -z          the input is gz compressed (gzip) [%s]\n",
                    (FMAP_FILE_GZ_COMPRESSION == opt->input_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -s INT      use shared memory with the following key [%d]\n", opt->shm_key);
  fmap_file_fprintf(fmap_file_stderr, "         -v          print verbose progress information\n");
  fmap_file_fprintf(fmap_file_stderr, "         -h          print this message\n");
  fmap_file_fprintf(fmap_file_stderr, "\n");

  free(reads_format);

  return 1;
}

static fmap_sfferr_opt_t *
fmap_sfferr_opt_init()
{
  fmap_sfferr_opt_t *opt = NULL;

  opt = fmap_calloc(1, sizeof(fmap_sfferr_opt_t), "opt");

  // program defaults
  opt->fn_fasta = opt->fn_sff = opt->fn_sam = NULL;
  opt->reads_format = FMAP_READS_FORMAT_UNKNOWN;
  opt->rand_sample_num = 1.0;
  opt->input_compr = FMAP_FILE_NO_COMPRESSION;
  opt->shm_key = 0;

  return opt;
}

static void
fmap_sfferr_opt_destroy(fmap_sfferr_opt_t *opt)
{

  free(opt->fn_fasta);
  free(opt->fn_sff);
  free(opt->fn_sam);
  free(opt);
}

int 
fmap_sfferr_main(int argc, char *argv[])
{
  int c;
  fmap_sfferr_opt_t *opt = NULL;

  srand48(0); // random seed

  opt = fmap_sfferr_opt_init();
  opt->argv = argv;
  opt->argc = argc;

  while((c = getopt(argc, argv, "f:r:F:S:R:S:jzs:vh")) >= 0) {
      switch(c) {
        case 'f':
          opt->fn_fasta = fmap_strdup(optarg); break;
        case 'r':
          opt->fn_sff = fmap_strdup(optarg); 
          fmap_get_reads_file_format_from_fn_int(opt->fn_sff, &opt->reads_format, &opt->input_compr);
          break;
        case 'F':
          opt->reads_format = fmap_get_reads_file_format_int(optarg); break;
        case 'S':
          opt->fn_sam = fmap_strdup(optarg); break;
        case 'R':
          opt->rand_sample_num = atof(optarg); break;
        case 's':
          opt->shm_key = atoi(optarg); break;
        case 'j':
          opt->input_compr = FMAP_FILE_BZ2_COMPRESSION; 
          fmap_get_reads_file_format_from_fn_int(opt->fn_sff, &opt->reads_format, &opt->input_compr);
          break;
        case 'z':
          opt->input_compr = FMAP_FILE_GZ_COMPRESSION; 
          fmap_get_reads_file_format_from_fn_int(opt->fn_sff, &opt->reads_format, &opt->input_compr);
          break;
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
          fmap_error("option -f and option -s may not be specified together", Exit, CommandLineArgument)
          ;
      }
      if(NULL == opt->fn_sff) {
          fmap_error("option -r must be specified", Exit, CommandLineArgument);
      }
      if(NULL == opt->fn_sam) {
          fmap_error("option -S must be specified", Exit, CommandLineArgument);
      }
      if(FMAP_READS_FORMAT_UNKNOWN == opt->reads_format) {
          fmap_error("the reads format (-r) was unrecognized", Exit, CommandLineArgument);
      }

      fmap_error_cmd_check_int(opt->rand_sample_num, 0, 1, "-R");
  }

  fmap_sfferr_core(opt);

  fmap_progress_print2("terminating successfully");

  return 0;
}

#endif /* HAVE_SAMTOOLS */
