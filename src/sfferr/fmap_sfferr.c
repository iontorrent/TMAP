#include <stdlib.h>
#include <config.h>
#include <math.h>
#include <stdint.h>
#ifdef HAVE_SAMTOOLS
#include <sam.h>
#endif
#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "../util/fmap_progress.h"
#include "../index/fmap_refseq.h"
#include "../server/fmap_shm.h"
#include "../seq/fmap_seq.h"
#include "../io/fmap_seq_io.h"
#include "fmap_sfferr.h"

#ifdef HAVE_SAMTOOLS

static void 
fmap_sfferr_core(fmap_sfferr_opt_t *opt)
{
  fmap_refseq_t *refseq = NULL;
  fmap_shm_t *shm = NULL;
  samfile_t *fp_sam= NULL;
  fmap_file_t *fp_sff = NULL;
  fmap_sff_io_t *sffio = NULL;
  bam1_t *bam = NULL;
  fmap_sff_t *sff = NULL;
  char *name_bam_prev = NULL;

  if(0 == opt->shm_key) {
      fmap_progress_print("reading in reference data");
      refseq = fmap_refseq_read(opt->fn_fasta, 0);
      fmap_progress_print2("reference data read in");
  }
  else {
      fmap_progress_print("retrieving reference data from shared memory");
      shm = fmap_shm_init(opt->shm_key, 0, 0);
      if(NULL == (refseq = fmap_refseq_shm_unpack(fmap_shm_get_buffer(shm, FMAP_SHM_LISTING_REFSEQ)))) {
          fmap_error("the packed reference sequence was not found in shared memory", Exit, SharedMemoryListing);
      }
      fmap_progress_print2("reference data retrieved from shared memory");
  }

  // open the input files
  fp_sff = fmap_file_fopen(opt->fn_sff, "rb", opt->input_compr);
  sffio = fmap_sff_io_init(fp_sff);
  fp_sam = samopen(opt->fn_sam, (0 == opt->is_bam) ? "r" : "rb", NULL);

  bam = bam_init1();
  sff = fmap_sff_init();
  while(0 < samread(fp_sam, bam) 
        && 0 < fmap_sff_io_read(sffio, sff)) {
      // check that the alignment does not have two of the same records
      if(NULL != name_bam_prev && 0 == strcmp(bam1_qname(bam), name_bam_prev)) {
          fmap_error("two SAM/BAM records with the same read name found", Exit, OutOfRange); 
      }
      free(name_bam_prev);
      name_bam_prev = fmap_strdup(bam1_qname(bam));

      // check that the read names between the SAM/BAM and SFF match
      if(0 != strcmp(bam1_qname(bam), sff->rheader->name->s)) {
          fmap_file_fprintf(fmap_file_stderr, "SAM/BAM: \"%s\"\nSFF: \"%s\"\n",
                            bam1_qname(bam),
                            sff->rheader->name->s);
          fmap_error("the read names from the SAM/BAM does not match the SFF", Exit, OutOfRange);
      }

      // process
      // TODO

      // free structures
      bam_destroy1(bam);
      bam = bam_init1();
  }
  // free structures
  free(name_bam_prev);
  bam_destroy1(bam);
  fmap_sff_destroy(sff);

  // close the input files
  fmap_sff_io_destroy(sffio);
  fmap_file_fclose(fp_sff);
  samclose(fp_sam);
}

static int 
usage(fmap_sfferr_opt_t *opt)
{
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Usage: %s sfferr [options]", PACKAGE);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (optional):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -f FILE     the FASTA reference file name [%s]\n", opt->fn_fasta);
  fmap_file_fprintf(fmap_file_stderr, "         -r FILE     the sff file name [%s]\n", opt->fn_sff);
  fmap_file_fprintf(fmap_file_stderr, "         -S FILE     the SAM/BAM file name [%s]\n", opt->fn_sam);
  fmap_file_fprintf(fmap_file_stderr, "         -b          the input is a BAM file [%s]\n", (0 == opt->is_bam) ? "false" : "true");
  fmap_file_fprintf(fmap_file_stderr, "         -R FLOAT    the fraction of reads to sample [%.2lf]\n", opt->rand_sample_num);
  fmap_file_fprintf(fmap_file_stderr, "         -j          the input is bz2 compressed (bzip2) [%s]\n",
                    (FMAP_FILE_BZ2_COMPRESSION == opt->input_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -z          the input is gz compressed (gzip) [%s]\n",
                    (FMAP_FILE_GZ_COMPRESSION == opt->input_compr) ? "true" : "false");
  fmap_file_fprintf(fmap_file_stderr, "         -s INT      use shared memory with the following key [%d]\n", opt->shm_key);
  fmap_file_fprintf(fmap_file_stderr, "         -v          print verbose progress information\n");
  fmap_file_fprintf(fmap_file_stderr, "         -h          print this message\n");
  fmap_file_fprintf(fmap_file_stderr, "\n");

  return 1;
}

static fmap_sfferr_opt_t *
fmap_sfferr_opt_init()
{
  fmap_sfferr_opt_t *opt = NULL;

  opt = fmap_calloc(1, sizeof(fmap_sfferr_opt_t), "opt");

  // program defaults
  opt->fn_fasta = opt->fn_sff = opt->fn_sam = NULL;
  opt->is_bam = 0;
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

  while((c = getopt(argc, argv, "f:r:S:b:R:jzs:vh")) >= 0) {
      switch(c) {
        case 'f':
          opt->fn_fasta = fmap_strdup(optarg); break;
        case 'r':
          opt->fn_sff = fmap_strdup(optarg); 
          break;
        case 'S':
          opt->fn_sam = fmap_strdup(optarg); break;
        case 'R':
          opt->rand_sample_num = atof(optarg); break;
        case 's':
          opt->shm_key = atoi(optarg); break;
        case 'j':
          opt->input_compr = FMAP_FILE_BZ2_COMPRESSION; 
          break;
        case 'z':
          opt->input_compr = FMAP_FILE_GZ_COMPRESSION; 
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
      if(NULL == opt->fn_sff) { // do we also check the SFF format
          fmap_error("option -r must be specified", Exit, CommandLineArgument);
      }
      if(NULL == opt->fn_sam) {
          fmap_error("option -S must be specified", Exit, CommandLineArgument);
      }

      fmap_error_cmd_check_int(opt->rand_sample_num, 0, 1, "-R");
  }

  fmap_sfferr_core(opt);

  fmap_sfferr_opt_destroy(opt);
  
  fmap_progress_print2("terminating successfully");

  return 0;
}

#endif /* HAVE_SAMTOOLS */
