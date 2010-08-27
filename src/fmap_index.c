#include <stdlib.h>
#include <string.h>
#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_refseq.h"
#include "fmap_bwt_gen.h"
#include "fmap_bwt.h"
#include "fmap_sa.h"
#include "fmap_progress.h"
#include "fmap_index.h"

// TODO
// - first-level hashing
// - BWT of the reverse sequence, if useful
// - Timing data structure ?
// - standard messaging ?

static void fmap_index_core(fmap_index_opt_t *opt)
{
  uint64_t ref_len = 0;

  // pack the reference sequence
  ref_len = fmap_refseq_fasta2pac(opt->fn_fasta, FMAP_FILE_NO_COMPRESSION);

  // check returned genome size
  if(opt->is_large < 0) {
      if(FMAP_INDEX_TOO_BIG_GENOME <= ref_len) { // too big (2^32 - 1)!
          fmap_error("Reference sequence too large", Exit, OutOfRange);
      }
      else if(FMAP_INDEX_LARGE_GENOME <= ref_len) { 
          opt->is_large = 1;
      }
      else {
          opt->is_large = 0;
      }
  }

  // create the bwt 
  fmap_bwt_pac2bwt(opt->fn_fasta, opt->is_large, opt->occ_interval);

  // create the suffix array
  fmap_sa_bwt2sa(opt->fn_fasta, opt->sa_interval);
}

static int usage(fmap_index_opt_t *opt)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: %s index [optionsn", PACKAGE);
  fprintf(stderr, "\n");
  fprintf(stderr, "Options (required):\n");
  fprintf(stderr, "         -f FILE     the FASTA file name to index\n");
  fprintf(stderr, "Options (optional):\n");
  fprintf(stderr, "         -o INT      the occurrence interval [%d]\n", opt->occ_interval);
  fprintf(stderr, "         -i INT      the suffix array interval [%d]\n", opt->sa_interval);
  fprintf(stderr, "         -a STRING   override BWT construction algorithm:\n");
  fprintf(stderr, "                     \t\"bwtsw\" (large genomes)\n");
  fprintf(stderr, "                     \t\"is\" (short genomes)\n");
  fprintf(stderr, "         -h          print this message\n");
  fprintf(stderr, "         -v          print verbose progress information\n");
  fprintf(stderr, "\n");
  return 1;
}

int fmap_index(int argc, char *argv[])
{
  int c;
  fmap_index_opt_t opt;

  opt.fn_fasta = NULL;
  opt.occ_interval = FMAP_BWT_OCC_INTERVAL; 
  opt.sa_interval = FMAP_SA_INTERVAL; 
  opt.is_large = -1;

  // Set output progress
  fmap_progress_set_command(argv[0]);

  while((c = getopt(argc, argv, "f:o:i:a:hv")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = fmap_strdup(optarg); break;
        case 'o':
          opt.occ_interval = atoi(optarg); break;
        case 'i':
          opt.sa_interval = atoi(optarg); break;
        case 'a':
          if(0 == strcmp("is", optarg)) opt.is_large = 0;
          else if(0 == strcmp("bwtsw", optarg)) opt.is_large = 1;
          else fmap_error("Option -a value not correct", Exit, CommandLineArgument); 
          break; 
        case 'v':
          fmap_progress_set_verbosity(1); break;
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
  if(0 < opt.occ_interval && 0 != (opt.occ_interval % 16)) {
      fmap_error("option -o out of range", Exit, CommandLineArgument);
  }
  if(0 < opt.sa_interval && 0 != (opt.sa_interval % 2)) {
      fmap_error("option -i out of range", Exit, CommandLineArgument);
  }

  fmap_index_core(&opt);

  free(opt.fn_fasta);

  return 0;
}
