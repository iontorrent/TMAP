#include <stdlib.h>
#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_refseq.h"
#include "fmap_bwt_gen.h"
#include "fmap_bwt.h"
#include "fmap_sa.h"
#include "fmap_index.h"

// TODO
// - first-level hashing
// - BWT of the reverse sequence, if useful
// - Timing data structure ?
// - standard messaging ?

static void fmap_index_core(fmap_index_opt_t *opt)
{
  uint64_t ref_len = 0;
  int is_large = 0;

  // pack the reference sequence
  ref_len = fmap_refseq_fasta2pac(opt->fn_fasta, FMAP_FILE_NO_COMPRESSION);

  // check returned genome size
  if(FMAP_INDEX_TOO_BIG_GENOME <= ref_len) { // too big (2^32 - 1)!
      fmap_error("Reference sequence too large", Exit, OutOfRange);
  }
  else if(FMAP_INDEX_LARGE_GENOME <= ref_len) { 
      is_large = 1;
  }

  // create the bwt 
  fmap_bwt_pac2bwt(opt->fn_fasta, is_large, opt->occ_interval);

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
  fprintf(stderr, "         -h          print this message\n");
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

  while((c = getopt(argc, argv, "f:o:i:h")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = fmap_strdup(optarg); break;
        case 'o':
          opt.occ_interval = atoi(optarg); break;
        case 'i':
          opt.sa_interval = atoi(optarg); break;
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
