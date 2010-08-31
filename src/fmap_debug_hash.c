#include <stdlib.h>
#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_refseq.h"
#include "fmap_bwt_gen.h"
#include "fmap_bwt.h"
#include "fmap_sa.h"
#include "fmap_seq.h"
#include "fmap_definitions.h"
#include "fmap_debug_hash.h"

static void 
fmap_debug_hash_core(fmap_debug_hash_opt_t *opt)
{
  fmap_bwt_t *bwt=NULL;
  bwt = fmap_bwt_read(opt->fn_fasta, 0);
  fmap_bwt_gen_hash(bwt, opt->hash_width);
  fmap_bwt_destroy(bwt);
}

static int 
usage(fmap_debug_hash_opt_t *opt)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: %s hash [optionsn", PACKAGE);
  fprintf(stderr, "\n");
  fprintf(stderr, "Optiosn (required):\n");
  fprintf(stderr, "         -f FILE     the FASTA reference file name\n");
  fprintf(stderr, "         -w INT      the bwt hash width\n");
  fprintf(stderr, "Options (optional):\n");
  fprintf(stderr, "         -h          print this message\n");
  fprintf(stderr, "\n");
  return 1;
}

int 
fmap_debug_hash(int argc, char *argv[])
{
  int c;
  fmap_debug_hash_opt_t opt;

  opt.fn_fasta = NULL;
  opt.hash_width = 0;

  while((c = getopt(argc, argv, "f:w:h")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = fmap_strdup(optarg); break;
        case 'w':
          opt.hash_width = atoi(optarg); break;
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
  if(opt.hash_width <= 0) {
      fmap_error("required option -w", Exit, CommandLineArgument);
  }

  fmap_debug_hash_core(&opt);

  free(opt.fn_fasta);

  return 0;
}
