/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "tmap_refseq.h"
#include "tmap_bwt_gen.h"
#include "tmap_bwt.h"
#include "tmap_sa.h"
#include "tmap_index.h"

// TODO
// - first-level hashing
// - BWT of the reverse sequence, if useful
// - Timing data structure ?
// - standard messaging ?

static void tmap_index_core(tmap_index_opt_t *opt)
{
  uint64_t ref_len = 0;

  // pack the reference sequence
  ref_len = tmap_refseq_fasta2pac(opt->fn_fasta, TMAP_FILE_NO_COMPRESSION);
      
  if(TMAP_INDEX_TOO_BIG_GENOME <= ref_len) { // too big (2^32 - 1)!
      tmap_error("Reference sequence too large", Exit, OutOfRange);
  }

  // check returned genome size
  if(opt->is_large < 0) {
      if(TMAP_INDEX_LARGE_GENOME <= ref_len) { 
          opt->is_large = 1;
          tmap_progress_print("defaulting to \"bwtsw\" bwt construction algorithm");
      }
      else {
          opt->is_large = 0;
          tmap_progress_print("defaulting to \"is\" bwt construction algorithm");
      }
  }

  // create the bwt 
  tmap_bwt_pac2bwt(opt->fn_fasta, opt->is_large, opt->occ_interval, opt->hash_width);

  // create the suffix array
  tmap_sa_bwt2sa(opt->fn_fasta, opt->sa_interval);
}

static int usage(tmap_index_opt_t *opt)
{
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Usage: %s index [options]", PACKAGE);
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Options (required):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -f FILE     the FASTA file name to index\n");
  tmap_file_fprintf(tmap_file_stderr, "Options (optional):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -o INT      the occurrence interval [%d]\n", opt->occ_interval);
  tmap_file_fprintf(tmap_file_stderr, "         -w INT      the k-mer occurrence hash width [%d]\n", opt->hash_width);
  tmap_file_fprintf(tmap_file_stderr, "         -i INT      the suffix array interval [%d]\n", opt->sa_interval);
  tmap_file_fprintf(tmap_file_stderr, "         -a STRING   override BWT construction algorithm:\n");
  tmap_file_fprintf(tmap_file_stderr, "                     \t\"bwtsw\" (large genomes)\n");
  tmap_file_fprintf(tmap_file_stderr, "                     \t\"is\" (short genomes)\n");
  tmap_file_fprintf(tmap_file_stderr, "         --version   print the index format that will be created and exit\n");
  tmap_file_fprintf(tmap_file_stderr, "         -v          print verbose progress information\n");
  tmap_file_fprintf(tmap_file_stderr, "         -h          print this message\n");
  tmap_file_fprintf(tmap_file_stderr, "\n");
  return 1;
}

int tmap_index(int argc, char *argv[])
{
  int c;
  tmap_index_opt_t opt;

  opt.fn_fasta = NULL;
  opt.occ_interval = TMAP_BWT_OCC_INTERVAL; 
  opt.hash_width = TMAP_BWT_HASH_WIDTH;
  opt.sa_interval = TMAP_SA_INTERVAL; 
  opt.is_large = -1;
      
  if(2 == argc && 0 == strcmp("--version", argv[1])) {
      tmap_file_stdout = tmap_file_fdopen(fileno(stdout), "wb", TMAP_FILE_NO_COMPRESSION);
      tmap_file_fprintf(tmap_file_stdout, "%s\n", tmap_refseq_get_version_format(PACKAGE_VERSION));
      tmap_file_fclose(tmap_file_stdout);
      return 1;
  }

  while((c = getopt(argc, argv, "f:o:i:w:a:hv")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = tmap_strdup(optarg); break;
        case 'o':
          opt.occ_interval = atoi(optarg); break;
        case 'i':
          opt.sa_interval = atoi(optarg); break;
        case 'w':
          opt.hash_width = atoi(optarg); break;
        case 'a':
          if(0 == strcmp("is", optarg)) opt.is_large = 0;
          else if(0 == strcmp("bwtsw", optarg)) opt.is_large = 1;
          else tmap_error("Option -a value not correct", Exit, CommandLineArgument); 
          break; 
        case 'v':
          tmap_progress_set_verbosity(1); break;
        case 'h':
        default:
          return usage(&opt);
      }
  }

  if(argc != optind || 1 == argc) {
      return usage(&opt);
  }
  if(NULL == opt.fn_fasta) {
      tmap_error("required option -f", Exit, CommandLineArgument);
  }
  if(0 < opt.occ_interval && 0 != (opt.occ_interval % 16)) {
      tmap_error("option -o out of range", Exit, CommandLineArgument);
  }
  if(opt.hash_width <= 0) {
      tmap_error("option -w out of range", Exit, CommandLineArgument);
  }
  if(0 < opt.sa_interval && 0 != (opt.sa_interval % 2)) {
      tmap_error("option -i out of range", Exit, CommandLineArgument);
  }

  tmap_index_core(&opt);

  free(opt.fn_fasta);

  return 0;
}
