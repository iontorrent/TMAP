/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_definitions.h"
#include "../server/tmap_shm.h"
#include "tmap_refseq.h"
#include "tmap_bwt_gen.h"
#include "tmap_bwt.h"
#include "tmap_sa.h"
#include "tmap_index.h"
#include "tmap_index_size.h"

static size_t
tmap_index_size_get(const char *fn_fasta)
{
  size_t s = 0;

  s += tmap_refseq_shm_read_num_bytes(fn_fasta);
  s += tmap_bwt_shm_read_num_bytes(fn_fasta);
  s += tmap_sa_shm_read_num_bytes(fn_fasta);

  return s;
}

static int 
usage(tmap_index_size_opt_t *opt)
{
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Usage: %s indexsize [options]", PACKAGE);
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Options (required):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -f FILE     the FASTA file name to index\n");
  tmap_file_fprintf(tmap_file_stderr, "         -h          print this message\n");
  tmap_file_fprintf(tmap_file_stderr, "\n");
  return 1;
}

int 
tmap_index_size(int argc, char *argv[])
{
  int c;
  tmap_index_size_opt_t opt;

  opt.fn_fasta = NULL;
      
  while((c = getopt(argc, argv, "f:h")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = tmap_strdup(optarg); break;
        case 'h':
          return usage(&opt);
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

  tmap_progress_set_verbosity(1); 
  tmap_progress_print2("index size (in bytes) = %llu", (unsigned long long int)tmap_index_size_get(opt.fn_fasta));

  free(opt.fn_fasta);

  return 0;
}
