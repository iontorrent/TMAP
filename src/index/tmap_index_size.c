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
tmap_index_size_print(const char *fn_fasta, uint64_t len, int32_t occ_interval, int32_t hash_width, int32_t sa_interval)
{
  size_t s, x, y, z;

  if(NULL != fn_fasta) { // exact
      x = tmap_refseq_shm_read_num_bytes(fn_fasta);
      y = tmap_bwt_shm_read_num_bytes(fn_fasta);
      z = tmap_sa_shm_read_num_bytes(fn_fasta);
  }
  else { // approx
      tmap_progress_print("the sizes will be approximated");
      x = tmap_refseq_approx_num_bytes(len);
      y = tmap_bwt_approx_num_bytes(2*len, occ_interval, hash_width);
      z = tmap_sa_approx_num_bytes(2*len, sa_interval);
  }


  s = x + y + z;

  tmap_progress_print2("reference size (in bytes) = %llu", (unsigned long long int)x);
  tmap_progress_print2("bwt size (in bytes) = %llu", (unsigned long long int)y);
  tmap_progress_print2("sa size (in bytes) = %llu", (unsigned long long int)z);
  tmap_progress_print2("total index size (in bytes) = %llu", (unsigned long long int)s);

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
  tmap_file_fprintf(tmap_file_stderr, "         -o INT      the occurrence interval (use %d, %d, %d, ...) [%d]\n",
                    TMAP_BWT_OCC_MOD, TMAP_BWT_OCC_MOD*2, TMAP_BWT_OCC_MOD*3, opt->occ_interval);
  tmap_file_fprintf(tmap_file_stderr, "         -w INT      the k-mer occurrence hash width [%d]\n", opt->hash_width);
  tmap_file_fprintf(tmap_file_stderr, "         -i INT      the suffix array interval (use 1, 2, 4, ...)[%d]\n", opt->sa_interval);
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
  opt.len = 0;
  opt.occ_interval = TMAP_BWT_OCC_INTERVAL;
  opt.hash_width = INT32_MAX;
  opt.sa_interval = TMAP_SA_INTERVAL;
      
  while((c = getopt(argc, argv, "f:l:o:i:w:h")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = tmap_strdup(optarg); break;
        case 'l':
          opt.len = atol(optarg); break;
        case 'o':
          opt.occ_interval = atoi(optarg); break;
        case 'i':
          opt.sa_interval = atoi(optarg); break;
        case 'w':
          opt.hash_width = atoi(optarg); break;
        case 'h':
          return usage(&opt);
        default:
          return usage(&opt);
      }
  }

  if(argc != optind || 1 == argc) {
      return usage(&opt);
  }
  if(NULL == opt.fn_fasta && 0 == opt.len) {
      tmap_error("option -f or option -l is required", Exit, CommandLineArgument);
  }

  tmap_progress_set_verbosity(1); 
  tmap_index_size_print(opt.fn_fasta, opt.len, opt.occ_interval, opt.hash_width, opt.sa_interval);

  free(opt.fn_fasta);

  return 0;
}
