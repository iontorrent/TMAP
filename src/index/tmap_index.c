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

tmap_index_t*
tmap_index_init(const char *fn_fasta, key_t shm_key)
{
  tmap_index_t *index = NULL;

  index = tmap_calloc(1, sizeof(tmap_index_t), "index");

  index->shm_key = shm_key;

  // get the reference information
  if(0 == index->shm_key) {
      tmap_progress_print("reading in reference data");
      index->refseq = tmap_refseq_read(fn_fasta);
      index->bwt = tmap_bwt_read(fn_fasta);
      index->sa = tmap_sa_read(fn_fasta);
      tmap_progress_print2("reference data read in");
  }
  else {
      tmap_progress_print("retrieving reference data from shared memory");
      index->shm = tmap_shm_init(index->shm_key, 0, 0);
      if(NULL == (index->refseq = tmap_refseq_shm_unpack(tmap_shm_get_buffer(index->shm, TMAP_SHM_LISTING_REFSEQ)))) {
          tmap_error("the packed reference sequence was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (index->bwt = tmap_bwt_shm_unpack(tmap_shm_get_buffer(index->shm, TMAP_SHM_LISTING_BWT)))) {
          tmap_error("the BWT string was not found in shared memory", Exit, SharedMemoryListing);
      }
      if(NULL == (index->sa = tmap_sa_shm_unpack(tmap_shm_get_buffer(index->shm, TMAP_SHM_LISTING_SA)))) {
          tmap_error("the SA was not found in shared memory", Exit, SharedMemoryListing);
      }
      tmap_progress_print2("reference data retrieved from shared memory");
  }

  if((index->refseq->len << 1) != index->bwt->seq_len) {
      tmap_error("refseq and bwt lengths do not match", Exit, OutOfRange);
  }
  if((index->refseq->len << 1) != index->sa->seq_len) {
      tmap_error("refseq and sa lengths do not match", Exit, OutOfRange);
  }
  
  return index;
}

void
tmap_index_destroy(tmap_index_t *index)
{
  tmap_refseq_destroy(index->refseq);
  tmap_bwt_destroy(index->bwt);
  tmap_sa_destroy(index->sa);
  if(0 < index->shm_key) {
      tmap_shm_destroy(index->shm, 0);
  }
  free(index);
}

static void tmap_index_core(tmap_index_opt_t *opt)
{
  uint64_t ref_len = 0;

  // pack the reference sequence
  ref_len = tmap_refseq_fasta2pac(opt->fn_fasta, TMAP_FILE_NO_COMPRESSION, 0);
      
  if(TMAP_INDEX_TOO_BIG_GENOME <= ref_len) { // too big (2^32 - 1)!
      tmap_error("Reference sequence too large", Exit, OutOfRange);
  }

  // check returned genome size
  if(opt->is_large < 0) {
      if(TMAP_INDEX_LARGE_GENOME <= ref_len) { 
          opt->is_large = 1;
          tmap_progress_print("defaulting to \"bwtsw\" BWT construction algorithm");
      }
      else {
          opt->is_large = 0;
          tmap_progress_print("defaulting to \"is\" BWT construction algorithm");
      }
  }

  // create the bwt 
  tmap_bwt_pac2bwt(opt->fn_fasta, opt->is_large, opt->occ_interval, opt->hash_width);

  // create the suffix array
  tmap_sa_bwt2sa(opt->fn_fasta, opt->sa_interval);

  // pack the reference sequence
  ref_len = tmap_refseq_fasta2pac(opt->fn_fasta, TMAP_FILE_NO_COMPRESSION, 1);
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
      return 0;
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
  if(opt.occ_interval < 16 || 0 != (opt.occ_interval % 16)) {
      tmap_error("option -o out of range", Exit, CommandLineArgument);
  }
  if(opt.hash_width < 0) {
      tmap_error("option -w out of range", Exit, CommandLineArgument);
  }
  if(opt.sa_interval <= 0 || (1 < opt.sa_interval && 0 != (opt.sa_interval % 2))) {
      tmap_error("option -i out of range", Exit, CommandLineArgument);
  }

  tmap_index_core(&opt);

  free(opt.fn_fasta);

  return 0;
}
