/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/shm.h>
#include <time.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_rand.h"
#include "../io/tmap_file.h"
#include "tmap_refseq.h"
#include "tmap_bwt_gen.h"
#include "tmap_bwt.h"
#include "tmap_sa.h"
#include "tmap_index.h"
#include "tmap_bwt_match.h"
#include "tmap_index_speed.h"

static int32_t  
tmap_index_speed_test(tmap_index_t *index, clock_t *total_clock, tmap_index_speed_opt_t *opt)
{
  tmap_rand_t *rand;
  int32_t i, j, num_found, found;
  tmap_bwt_int_t k, l;
  tmap_bwt_int_t pacpos;
  uint8_t *seq = NULL;
  tmap_bwt_match_occ_t cur;
  clock_t start_clock = 0;
  tmap_bwt_int_t cntk[4], cntl[4], ok, ol;

  rand = tmap_rand_init(13);
  if(0 == opt->func) seq = tmap_malloc(opt->kmer_length * sizeof(uint8_t), "seq");

  num_found = 0;
  for(i=0;i<opt->kmer_num;i++) {
      if(0 < i && 0 == (i % 1000000)) {
          tmap_progress_print2("processed %d kmers", i);
      }
      if(0 == opt->func) {
          if(tmap_rand_get(rand) < opt->rand_frac) {
              for(j=0;j<opt->kmer_length;j++) {
                  seq[j] = (uint8_t)(tmap_rand_get(rand) * 4);
              }
          }
          else {
              // get a position (one-based)
              pacpos = 1 + (tmap_bwt_int_t)(tmap_rand_get(rand) * (index->refseq->len - opt->kmer_length + 1));
              if(0 == tmap_refseq_subseq(index->refseq, pacpos, opt->kmer_length, seq)) {
                  i--;
                  continue;
              }
          }
          found = 0;
          start_clock = clock();
          if(0 < tmap_bwt_match_exact_reverse(index->bwt, opt->kmer_length, seq, &cur)) {
              if(0 <= opt->enum_max_hits && (cur.l - cur.k + 1) <= opt->enum_max_hits) {
                  for(k=cur.k;k<=cur.l;k++) {
                      // retrieve the packed position
                      pacpos = index->bwt->seq_len - tmap_sa_pac_pos(index->sa, index->bwt, k) - opt->kmer_length + 1;
                  }
              }
              found = 1;
          }
          (*total_clock) += clock() - start_clock;
          if(0 < found) {
              num_found++;
          }
      }
      else {
          k = tmap_rand_int(rand) % index->bwt->seq_len;
          l = k; //l = tmap_rand_int(rand) % index->bwt->seq_len;
          j = tmap_rand_int(rand) % 4;
          start_clock = clock();
          switch(opt->func) {
            case 1:
              k = tmap_bwt_occ(index->bwt, k, j);
              break;
            case 2:
              tmap_bwt_occ4(index->bwt, k, cntk);
              break;
            case 3:
              tmap_bwt_2occ(index->bwt, k, l, j, &ok, &ol);
              break;
            case 4:
              tmap_bwt_2occ4(index->bwt, k, l, cntk, cntl);
              break;
            case 5:
              k = tmap_sa_pac_pos(index->sa, index->bwt, k);
              break;
            default:
              tmap_bug();
          }
          (*total_clock) += clock() - start_clock;
      }
  }
  tmap_progress_print2("processed %d kmers", i);

  tmap_rand_destroy(rand);
  if(0 == opt->func) free(seq);

  return num_found;
}

static void 
tmap_index_speed_core(tmap_index_speed_opt_t *opt)
{
  tmap_index_t *index = NULL;
  clock_t total_clock= 0;
  time_t start_time, end_time;
  int32_t num_found;

  // read in the index
  index = tmap_index_init(opt->fn_fasta, opt->shm_key);

  // modify the hash width
  if(0 <= opt->hash_width) {
      if(index->bwt->hash_width < opt->hash_width) {
          tmap_error("the hash width too large (-w)", Exit, CommandLineArgument);
      }
      index->bwt->hash_width = opt->hash_width;
  }

  // clock on
  start_time = time(NULL);

  // run the speed test
  num_found = tmap_index_speed_test(index, &total_clock, opt);

  // clock off
  end_time = time(NULL);

  // free
  tmap_index_destroy(index);

  // print the results
  tmap_progress_print2("[tmap index] found %d out of %d (%.2f%%)", num_found, opt->kmer_num, (100.0 * num_found) / opt->kmer_num);
  tmap_progress_print2("[tmap index] wallclock time: %d seconds", (int)difftime(end_time,start_time));
  tmap_progress_print2("[tmap index] index lookup cpu cycle time: %.2f seconds", (float)(total_clock) / CLOCKS_PER_SEC);
}

static int 
usage(tmap_index_speed_opt_t *opt)
{
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Usage: %s index [options]", PACKAGE);
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Options (required):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -f FILE     the FASTA file name to test [%s]\n", 
                    (NULL == opt->fn_fasta) ? "not using" : opt->fn_fasta);
  tmap_file_fprintf(tmap_file_stderr, "         -k INT      use shared memory with the following key [%d]\n", opt->shm_key);
  tmap_file_fprintf(tmap_file_stderr, "Options (optional):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -w INT      the new index hash width [%d]\n", opt->hash_width);
  tmap_file_fprintf(tmap_file_stderr, "         -N INT      the number of kmers to simulate [%d]\n", opt->kmer_num);
  tmap_file_fprintf(tmap_file_stderr, "         -F INT      the function to test [%d]:\n", opt->func);
  tmap_file_fprintf(tmap_file_stderr, "                         0 - k-mer lookups and position retrieval\n");
  tmap_file_fprintf(tmap_file_stderr, "                         1 - tmap_bwt_occ\n");
  tmap_file_fprintf(tmap_file_stderr, "                         2 - tmap_bwt_occ4\n");
  tmap_file_fprintf(tmap_file_stderr, "                         3 - tmap_bwt_2occ\n");
  tmap_file_fprintf(tmap_file_stderr, "                         4 - tmap_bwt_2occ4\n");
  tmap_file_fprintf(tmap_file_stderr, "                         5 - tmap_sa_pac_pos\n");
  tmap_file_fprintf(tmap_file_stderr, "         -e INT      the maximum number of hits to enumerate with -F 0 (-1 for unlimited, 0 to disable) [%d]\n", opt->enum_max_hits);
  tmap_file_fprintf(tmap_file_stderr, "         -K INT      the kmer length to simulate with -F 0 [%d]\n", opt->kmer_length);
  tmap_file_fprintf(tmap_file_stderr, "         -R FLOAT    the fraction of random kmers with -F 0 [%.2lf]\n", opt->rand_frac);
  tmap_file_fprintf(tmap_file_stderr, "         -v          print verbose progress information\n");
  tmap_file_fprintf(tmap_file_stderr, "         -h          print this message\n");
  tmap_file_fprintf(tmap_file_stderr, "\n");
  return 1;
}

int 
tmap_index_speed(int argc, char *argv[])
{
  int c;
  tmap_index_speed_opt_t opt;

  opt.fn_fasta = NULL;
  opt.hash_width = -1;
  opt.enum_max_hits = 1024;
  opt.kmer_length = 12;
  opt.kmer_num = 100000;
  opt.rand_frac = 0.0;
  opt.shm_key = 0;
  opt.func = 0;
      
  while((c = getopt(argc, argv, "f:k:w:e:K:N:R:F:hv")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = tmap_strdup(optarg); break;
        case 'k':
          opt.shm_key = atoi(optarg); break;
        case 'w':
          opt.hash_width = atoi(optarg); break;
        case 'e':
          opt.enum_max_hits = atoi(optarg); break;
        case 'K':
          opt.kmer_length = atoi(optarg); break;
        case 'N':
          opt.kmer_num = atoi(optarg); break;
        case 'R':
          opt.rand_frac = atof(optarg); break;
        case 'v':
          tmap_progress_set_verbosity(1); break;
        case 'F':
          opt.func = atoi(optarg); break;
        case 'h':
        default:
          return usage(&opt);
      }
  }

  if(argc != optind || 1 == argc) {
      return usage(&opt);
  }
  if(NULL == opt.fn_fasta && 0 == opt.shm_key) {
      tmap_error("required option -f or -k", Exit, CommandLineArgument);
  }
  else if(NULL != opt.fn_fasta && 0 != opt.shm_key) {
      tmap_error("both option -f or -k cannot be used", Exit, CommandLineArgument);
  }
  if(opt.kmer_length <= 0) {
      tmap_error("the kmer length to simulate must be greater than zero (-K)", Exit, CommandLineArgument);
  }
  if(opt.kmer_num <= 0) {
      tmap_error("the numer of kmers to simulate must be greater than zero (-N)", Exit, CommandLineArgument);
  }
  if(opt.rand_frac < 0 || 1 < opt.rand_frac) {
      tmap_error("the option -R must be between 0 and 1", Exit, CommandLineArgument);
  }
  if(opt.func < 0 || 5 < opt.func) {
      tmap_error("the option -F must be between 0 and 5", Exit, CommandLineArgument);
  }

  tmap_index_speed_core(&opt);

  free(opt.fn_fasta);

  return 0;
}
