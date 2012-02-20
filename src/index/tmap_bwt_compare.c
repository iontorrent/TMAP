#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <config.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../io/tmap_file.h"
#include "tmap_bwt.h"
#include "tmap_bwt_match.h"
#include "tmap_bwt_compare.h"

void 
tmap_bwt_compare_core2(tmap_bwt_t *bwt[2], int32_t length, int32_t print_msg, int32_t warn)
{
  uint8_t *seqs[2] = {NULL,NULL};
  char *str = NULL;
  int32_t i, asymmetric, k, m;
  uint64_t hash_j;
  int64_t sum, j;
  tmap_bwt_match_occ_t sa[2];
  tmap_bwt_int_t n[2][2];

  for(i=1;i<=length;i++) {
      seqs[0] = tmap_calloc(i, sizeof(uint8_t), "seqs[0]");
      seqs[1] = tmap_calloc(i, sizeof(uint8_t), "seqs[1]");
      str = tmap_calloc(i+1, sizeof(char), "str");
      for(j=0;j<i;j++) {
          seqs[1][j] = 3;
      }

      asymmetric = 0;
      j = 0;
      hash_j = sum = 0;
      while(1) {
          if(i == j) {
              for(k=0;k<i;k++) {
                  seqs[1][k] = 3 - seqs[0][i-k-1];
              }
              for(k=0;k<2;k++) {
                  for(m=0;m<2;m++) {
                      n[m][k] = tmap_bwt_match_exact_reverse(bwt[m], i, seqs[k], &sa[m]);
                  }
                  if(n[0][k] != n[1][k] || sa[0].k != sa[1].k || sa[0].l != sa[1].l || sa[0].hi != sa[1].hi || sa[0].offset != sa[1].offset) {
                      tmap_progress_print2("BWTs did not match");
                      tmap_progress_print2("n=[%llu,%llu]", n[0][k], n[1][k]);
                      tmap_progress_print2("k=[%llu,%llu]", sa[0].k, sa[1].k);
                      tmap_progress_print2("l=[%llu,%llu]", sa[0].l, sa[1].l);
                      tmap_progress_print2("hi=[%llu,%llu]", sa[0].hi, sa[1].hi);
                      tmap_progress_print2("offset=[%llu,%llu]", sa[0].offset, sa[1].offset);
                      tmap_bug();
                  }
              }
              for(k=0;k<2;k++) {
                  // use m == 0 && k = 0
                  if(0 == k) {
                      if(0 < n[0][k] && TMAP_BWT_INT_MAX != sa[0].k && sa[0].k <= sa[0].l) {
                          sum += n[0][k];
                      }
                  }
              }
              if(0 == asymmetric && n[0][0] != n[0][1]) {
                  asymmetric = 1;
                  //fprintf(stderr, "n[0][0]=%u n[0][1]=%u\n", n[0][0], n[0][1]);
                  tmap_error("Asymmetry found", Warn, OutOfRange);
              }

              j--;
              while(0 <= j && 3 == seqs[0][j]) {
                  seqs[0][j] = 0;
                  hash_j >>= 2;
                  j--;
              }
              if(j < 0) break;
              seqs[0][j]++;
              hash_j++;
              j++;
          }
          else {
              hash_j <<= 2;
              j++;
          }
      }

      free(seqs[0]);
      free(seqs[1]);
      free(str);

      j = (sum == (bwt[0]->seq_len - i + 1)) ? 0 : 1; // j==1 on fail
      if(1 == print_msg) {
          if(0 == j) tmap_progress_print2("%d-mer validation passed", i);
          else tmap_progress_print2("%d-mer validation failed: observed (%llu) != expected (%llu)\n", i, sum, bwt[0]->seq_len - i + 1);
      }
      if(0 == warn && 1 == j) {
          tmap_error("inconsistency found in the BWT", Exit, OutOfRange);
      }
  }
}

void 
tmap_bwt_compare_core(char *fn_fasta[], int32_t length, int32_t use_hash)
{
  tmap_bwt_t *bwt[2];
  int32_t i, hash_width[2] = {0,0};

  for(i=0;i<2;i++) {
      bwt[i] = tmap_bwt_read(fn_fasta[i]);
      if(0 == use_hash) {
          hash_width[i] = bwt[i]->hash_width;
          bwt[i]->hash_width = 0;
      }
  }

  tmap_bwt_compare_core2(bwt, length, 1, 1);

  for(i=0;i<2;i++) {
      if(0 == use_hash) {
          bwt[i]->hash_width = hash_width[i];
      }
      tmap_bwt_destroy(bwt[i]);
  }
}

int 
tmap_bwt_compare(int argc, char *argv[])
{
  int c, length = 12, use_hash = 0;

  tmap_progress_set_verbosity(1); 

  while ((c = getopt(argc, argv, "l:pH")) >= 0) {
      switch (c) {
        case 'l': length = atoi(optarg); break;
        case 'H': use_hash = 1; break;
        default: return 1;
      }
  }

  if (optind + 2 != argc) {
      tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [options] <in1.fasta> < in2.fasta>\n", PACKAGE, argv[0]);
      tmap_file_fprintf(tmap_file_stderr, "Options:\n");
      tmap_file_fprintf(tmap_file_stderr, "         -l INT    the kmer length to compare\n");
      tmap_file_fprintf(tmap_file_stderr, "         -H        use the hash to compute the SA intervals\n");
      tmap_file_fprintf(tmap_file_stderr, "\n");
      return 1;
  }

  tmap_bwt_compare_core(argv + optind, length, use_hash);

  return 0;
}
