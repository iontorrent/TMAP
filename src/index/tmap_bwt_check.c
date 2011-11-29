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
#include "../io/tmap_file.h"
#include "tmap_bwt.h"
#include "tmap_bwt_match.h"


void 
tmap_bwa_check_core(const char *fn_fasta, int32_t length, int32_t print_sa, int32_t use_hash)
{
  tmap_bwt_t *bwt;
  uint8_t *seqs[2] = {NULL,NULL};
  int32_t i, n[2], hash_width = 0, asymmetric = 0, k, l;
  uint64_t hash_j;
  int64_t sum, j;
  tmap_bwt_match_occ_t sa;

  bwt = tmap_bwt_read(fn_fasta);

  if(0 == use_hash) {
      hash_width = bwt->hash_width;
      bwt->hash_width = 0;
  }
  
  /*
  seqs[0] = tmap_calloc(2, sizeof(uint8_t), "seqs[0]");
  seqs[0][0] = 0;
  seqs[0][1] = 1;
  n = tmap_bwt_match_exact(bwt, 2, seqs[0], &sa);
  int k;
  for(k=0;k<2;k++) {
      fputc("ACGTN"[seqs[0][k]], tmap_file_stderr);
  }
  if(0 < n && sa.k <= sa.l) {
      tmap_file_fprintf(tmap_file_stderr, " %llu %llu %d\n", sa.k, sa.l, n);
  }
  else {
      tmap_file_fprintf(tmap_file_stderr, " NA NA NA\n");
  }
  free(seqs[0]);
  */

  for(i=1;i<=length;i++) {
      seqs[0] = tmap_calloc(i, sizeof(uint8_t), "seqs[0]");
      seqs[1] = tmap_calloc(i, sizeof(uint8_t), "seqs[1]");
      for(j=0;j<i;j++) {
          seqs[1][j] = 3;
      }

      j = 0;
      hash_j = sum = 0;
      while(1) {
          if(i == j) {
              for(k=0;k<i;k++) {
                  seqs[1][k] = 3 - seqs[0][i-k-1];
              }
              for(k=0;k<2;k++) {
                  n[k] = tmap_bwt_match_exact(bwt, i, seqs[k], &sa);
                  if(0 == k) {
                      if(0 < n[k] && sa.k <= sa.l) {
                          sum += n[k];
                      }
                      if(1 == print_sa) {
                          for(l=0;l<i;l++) {
                              tmap_file_fprintf(tmap_file_stderr, "%c", "ACGTN"[seqs[k][l]]);
                          }
                          if(0 < n[k] && sa.k <= sa.l) {
                              tmap_file_fprintf(tmap_file_stderr, "\t%llu\t%llu\t%d\n", sa.k, sa.l, n[k]);
                          }
                          else {
                              tmap_file_fprintf(tmap_file_stderr, "\tNA\tNA\tNA\n");
                          }
                      }
                  }
              }
              if(0 == asymmetric && n[0] != n[1]) {
                  asymmetric = 1;
                  tmap_error("Asymmetriy found", Warn, OutOfRange);
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

      tmap_file_fprintf(tmap_file_stderr, "kmer:%d sum=%llu expected sum:%llu difference:%lld matched:%s\n",
              i,
              sum,
              bwt->seq_len - i + 1,
              sum - (bwt->seq_len - i + 1),
              sum == (bwt->seq_len - i + 1) ? "true" : "false");
  }

  if(0 == use_hash) {
      bwt->hash_width = hash_width;
  }
  tmap_bwt_destroy(bwt);
}

int 
tmap_bwt_check(int argc, char *argv[])
{
	int c, length = 12, print_sa = 0, use_hash = 0;

	while ((c = getopt(argc, argv, "l:pH")) >= 0) {
		switch (c) {
                  case 'l': length = atoi(optarg); break;
                  case 'p': print_sa = 1; break;
                  case 'H': use_hash = 1; break;
                  default: return 1;
		}
	}

	if (optind + 1 > argc) {
                tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [options] <in.fasta>\n", PACKAGE, argv[0]);
		tmap_file_fprintf(tmap_file_stderr, "Options:\n");
		tmap_file_fprintf(tmap_file_stderr, "         -l INT    the kmer length to check\n");
		tmap_file_fprintf(tmap_file_stderr, "         -p        print out the SA intervals for each kmer\n");
		tmap_file_fprintf(tmap_file_stderr, "         -H        use the hash to compute the SA intervals\n");
		tmap_file_fprintf(tmap_file_stderr, "\n");
		return 1;
	}
	tmap_bwa_check_core(argv[optind], length, print_sa, use_hash);
	return 0;
}
