#include <stdlib.h>
#include <stdint.h>

#include "../util/fmap_definitions.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_error.h"
#include "../sw/fmap_sw.h"
#include "fmap_sam2fs_aux.h"

#define FMAP_SAM2FS_AUX_GAP_O 1
#define FMAP_SAM2FS_AUX_GAP_E 5

typedef struct {
    int32_t *qseq;
    int32_t *tseq;
    char *aln;
    int32_t len, mem;
} fmap_sam2fs_aux_aln_t;

static fmap_sam2fs_aux_aln_t * 
fmap_sam2fs_aux_aln_init()
{
  return fmap_calloc(1, sizeof(fmap_sam2fs_aux_aln_t), "return");
}

static void
fmap_sam2fs_aux_aln_destroy(fmap_sam2fs_aux_aln_t *aln)
{
  free(aln->qseq);
  free(aln->tseq);
  free(aln->aln);
  free(aln);
}

static void
fmap_sam2fs_aux_aln_add(fmap_sam2fs_aux_aln_t *aln, int32_t read_n, int32_t ref_n)
{
  if(aln->mem <= aln->len) {
      aln->mem = aln->len + 1;
      fmap_roundup32(aln->mem);
      aln->qseq = fmap_realloc(aln->qseq, sizeof(int32_t) * aln->mem, "aln->qseq");
      aln->tseq = fmap_realloc(aln->tseq, sizeof(int32_t) * aln->mem, "aln->tseq");
      aln->aln = fmap_realloc(aln->aln, sizeof(char) * aln->mem, "aln->aln");
  }
  aln->qseq[aln->len] = read_n;
  aln->tseq[aln->len] = ref_n;
  if(-1 == read_n) {
      aln->aln[aln->len] = '-';
  }
  else if(-1 == ref_n) {
      aln->aln[aln->len] = '+';
  }
  else if(read_n == ref_n) {
      aln->aln[aln->len] = '|';
  }
  else {
      aln->aln[aln->len] = ' ';
  }
  aln->len++;
}

enum {
    FMAP_SAM2FS_AUX_FROM_M = 0,
    FMAP_SAM2FS_AUX_FROM_I = 1,
    FMAP_SAM2FS_AUX_FROM_D = 2,
    FMAP_SAM2FS_AUX_FROM_S = 3
};

typedef struct {
    int32_t match_score;
    int32_t match_from;
    int32_t ins_score;
    int32_t ins_from;
    int32_t del_score;
    int32_t del_from;
} fmap_sam2fs_aux_dpcell_t;

typedef struct {
    uint8_t *flow;
    int32_t l, m;
} fmap_sam2fs_aux_flow_t;

static fmap_sam2fs_aux_flow_t *
fmap_sam2fs_aux_flow_init()
{
  return fmap_calloc(1, sizeof(fmap_sam2fs_aux_flow_t), "return");
}

static void
fmap_sam2fs_aux_flow_convert(fmap_sam2fs_aux_flow_t *a, uint8_t *seq, int32_t len, uint8_t *flow_order)
{
  int32_t i, k, l;
  
  a->l = 0;
  a->m = 0;

  i=k=0;
  while(i < len) {
      // empty flow
      while(flow_order[k] != seq[i]) {
          if(a->m <= a->l) {
              a->m = a->l + 1;
              fmap_roundup32(a->m);
              a->flow = fmap_realloc(a->flow, a->m * sizeof(uint8_t), "a->flow");
          }
          a->flow[a->l] = 0;
          a->l++;
          k = (k+1) & 3; // % 4
      }
      // current flow
      l = i+1;
      while(l < len && flow_order[k] == seq[l]) {
          l++;
      }
      if(a->m <= a->l) {
          a->m = a->l + 1;
          fmap_roundup32(a->m);
          a->flow = fmap_realloc(a->flow, a->m * sizeof(uint8_t), "a->flow");
      }
      a->flow[a->l] = l-i;
      a->l++;
      k = (k+1) & 3; // % 4
      i = l;
  }
}

static void
fmap_sam2fs_aux_flow_destroy(fmap_sam2fs_aux_flow_t *a)
{
  free(a->flow);
  free(a);
}

void
fmap_sam2fs_aux_flow_align(fmap_file_t *fp, uint8_t *qseq, int32_t qseq_len, uint8_t *tseq, int32_t tseq_len, uint8_t *flow_order)
{
  int32_t i, j, k;
  int32_t ctype;
  int32_t best_i, best_j, best_score, best_ctype;
  fmap_sam2fs_aux_flow_t *f_qseq=NULL, *f_tseq=NULL;
  fmap_sam2fs_aux_dpcell_t **dp=NULL;
  fmap_sam2fs_aux_aln_t *aln;

  // init
  f_qseq = fmap_sam2fs_aux_flow_init();
  f_tseq = fmap_sam2fs_aux_flow_init();

  // convert bases to flow space
  fmap_sam2fs_aux_flow_convert(f_qseq, qseq, qseq_len, flow_order);
  fmap_sam2fs_aux_flow_convert(f_tseq, tseq, tseq_len, flow_order);
  
  /*
  // target
  for(i=0;i<tseq_len;i++) {
      fputc("ACGT"[tseq[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<f_tseq->l;i++) {
      if(0 < i) fputc(',', stderr);
      fprintf(stderr, "%d", f_tseq->flow[i]);
  }
  fputc('\n', stderr);
  // query
  for(i=0;i<qseq_len;i++) {
      fputc("ACGT"[qseq[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<f_qseq->l;i++) {
      if(0 < i) fputc(',', stderr);
      fprintf(stderr, "%d", f_qseq->flow[i]);
  }
  fputc('\n', stderr);
  */

  // init dp matrix
  dp = fmap_malloc(sizeof(fmap_sam2fs_aux_dpcell_t*)*(1+f_qseq->l), "dp");
  for(i=0;i<=f_qseq->l;i++) {
      dp[i] = fmap_malloc(sizeof(fmap_sam2fs_aux_dpcell_t)*(1+f_tseq->l), "dp[i]");
  }
  for(i=0;i<=f_qseq->l;i++) {
      for(j=0;j<=f_tseq->l;j++) {
          dp[i][j].match_score = dp[i][j].ins_score = dp[i][j].del_score = FMAP_SW_MINOR_INF;
          dp[i][j].match_from = dp[i][j].ins_from = dp[i][j].del_from = FMAP_SAM2FS_AUX_FROM_S;
      }
  }

  // init start cells
  dp[0][0].match_score = 0; // global
  for(i=k=1;i<=f_qseq->l;i++) {
      dp[i][0].ins_score = 0 - FMAP_SAM2FS_AUX_GAP_O - (k-1)*FMAP_SAM2FS_AUX_GAP_E;
      dp[i][0].ins_from = ((k==1) ? FMAP_SAM2FS_AUX_FROM_M : FMAP_SAM2FS_AUX_FROM_I);
      if(0 == (i%4)) {
          k++;
      }
  }
  for(j=k=1;j<=f_tseq->l;j++) {
      dp[0][j].del_score = 0 - FMAP_SAM2FS_AUX_GAP_O - (k-1)*FMAP_SAM2FS_AUX_GAP_E;
      dp[0][j].del_from = ((1 == k) ? FMAP_SAM2FS_AUX_FROM_M : FMAP_SAM2FS_AUX_FROM_D);
      if(0 == (j%4)) {
          k++;
      }
  }
  // align
  for(i=1;i<=f_qseq->l;i++) {
      int32_t i_from = ((i<4) ? 0 : (i-4)); // don't go before the first row
      for(j=1;j<=f_tseq->l;j++) {
          int32_t j_from = ((j<4) ? 0 : (j-4)); // don't go before the first column

          // horizontal
          if(dp[i][j_from].del_score - FMAP_SAM2FS_AUX_GAP_E < dp[i][j_from].match_score - FMAP_SAM2FS_AUX_GAP_O) {
              dp[i][j].del_score = dp[i][j_from].match_score - FMAP_SAM2FS_AUX_GAP_O;
              dp[i][j].del_from = FMAP_SAM2FS_AUX_FROM_M;
          }
          else {
              dp[i][j].del_score = dp[i][j_from].del_score - FMAP_SAM2FS_AUX_GAP_E;
              dp[i][j].del_from = FMAP_SAM2FS_AUX_FROM_D;
          }

          // vertical
          if(dp[i_from][j].ins_score - FMAP_SAM2FS_AUX_GAP_E < dp[i_from][j].match_score - FMAP_SAM2FS_AUX_GAP_O) {
              dp[i][j].ins_score = dp[i_from][j].match_score - FMAP_SAM2FS_AUX_GAP_O;
              dp[i][j].ins_from = FMAP_SAM2FS_AUX_FROM_M;
          }
          else {
              dp[i][j].ins_score = dp[i_from][j].ins_score - FMAP_SAM2FS_AUX_GAP_E;
              dp[i][j].ins_from = FMAP_SAM2FS_AUX_FROM_I;
          }

          // diagonal
          int32_t s = ((f_qseq->flow[i-1] < f_tseq->flow[j-1]) ? (f_tseq->flow[j-1]-f_qseq->flow[i-1]) : (f_qseq->flow[i-1]-f_tseq->flow[j-1]));
          if(dp[i-1][j-1].ins_score < dp[i-1][j-1].match_score) {
              if(dp[i-1][j-1].del_score < dp[i-1][j-1].match_score) {
                  dp[i][j].match_score = dp[i-1][j-1].match_score - s;
                  dp[i][j].match_from = FMAP_SAM2FS_AUX_FROM_M;
              }
              else {
                  dp[i][j].match_score = dp[i-1][j-1].del_score - s;
                  dp[i][j].match_from = FMAP_SAM2FS_AUX_FROM_D;
              }
          }
          else {
              if(dp[i-1][j-1].del_score < dp[i-1][j-1].ins_score) {
                  dp[i][j].match_score = dp[i-1][j-1].ins_score - s;
                  dp[i][j].match_from = FMAP_SAM2FS_AUX_FROM_I;
              }
              else {
                  dp[i][j].match_score = dp[i-1][j-1].del_score - s;
                  dp[i][j].match_from = FMAP_SAM2FS_AUX_FROM_D;
              }
          }
      }
  }

  // Get best scoring cell
  best_score = FMAP_SW_MINOR_INF; 
  best_ctype = FMAP_SAM2FS_AUX_FROM_S;
  best_i = -1;
  best_j = -1;

  // last four flows in the last row
  for(k=0;k<4;k++) {
      i = f_qseq->l - k;
      if(i <= 0) {
          break;
      }
      if(best_score < dp[i][f_tseq->l].del_score) {
          best_i = i;
          best_j = f_tseq->l;
          best_score = dp[i][f_tseq->l].del_score;
          best_ctype = FMAP_SAM2FS_AUX_FROM_D;
      }
      if(best_score < dp[i][f_tseq->l].ins_score) {
          best_i = i;
          best_j = f_tseq->l;
          best_score = dp[i][f_tseq->l].ins_score;
          best_ctype = FMAP_SAM2FS_AUX_FROM_I;
      }
      if(best_score < dp[i][f_tseq->l].match_score) {
          best_i = i;
          best_j = f_tseq->l;
          best_score = dp[i][f_tseq->l].match_score;
          best_ctype = FMAP_SAM2FS_AUX_FROM_M;
      }
  }
  // last col, last four flow
  for(k=0;k<4;k++) {
      j = f_tseq->l - k;
      if(j <= 0) {
          break;
      }
      if(best_score < dp[f_qseq->l][j].del_score) {
          best_i = f_qseq->l;
          best_j = j;
          best_score = dp[f_qseq->l][j].del_score;
          best_ctype = FMAP_SAM2FS_AUX_FROM_D;
      }
      if(best_score < dp[f_qseq->l][j].ins_score) {
          best_i = f_qseq->l;
          best_j = j;
          best_score = dp[f_qseq->l][j].ins_score;
          best_ctype = FMAP_SAM2FS_AUX_FROM_I;
      }
      if(best_score < dp[f_qseq->l][j].match_score) {
          best_i = f_qseq->l;
          best_j = j;
          best_score = dp[f_qseq->l][j].match_score;
          best_ctype = FMAP_SAM2FS_AUX_FROM_M;
      }
  }

  i = best_i;
  j = best_j;
  ctype = best_ctype;
  
  aln = fmap_sam2fs_aux_aln_init();

  // add in alignment end padding
  for(k=f_qseq->l;i<k;k--) { // insertion
      fmap_sam2fs_aux_aln_add(aln, f_qseq->flow[k-1], -1);
  }
  for(k=f_tseq->l;j<k;k--) { // deletion
      fmap_sam2fs_aux_aln_add(aln, -1, f_tseq->flow[k-1]);
  }

  // trace path back
  while(0 < i || 0 < j) {
      int32_t next_ctype = -1;
      //fprintf(stderr, "i=%d j=%d ctype=%d\n", i, j, ctype);

      switch(ctype) {
        case FMAP_SAM2FS_AUX_FROM_M:
          next_ctype = dp[i][j].match_from;
          fmap_sam2fs_aux_aln_add(aln, f_qseq->flow[i-1], f_tseq->flow[j-1]);
          i--;
          j--;
          break;
        case FMAP_SAM2FS_AUX_FROM_I:
          next_ctype = dp[i][j].ins_from;
          for(k=0;k<4;k++) {
              fmap_sam2fs_aux_aln_add(aln, f_qseq->flow[i-1], -1);
              i--;
          }
          break;
        case FMAP_SAM2FS_AUX_FROM_D:
          next_ctype = dp[i][j].del_from;
          for(k=0;k<4;k++) {
              fmap_sam2fs_aux_aln_add(aln, -1, f_tseq->flow[j-1]);
              j--;
          }
          break;
        default:
          fmap_error(NULL, Exit, OutOfRange);
      }

      ctype = next_ctype;
  }

  // print
  for(i=aln->len-1;0<=i;i--) {
      if(0 <= aln->tseq[i]) fmap_file_fprintf(fp, "%d", aln->tseq[i]);
      else fmap_file_fprintf(fp, "-");
      if(0 < i) fmap_file_fprintf(fp, ","); 
  }
  fmap_file_fprintf(fp, "\t"); 
  for(i=aln->len-1;0<=i;i--) {
      fmap_file_fprintf(fp, "%c", aln->aln[i]);
      if(0 < i) fmap_file_fprintf(fp, ","); 
  }
  fmap_file_fprintf(fp, "\t"); 
  for(i=aln->len-1;0<=i;i--) {
      if(0 <= aln->qseq[i]) fmap_file_fprintf(fp, "%d", aln->qseq[i]);
      else fmap_file_fprintf(fp, "-");
      if(0 < i) fmap_file_fprintf(fp, ","); 
  }

  // destroy
  for(i=0;i<=f_qseq->l;i++) {
      free(dp[i]);
  }
  free(dp);
  fmap_sam2fs_aux_flow_destroy(f_qseq);
  fmap_sam2fs_aux_flow_destroy(f_tseq);
  fmap_sam2fs_aux_aln_destroy(aln);
}
