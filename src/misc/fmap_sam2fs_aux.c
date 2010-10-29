#include <stdlib.h>
#include <stdint.h>

#include "../util/fmap_definitions.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_error.h"
#include "../sw/fmap_sw.h"
#include "fmap_sam2fs_aux.h"

// do not use
#define FMAP_SAM2FS_AUX_GAP_O 4
#define FMAP_SAM2FS_AUX_GAP_E 4

// added gap penalty
#define FMAP_SAM2FS_AUX_GAP 0

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
fmap_sam2fs_aux_flow_convert(fmap_sam2fs_aux_flow_t *a, uint8_t *seq, int32_t len, 
                             uint8_t *flow_order, uint8_t *qseq, int32_t qseq_len)
{
  int32_t i, k, l, next_i;

  i = k = a->l = a->m = 0;

  // move to the first flow
  if(NULL != qseq) {
      // we skip the flows not found in the qseq
      while('-' == qseq[i]) { // this could move beyond the end if all gaps
          if(qseq_len <= i) {
              fmap_error(NULL, Exit, OutOfRange);
          }
          i++;
      }
      while(flow_order[k] != qseq[i]) {
          k = (k+1) & 3; // % 4
      }
  }
  else {
      while('-' == seq[i]) { // this could move beyond the end if all gaps
          if(len <= i) {
              fmap_error(NULL, Exit, OutOfRange);
          }
          i++;
      }
      while(flow_order[k] != seq[i]) {
          k = (k+1) & 3; // % 4
      }
  }

  i=0;
  while(i < len) {
      int32_t before = i;

      // move beyond the initial gaps
      while(i < len && '-' == seq[i]) {
          i++;
      }
      if(len <= i) break;

      if(3 < seq[i]) {
          fprintf(stderr, "i=%d seq[i]=%d len=%d NULL=%d\n", i, seq[i], len, (NULL==qseq) ? 1 : 0);
          fmap_error(NULL, Exit, OutOfRange);
      }

      // skip over empty flow
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

      // get the number of bases in current flow
      next_i = i+1;
      l = 1;
      while(next_i < len 
            && (flow_order[k] == seq[next_i] || '-' == seq[next_i])) {
          if(flow_order[k] == seq[next_i]) {
              l++;
          }
          next_i++;
      }
      if(a->m <= a->l) { // more memory?
          a->m = a->l + 1;
          fmap_roundup32(a->m);
          a->flow = fmap_realloc(a->flow, a->m * sizeof(uint8_t), "a->flow");
      }
      a->flow[a->l] = l;
      a->l++;
      k = (k+1) & 3; // % 4
      i = next_i;
      if(i <= before) {
          fmap_error(NULL, Exit, OutOfRange);
      }
      // move beyond the gaps
      while(i < len && '-' == seq[i]) {
          i++;
      }
  }

  if(NULL != qseq) {
      k = (k+3) & 3; // move back one flow

      // get the last base in qseq
      i = qseq_len-1;
      while('-' == qseq[i]) {
          i--;
          if(i < 0) {
              fmap_error(NULL, Exit, OutOfRange);
          }
      }

      // add empty flows
      while(flow_order[k] != qseq[i]) {
          if(a->m <= a->l) {
              a->m = a->l + 1;
              fmap_roundup32(a->m);
              a->flow = fmap_realloc(a->flow, a->m * sizeof(uint8_t), "a->flow");
          }
          a->flow[a->l] = 0;
          a->l++;
          k = (k+1) & 3; // % 4
      }
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
  int32_t gap_sum_i, gap_sum_j;
  fmap_sam2fs_aux_flow_t *f_qseq=NULL, *f_tseq=NULL;
  fmap_sam2fs_aux_dpcell_t **dp=NULL;
  fmap_sam2fs_aux_aln_t *aln;

  // init
  f_qseq = fmap_sam2fs_aux_flow_init();
  f_tseq = fmap_sam2fs_aux_flow_init();

  // convert bases to flow space
  // NOTE: this will pad empty flows at the beginning and the end of the
  // reference
  fmap_sam2fs_aux_flow_convert(f_qseq, qseq, qseq_len, flow_order, NULL, 0);
  fmap_sam2fs_aux_flow_convert(f_tseq, tseq, tseq_len, flow_order, qseq, qseq_len);

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
  for(i=k=1,gap_sum_i=FMAP_SAM2FS_AUX_GAP;i<=f_qseq->l;i++) {
      if(4 < i) gap_sum_i -= f_qseq->flow[i-5]; // substract off
      gap_sum_i += f_qseq->flow[i-1]; // add current flow
      dp[i][0].ins_score = 0 - gap_sum_i;
      //dp[i][0].ins_score = 0 - FMAP_SAM2FS_AUX_GAP_O - (k-1)*FMAP_SAM2FS_AUX_GAP_E;
      dp[i][0].ins_from = ((k==1) ? FMAP_SAM2FS_AUX_FROM_M : FMAP_SAM2FS_AUX_FROM_I);
      if(0 == (i%4)) {
          k++;
      }
  }
  for(j=k=1,gap_sum_j=FMAP_SAM2FS_AUX_GAP;j<=f_tseq->l;j++) {
      if(4 < j) gap_sum_j -= f_tseq->flow[j-5]; // subtract off
      gap_sum_j += f_tseq->flow[j-1];
      dp[0][j].del_score = 0 - gap_sum_j;
      //dp[0][j].del_score = 0 - FMAP_SAM2FS_AUX_GAP_O - (k-1)*FMAP_SAM2FS_AUX_GAP_E;
      dp[0][j].del_from = ((1 == k) ? FMAP_SAM2FS_AUX_FROM_M : FMAP_SAM2FS_AUX_FROM_D);
      if(0 == (j%4)) {
          k++;
      }
  }
  // align
  gap_sum_i = gap_sum_j = FMAP_SAM2FS_AUX_GAP;
  for(i=1;i<=f_qseq->l;i++) {
      int32_t i_from = ((i<4) ? 0 : (i-4)); // don't go before the first row

      if(4 < i) gap_sum_i -= f_qseq->flow[i-5]; // substract off
      gap_sum_i += f_qseq->flow[i-1]; // add current flow

      for(j=1;j<=f_tseq->l;j++) {
          int32_t j_from = ((j<4) ? 0 : (j-4)); // don't go before the first column

          // gap sum
          if(4 < j) gap_sum_j -= f_tseq->flow[j-5]; // subtract off
          gap_sum_j += f_tseq->flow[j-1];

          // horizontal
          /*
             if(dp[i][j_from].del_score - FMAP_SAM2FS_AUX_GAP_E < dp[i][j_from].match_score - FMAP_SAM2FS_AUX_GAP_O) {
             dp[i][j].del_score = dp[i][j_from].match_score - FMAP_SAM2FS_AUX_GAP_O;
             dp[i][j].del_from = FMAP_SAM2FS_AUX_FROM_M;
             }
             else {
             dp[i][j].del_score = dp[i][j_from].del_score - FMAP_SAM2FS_AUX_GAP_E;
             dp[i][j].del_from = FMAP_SAM2FS_AUX_FROM_D;
             }
             */
          if(dp[i][j_from].del_score - gap_sum_j < dp[i][j_from].match_score - gap_sum_j) {
              dp[i][j].del_score = dp[i][j_from].match_score - gap_sum_j;
              dp[i][j].del_from = FMAP_SAM2FS_AUX_FROM_M;
          }
          else {
              dp[i][j].del_score = dp[i][j_from].del_score - gap_sum_j;
              dp[i][j].del_from = FMAP_SAM2FS_AUX_FROM_D;
          } 

          // vertical
          /*
             if(dp[i_from][j].ins_score - FMAP_SAM2FS_AUX_GAP_E < dp[i_from][j].match_score - FMAP_SAM2FS_AUX_GAP_O) {
             dp[i][j].ins_score = dp[i_from][j].match_score - FMAP_SAM2FS_AUX_GAP_O;
             dp[i][j].ins_from = FMAP_SAM2FS_AUX_FROM_M;
             }
             else {
             dp[i][j].ins_score = dp[i_from][j].ins_score - FMAP_SAM2FS_AUX_GAP_E;
             dp[i][j].ins_from = FMAP_SAM2FS_AUX_FROM_I;
             }
             */
          if(dp[i_from][j].ins_score - gap_sum_i < dp[i_from][j].match_score - gap_sum_i) {
              dp[i][j].ins_score = dp[i_from][j].match_score - gap_sum_i;
              dp[i][j].ins_from = FMAP_SAM2FS_AUX_FROM_M;
          }
          else {
              dp[i][j].ins_score = dp[i_from][j].ins_score - gap_sum_i;
              dp[i][j].ins_from = FMAP_SAM2FS_AUX_FROM_I;
          }

          // diagonal
          int32_t s = ((f_qseq->flow[i-1] < f_tseq->flow[j-1]) ? (f_tseq->flow[j-1]-f_qseq->flow[i-1]) : (f_qseq->flow[i-1]-f_tseq->flow[j-1]));
          if(dp[i-1][j-1].ins_score <= dp[i-1][j-1].match_score) {
              if(dp[i-1][j-1].del_score <= dp[i-1][j-1].match_score) {
                  dp[i][j].match_score = dp[i-1][j-1].match_score - s;
                  dp[i][j].match_from = FMAP_SAM2FS_AUX_FROM_M;
              }
              else {
                  dp[i][j].match_score = dp[i-1][j-1].del_score - s;
                  dp[i][j].match_from = FMAP_SAM2FS_AUX_FROM_D;
              }
          }
          else {
              if(dp[i-1][j-1].del_score <= dp[i-1][j-1].ins_score) {
                  dp[i][j].match_score = dp[i-1][j-1].ins_score - s;
                  dp[i][j].match_from = FMAP_SAM2FS_AUX_FROM_I;
              }
              else {
                  dp[i][j].match_score = dp[i-1][j-1].del_score - s;
                  dp[i][j].match_from = FMAP_SAM2FS_AUX_FROM_D;
              }
          }
          /*
             fprintf(stderr, "i=%d j=%d M[%d,%d] I[%d,%d] D[%d,%d]\n",
             i, j,
             dp[i][j].match_score, dp[i][j].match_from,
             dp[i][j].ins_score, dp[i][j].ins_from,
             dp[i][j].del_score, dp[i][j].del_from);
             */
      }
  }

  // Get best scoring cell
  best_score = FMAP_SW_MINOR_INF; 
  best_ctype = FMAP_SAM2FS_AUX_FROM_S;
  best_i = -1;
  best_j = -1;

  /*
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
  */
  if(best_score < dp[f_qseq->l][f_tseq->l].del_score) {
      best_i = f_qseq->l;
      best_j = f_tseq->l;
      best_score = dp[f_qseq->l][f_tseq->l].del_score;
      best_ctype = FMAP_SAM2FS_AUX_FROM_D;
  }
  if(best_score < dp[f_qseq->l][f_tseq->l].ins_score) {
      best_i = f_qseq->l;
      best_j = f_tseq->l;
      best_score = dp[f_qseq->l][f_tseq->l].ins_score;
      best_ctype = FMAP_SAM2FS_AUX_FROM_I;
  }
  if(best_score < dp[f_qseq->l][f_tseq->l].match_score) {
      best_i = f_qseq->l;
      best_j = f_tseq->l;
      best_score = dp[f_qseq->l][j].match_score;
      best_ctype = FMAP_SAM2FS_AUX_FROM_M;
  }

  i = best_i;
  j = best_j;
  ctype = best_ctype;

  aln = fmap_sam2fs_aux_aln_init();

  /*
     fprintf(stderr, "f_qseq->l=%d f_tseq->l=%d\n",
     f_qseq->l, f_tseq->l);
     */

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
              if(i <= 0) break;
              fmap_sam2fs_aux_aln_add(aln, f_qseq->flow[i-1], -1);
              i--;
          }
          break;
        case FMAP_SAM2FS_AUX_FROM_D:
          next_ctype = dp[i][j].del_from;
          for(k=0;k<4;k++) {
              if(j <= 0) break;
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
