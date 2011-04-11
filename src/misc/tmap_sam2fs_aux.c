/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "../util/tmap_definitions.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_error.h"
#include "../sw/tmap_sw.h"
#include "tmap_sam2fs_aux.h"

// do not use
#define TMAP_SAM2FS_AUX_GAP_O 4
#define TMAP_SAM2FS_AUX_GAP_E 4

// added gap penalty
#define TMAP_SAM2FS_AUX_GAP 0

typedef struct {
    int32_t *qseq;
    int32_t *tseq;
    char *aln;
    int32_t len, mem;
} tmap_sam2fs_aux_aln_t;

static tmap_sam2fs_aux_aln_t * 
tmap_sam2fs_aux_aln_init()
{
  return tmap_calloc(1, sizeof(tmap_sam2fs_aux_aln_t), "return");
}

static void
tmap_sam2fs_aux_aln_destroy(tmap_sam2fs_aux_aln_t *aln)
{
  free(aln->qseq);
  free(aln->tseq);
  free(aln->aln);
  free(aln);
}

static void
tmap_sam2fs_aux_aln_add(tmap_sam2fs_aux_aln_t *aln, int32_t read_n, int32_t ref_n)
{
  if(aln->mem <= aln->len) {
      aln->mem = aln->len + 1;
      tmap_roundup32(aln->mem);
      aln->qseq = tmap_realloc(aln->qseq, sizeof(int32_t) * aln->mem, "aln->qseq");
      aln->tseq = tmap_realloc(aln->tseq, sizeof(int32_t) * aln->mem, "aln->tseq");
      aln->aln = tmap_realloc(aln->aln, sizeof(char) * aln->mem, "aln->aln");
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
    TMAP_SAM2FS_AUX_FROM_M = 0,
    TMAP_SAM2FS_AUX_FROM_I = 1,
    TMAP_SAM2FS_AUX_FROM_D = 2,
    TMAP_SAM2FS_AUX_FROM_S = 3
};

typedef struct {
    int32_t match_score;
    int32_t match_from;
    int32_t ins_score;
    int32_t ins_from;
    int32_t del_score;
    int32_t del_from;
} tmap_sam2fs_aux_dpcell_t;

typedef struct {
    uint8_t *flow;
    int32_t l, m;
} tmap_sam2fs_aux_flow_t;

static tmap_sam2fs_aux_flow_t *
tmap_sam2fs_aux_flow_init()
{
  return tmap_calloc(1, sizeof(tmap_sam2fs_aux_flow_t), "return");
}

static void
tmap_sam2fs_aux_flow_convert(tmap_sam2fs_aux_flow_t *a, uint8_t *seq, int32_t len, 
                             uint8_t *flow_order, int32_t flow_order_len, uint8_t *qseq, int32_t qseq_len)
{
  int32_t i, k, l, next_i;

  i = k = a->l = a->m = 0;

  // move to the first flow
  if(NULL != qseq) {
      // we skip the flows not found in the qseq
      while('-' == qseq[i]) { // this could move beyond the end if all gaps
          if(qseq_len <= i) {
              tmap_error(NULL, Exit, OutOfRange);
          }
          i++;
      }
      while(flow_order[k] != qseq[i]) {
          k = (k+1) % flow_order_len;
      }
  }
  else {
      while('-' == seq[i]) { // this could move beyond the end if all gaps
          if(len <= i) {
              tmap_error(NULL, Exit, OutOfRange);
          }
          i++;
      }
      while(flow_order[k] != seq[i]) {
          k = (k+1) % flow_order_len;
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

      /*
      if(3 < seq[i]) {
          fprintf(stderr, "i=%d seq[i]=%d len=%d NULL=%d\n", i, seq[i], len, (NULL==qseq) ? 1 : 0);
          tmap_error(NULL, Exit, OutOfRange);
      }
      */

      // skip over empty flow
      while(flow_order[k] != seq[i] && seq[i] <= 3) {
          if(a->m <= a->l) {
              a->m = a->l + 1;
              tmap_roundup32(a->m);
              a->flow = tmap_realloc(a->flow, a->m * sizeof(uint8_t), "a->flow");
          }
          a->flow[a->l] = 0;
          a->l++;
          k = (k+1) % flow_order_len;
      }

      // get the number of bases in current flow
      next_i = i+1;
      l = 1;
      while(next_i < len 
            && (flow_order[k] == seq[next_i] || '-' == seq[next_i] || 3 < seq[next_i])) {
          if(flow_order[k] == seq[next_i] || 3 < seq[i]) {
              l++;
          }
          next_i++;
      }
      if(a->m <= a->l) { // more memory?
          a->m = a->l + 1;
          tmap_roundup32(a->m);
          a->flow = tmap_realloc(a->flow, a->m * sizeof(uint8_t), "a->flow");
      }
      a->flow[a->l] = l;
      a->l++;
      k = (k+1) % flow_order_len;
      i = next_i;
      if(i <= before) {
          tmap_error(NULL, Exit, OutOfRange);
      }
      // move beyond the gaps
      while(i < len && '-' == seq[i]) {
          i++;
      }
  }

  if(NULL != qseq) {
      k = (k+flow_order_len-1) % flow_order_len;

      // get the last base in qseq
      i = qseq_len-1;
      while('-' == qseq[i]) {
          i--;
          if(i < 0) {
              tmap_error(NULL, Exit, OutOfRange);
          }
      }

      // add empty flows
      while(flow_order[k] != qseq[i]) {
          if(a->m <= a->l) {
              a->m = a->l + 1;
              tmap_roundup32(a->m);
              a->flow = tmap_realloc(a->flow, a->m * sizeof(uint8_t), "a->flow");
          }
          a->flow[a->l] = 0;
          a->l++;
          k = (k+1) % flow_order_len;
      }
  }
}

static void
tmap_sam2fs_aux_flow_destroy(tmap_sam2fs_aux_flow_t *a)
{
  free(a->flow);
  free(a);
}

tmap_sam2fs_aux_flow_order_t *
tmap_sam2fs_aux_flow_order_init(char *flow_order)
{
  int32_t i, j, k;
  tmap_sam2fs_aux_flow_order_t *f = NULL;

  f = tmap_calloc(1, sizeof(tmap_sam2fs_aux_flow_order_t), "f");

  f->flow_order_len = strlen(flow_order);
  f->flow_order = tmap_malloc(sizeof(uint8_t) * f->flow_order_len, "f->flow_order");
  f->jump_fwd = tmap_calloc(f->flow_order_len, sizeof(uint16_t), "f->jump_fwd");
  f->jump_rev = tmap_calloc(f->flow_order_len, sizeof(uint16_t), "f->jump_rev");

  for(i=0;i<f->flow_order_len;i++) {
      f->flow_order[i] = tmap_nt_char_to_int[(int)flow_order[i]];
  }

  for(i=0;i<f->flow_order_len;i++) {
      // find the next position in the flow order
      k = 1;
      j = (i + 1) % f->flow_order_len;
      while(f->flow_order[i] != f->flow_order[j] && i != j) {
          j = (j + 1) % f->flow_order_len;
          k++;
      }
      f->jump_fwd[i] = f->jump_rev[j] = k;
  }

  return f;
}

void
tmap_sam2fs_aux_flow_order_destroy(tmap_sam2fs_aux_flow_order_t *f)
{
  free(f->flow_order);
  free(f->jump_fwd);
  free(f->jump_rev);
  free(f);
}

// query - read
// target - reference
void
tmap_sam2fs_aux_flow_align(tmap_file_t *fp, uint8_t *qseq, int32_t qseq_len, 
                           uint8_t *tseq, int32_t tseq_len, 
                           tmap_sam2fs_aux_flow_order_t *flow_order,
                           int8_t strand, char sep)
{
  int32_t i, j, k;
  int32_t ctype;
  int32_t best_i, best_j, best_score, best_ctype;
  tmap_sam2fs_aux_flow_t *f_qseq=NULL, *f_tseq=NULL;
  tmap_sam2fs_aux_dpcell_t **dp=NULL;
  tmap_sam2fs_aux_aln_t *aln=NULL;
  int32_t *gap_sums_i=NULL, *gap_sums_j=NULL;

  // init
  f_qseq = tmap_sam2fs_aux_flow_init();
  f_tseq = tmap_sam2fs_aux_flow_init();

  if(1 == strand) { // reverse
      for(i=0;i<(qseq_len>>1);i++) {
          uint8_t tmp = qseq[i];
          qseq[i] = qseq[qseq_len-i-1];
          qseq[qseq_len-i-1] = tmp;
      }
      for(i=0;i<(tseq_len>>1);i++) {
          uint8_t tmp = tseq[i];
          tseq[i] = tseq[tseq_len-i-1];
          tseq[tseq_len-i-1] = tmp;
      }
      for(i=0;i<flow_order->flow_order_len>>1;i++) {
          uint8_t tmp = flow_order->flow_order[i];
          flow_order->flow_order[i] = flow_order->flow_order[flow_order->flow_order_len-i-1];
          flow_order->flow_order[flow_order->flow_order_len-i-1] = tmp;
      }
  }

  // convert bases to flow space
  // NOTE: this will pad empty flows at the beginning and the end of the
  // reference
  tmap_sam2fs_aux_flow_convert(f_qseq, qseq, qseq_len, flow_order->flow_order, flow_order->flow_order_len, NULL, 0);
  tmap_sam2fs_aux_flow_convert(f_tseq, tseq, tseq_len, flow_order->flow_order, flow_order->flow_order_len, qseq, qseq_len);

  /*
  // DEBUGGING
  fputc('\n', stderr);
  // flow order
  for(i=0;i<flow_order->flow_order_len;i++) {
      fputc("ACGT"[flow_order->flow_order[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<flow_order->flow_order_len;i++) {
      if(0 < i) fputc(',', stderr);
      fprintf(stderr, "%d", flow_order->jump_fwd[i]);
  }
  fputc('\n', stderr);
  for(i=0;i<flow_order->flow_order_len;i++) {
      if(0 < i) fputc(',', stderr);
      fprintf(stderr, "%d", flow_order->jump_rev[i]);
  }
  fputc('\n', stderr);
  // read - query
  for(i=0;i<qseq_len;i++) {
      fputc("ACGT"[qseq[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<f_qseq->l;i++) {
      if(0 < i) fputc(',', stderr);
      fprintf(stderr, "%d", f_qseq->flow[i]);
  }
  fputc('\n', stderr);
  // ref - target
  for(i=0;i<tseq_len;i++) {
      fputc("ACGT"[tseq[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<f_tseq->l;i++) {
      if(0 < i) fputc(',', stderr);
      fprintf(stderr, "%d", f_tseq->flow[i]);
  }
  fputc('\n', stderr);
  fputc('\n', stderr);
  */

  // init gap sums
  gap_sums_i = tmap_malloc(sizeof(int32_t)*f_qseq->l, "gap_sums_i");
  for(i=0;i<f_qseq->l;i++) {
      k = i % flow_order->flow_order_len;
      j = (i < flow_order->jump_rev[k]) ? 0 : (i - flow_order->jump_rev[k]);
      gap_sums_i[i] = 0;
      while(j <= i) {
          gap_sums_i[i] += f_qseq->flow[i];
          j++;
      }
  }
  gap_sums_j = tmap_malloc(sizeof(int32_t)*f_tseq->l, "gap_sums_j");
  for(i=0;i<f_tseq->l;i++) {
      k = i % flow_order->flow_order_len;
      j = (i < flow_order->jump_rev[k]) ? 0 : (i - flow_order->jump_rev[k]);
      gap_sums_j[i] = 0;
      while(j <= i) {
          gap_sums_j[i] += f_tseq->flow[i];
          j++;
      }
  }

  /*
  // DEBUGGING
  fputc('\n', stderr);
  for(i=0;i<f_qseq->l;i++) {
      if(0 < i) fputc(',', stderr);
      fprintf(stderr, "%d", gap_sums_i[i]);
  }
  fputc('\n', stderr);
  for(i=0;i<f_tseq->l;i++) {
      if(0 < i) fputc(',', stderr);
      fprintf(stderr, "%d", gap_sums_j[i]);
  }
  fputc('\n', stderr);
  fputc('\n', stderr);
  */

  // init dp matrix
  dp = tmap_malloc(sizeof(tmap_sam2fs_aux_dpcell_t*)*(1+f_qseq->l), "dp");
  for(i=0;i<=f_qseq->l;i++) {
      dp[i] = tmap_malloc(sizeof(tmap_sam2fs_aux_dpcell_t)*(1+f_tseq->l), "dp[i]");
  }
  for(i=0;i<=f_qseq->l;i++) {
      for(j=0;j<=f_tseq->l;j++) {
          dp[i][j].match_score = dp[i][j].ins_score = dp[i][j].del_score = TMAP_SW_MINOR_INF;
          dp[i][j].match_from = dp[i][j].ins_from = dp[i][j].del_from = TMAP_SAM2FS_AUX_FROM_S;
      }
  }

  // init start cells
  dp[0][0].match_score = 0; // global
  for(i=k=1;i<=f_qseq->l;i++) {
      dp[i][0].ins_score = 0 - gap_sums_i[i-1];
      //dp[i][0].ins_score = 0 - TMAP_SAM2FS_AUX_GAP_O - (k-1)*TMAP_SAM2FS_AUX_GAP_E;
      dp[i][0].ins_from = ((k==1) ? TMAP_SAM2FS_AUX_FROM_M : TMAP_SAM2FS_AUX_FROM_I);
      if(0 == (i%flow_order->flow_order_len)) {
          k++;
      }
  }
  for(j=k=1;j<=f_tseq->l;j++) {
      dp[0][j].del_score = 0 - gap_sums_j[j-1];
      //dp[0][j].del_score = 0 - TMAP_SAM2FS_AUX_GAP_O - (k-1)*TMAP_SAM2FS_AUX_GAP_E;
      dp[0][j].del_from = ((1 == k) ? TMAP_SAM2FS_AUX_FROM_M : TMAP_SAM2FS_AUX_FROM_D);
      if(0 == (j%flow_order->flow_order_len)) {
          k++;
      }
  }
  // align
  for(i=1;i<=f_qseq->l;i++) { // query
      //int32_t i_from = ((i<flow_order->flow_order_len) ? 0 : (i-flow_order->flow_order_len)); // don't go before the first row
      int32_t i_from = ((i<flow_order->jump_rev[i-1]) ? 0 : (i-flow_order->jump_rev[i-1])); // don't go before the first row

      for(j=1;j<=f_tseq->l;j++) { // target
          int32_t j_from = ((j<flow_order->jump_rev[j-1]) ? 0 : (j-flow_order->jump_rev[j-1])); // don't go before the first column

          // horizontal
          /*
             if(dp[i][j_from].del_score - TMAP_SAM2FS_AUX_GAP_E < dp[i][j_from].match_score - TMAP_SAM2FS_AUX_GAP_O) {
             dp[i][j].del_score = dp[i][j_from].match_score - TMAP_SAM2FS_AUX_GAP_O;
             dp[i][j].del_from = TMAP_SAM2FS_AUX_FROM_M;
             }
             else {
             dp[i][j].del_score = dp[i][j_from].del_score - TMAP_SAM2FS_AUX_GAP_E;
             dp[i][j].del_from = TMAP_SAM2FS_AUX_FROM_D;
             }
             */
          if(dp[i][j_from].del_score - gap_sums_j[j-1] < dp[i][j_from].match_score - gap_sums_j[j-1]) {
              dp[i][j].del_score = dp[i][j_from].match_score - gap_sums_j[j-1];
              dp[i][j].del_from = TMAP_SAM2FS_AUX_FROM_M;
          }
          else {
              dp[i][j].del_score = dp[i][j_from].del_score - gap_sums_j[j-1];
              dp[i][j].del_from = TMAP_SAM2FS_AUX_FROM_D;
          } 

          // vertical
          /*
             if(dp[i_from][j].ins_score - TMAP_SAM2FS_AUX_GAP_E < dp[i_from][j].match_score - TMAP_SAM2FS_AUX_GAP_O) {
             dp[i][j].ins_score = dp[i_from][j].match_score - TMAP_SAM2FS_AUX_GAP_O;
             dp[i][j].ins_from = TMAP_SAM2FS_AUX_FROM_M;
             }
             else {
             dp[i][j].ins_score = dp[i_from][j].ins_score - TMAP_SAM2FS_AUX_GAP_E;
             dp[i][j].ins_from = TMAP_SAM2FS_AUX_FROM_I;
             }
             */
          if(dp[i_from][j].ins_score - gap_sums_i[i-1] < dp[i_from][j].match_score - gap_sums_i[i-1]) {
              dp[i][j].ins_score = dp[i_from][j].match_score - gap_sums_i[i-1];
              dp[i][j].ins_from = TMAP_SAM2FS_AUX_FROM_M;
          }
          else {
              dp[i][j].ins_score = dp[i_from][j].ins_score - gap_sums_i[i-1];
              dp[i][j].ins_from = TMAP_SAM2FS_AUX_FROM_I;
          }

          // diagonal
          int32_t s = ((f_qseq->flow[i-1] < f_tseq->flow[j-1]) ? (f_tseq->flow[j-1]-f_qseq->flow[i-1]) : (f_qseq->flow[i-1]-f_tseq->flow[j-1]));
          if(dp[i-1][j-1].ins_score <= dp[i-1][j-1].match_score) {
              if(dp[i-1][j-1].del_score <= dp[i-1][j-1].match_score) {
                  dp[i][j].match_score = dp[i-1][j-1].match_score - s;
                  dp[i][j].match_from = TMAP_SAM2FS_AUX_FROM_M;
              }
              else {
                  dp[i][j].match_score = dp[i-1][j-1].del_score - s;
                  dp[i][j].match_from = TMAP_SAM2FS_AUX_FROM_D;
              }
          }
          else {
              if(dp[i-1][j-1].del_score <= dp[i-1][j-1].ins_score) {
                  dp[i][j].match_score = dp[i-1][j-1].ins_score - s;
                  dp[i][j].match_from = TMAP_SAM2FS_AUX_FROM_I;
              }
              else {
                  dp[i][j].match_score = dp[i-1][j-1].del_score - s;
                  dp[i][j].match_from = TMAP_SAM2FS_AUX_FROM_D;
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
  best_score = TMAP_SW_MINOR_INF; 
  best_ctype = TMAP_SAM2FS_AUX_FROM_S;
  best_i = -1;
  best_j = -1;

  /*
  // last four flows in the last row
  for(k=0;k<flow_order->flow_order_len;k++) {
  i = f_qseq->l - k;
  if(i <= 0) {
  break;
  }
  if(best_score < dp[i][f_tseq->l].del_score) {
  best_i = i;
  best_j = f_tseq->l;
  best_score = dp[i][f_tseq->l].del_score;
  best_ctype = TMAP_SAM2FS_AUX_FROM_D;
  }
  if(best_score < dp[i][f_tseq->l].ins_score) {
  best_i = i;
  best_j = f_tseq->l;
  best_score = dp[i][f_tseq->l].ins_score;
  best_ctype = TMAP_SAM2FS_AUX_FROM_I;
  }
  if(best_score < dp[i][f_tseq->l].match_score) {
  best_i = i;
  best_j = f_tseq->l;
  best_score = dp[i][f_tseq->l].match_score;
  best_ctype = TMAP_SAM2FS_AUX_FROM_M;
  }
  }
  // last col, last four flow
  for(k=0;k<flow_order->flow_order_len;k++) {
  j = f_tseq->l - k;
  if(j <= 0) {
  break;
  }
  if(best_score < dp[f_qseq->l][j].del_score) {
  best_i = f_qseq->l;
  best_j = j;
  best_score = dp[f_qseq->l][j].del_score;
  best_ctype = TMAP_SAM2FS_AUX_FROM_D;
  }
  if(best_score < dp[f_qseq->l][j].ins_score) {
  best_i = f_qseq->l;
  best_j = j;
  best_score = dp[f_qseq->l][j].ins_score;
  best_ctype = TMAP_SAM2FS_AUX_FROM_I;
  }
  if(best_score < dp[f_qseq->l][j].match_score) {
  best_i = f_qseq->l;
  best_j = j;
  best_score = dp[f_qseq->l][j].match_score;
  best_ctype = TMAP_SAM2FS_AUX_FROM_M;
  }
  }
  */
  /*
  fprintf(stderr, "best_score=%d MID{%d,%d,%d}\n",
          best_score,
          dp[f_qseq->l][f_tseq->l].match_score,
          dp[f_qseq->l][f_tseq->l].ins_score,
          dp[f_qseq->l][f_tseq->l].del_score);
          */
  if(best_score < dp[f_qseq->l][f_tseq->l].del_score) {
      best_i = f_qseq->l;
      best_j = f_tseq->l;
      best_score = dp[f_qseq->l][f_tseq->l].del_score;
      best_ctype = TMAP_SAM2FS_AUX_FROM_D;
  }
  if(best_score < dp[f_qseq->l][f_tseq->l].ins_score) {
      best_i = f_qseq->l;
      best_j = f_tseq->l;
      best_score = dp[f_qseq->l][f_tseq->l].ins_score;
      best_ctype = TMAP_SAM2FS_AUX_FROM_I;
  }
  if(best_score < dp[f_qseq->l][f_tseq->l].match_score) {
      best_i = f_qseq->l;
      best_j = f_tseq->l;
      best_score = dp[f_qseq->l][f_tseq->l].match_score;
      best_ctype = TMAP_SAM2FS_AUX_FROM_M;
  }

  /*
  fprintf(stderr, "best_i=%d best_j=%d best_score=%d best_ctype=%d\n",
          best_i, best_j, best_score, best_ctype);
          */

  i = best_i;
  j = best_j;
  ctype = best_ctype;

  aln = tmap_sam2fs_aux_aln_init();

  // add in alignment end padding
  for(k=f_qseq->l;i<k;k--) { // insertion
      tmap_sam2fs_aux_aln_add(aln, f_qseq->flow[k-1], -1);
  }
  for(k=f_tseq->l;j<k;k--) { // deletion
      tmap_sam2fs_aux_aln_add(aln, -1, f_tseq->flow[k-1]);
  }

  // trace path back
  while(0 < i || 0 < j) {
      int32_t next_ctype = -1;

      switch(ctype) {
        case TMAP_SAM2FS_AUX_FROM_M:
          next_ctype = dp[i][j].match_from;
          tmap_sam2fs_aux_aln_add(aln, f_qseq->flow[i-1], f_tseq->flow[j-1]);
          i--;
          j--;
          break;
        case TMAP_SAM2FS_AUX_FROM_I:
          next_ctype = dp[i][j].ins_from;
          for(k=0;k<flow_order->flow_order_len;k++) {
              if(i <= 0) break;
              tmap_sam2fs_aux_aln_add(aln, f_qseq->flow[i-1], -1);
              i--;
          }
          break;
        case TMAP_SAM2FS_AUX_FROM_D:
          next_ctype = dp[i][j].del_from;
          for(k=0;k<flow_order->flow_order_len;k++) {
              if(j <= 0) break;
              tmap_sam2fs_aux_aln_add(aln, -1, f_tseq->flow[j-1]);
              j--;
          }
          break;
        default:
          //fprintf(stderr, "i=%d j=%d ctype=%d\n", i, j, ctype);
          tmap_error(NULL, Exit, OutOfRange);
      }

      /*
      fprintf(stderr, "i=%d j=%d ctype=%d %c next=%d %c\n", 
              i, j, 
              ctype, "MIDS"[ctype], 
              next_ctype, "MIDS"[next_ctype]);
              */
      ctype = next_ctype;
  }

  if(1 == strand) { // reverse back
      for(i=0;i<(aln->len>>1);i++) {
          char c;
          // qseq
          c = aln->qseq[i];
          aln->qseq[i] = aln->qseq[aln->len-i-1];
          aln->qseq[aln->len-i-1] = c;
          // tseq
          c = aln->tseq[i];
          aln->tseq[i] = aln->tseq[aln->len-i-1];
          aln->tseq[aln->len-i-1] = c;
          // aln
          c = aln->aln[i];
          aln->aln[i] = aln->aln[aln->len-i-1];
          aln->aln[aln->len-i-1] = c;
      }
  }

  // print
  // read - query
  for(i=aln->len-1;0<=i;i--) {
      if(0 <= aln->qseq[i]) tmap_file_fprintf(fp, "%d", aln->qseq[i]);
      else tmap_file_fprintf(fp, "-");
      if(0 < i) tmap_file_fprintf(fp, ","); 
  }
  // match string
  tmap_file_fprintf(fp, "%c", sep); 
  for(i=aln->len-1;0<=i;i--) {
      tmap_file_fprintf(fp, "%c", aln->aln[i]);
      if(0 < i) tmap_file_fprintf(fp, ","); 
  }
  // ref - target
  tmap_file_fprintf(fp, "%c", sep); 
  for(i=aln->len-1;0<=i;i--) {
      if(0 <= aln->tseq[i]) tmap_file_fprintf(fp, "%d", aln->tseq[i]);
      else tmap_file_fprintf(fp, "-");
      if(0 < i) tmap_file_fprintf(fp, ","); 
  }
  
  if(1 == strand) { // reverse back
      for(i=0;i<(qseq_len>>1);i++) {
          uint8_t tmp = qseq[i];
          qseq[i] = qseq[qseq_len-i-1];
          qseq[qseq_len-i-1] = tmp;
      }
      for(i=0;i<(tseq_len>>1);i++) {
          uint8_t tmp = tseq[i];
          tseq[i] = tseq[tseq_len-i-1];
          tseq[tseq_len-i-1] = tmp;
      }
      for(i=0;i<2;i++) {
          uint8_t tmp = flow_order->flow_order[i];
          flow_order->flow_order[i] = flow_order->flow_order[flow_order->flow_order_len-i-1];
          flow_order->flow_order[flow_order->flow_order_len-i-1] = tmp;
      }
  }

  // destroy
  for(i=0;i<=f_qseq->l;i++) {
      free(dp[i]);
  }
  free(dp);
  tmap_sam2fs_aux_flow_destroy(f_qseq);
  tmap_sam2fs_aux_flow_destroy(f_tseq);
  tmap_sam2fs_aux_aln_destroy(aln);
  free(gap_sums_i);
  free(gap_sums_j);
}
