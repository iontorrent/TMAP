/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <config.h>

#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_definitions.h"
#include "../seq/tmap_sff.h"
#include "../seq/tmap_seq.h"
#include "../io/tmap_file.h"
#include "../index/tmap_refseq.h"
#include "tmap_sw.h"
#include "../map/tmap_map_util.h"
#include "tmap_fsw.h"

/* TODO/IDEAS
   - do we need a transition probability added for going from one flow to the next?  For example, if we have a few flows that are empty, there is a certain probability of observing this... or this considered in your flow probabilities.  Basically, longer base alignments should be better, but using more flows is good/bad?
   - TODO: optimize memory usage by pre-allocating the cells/scores
   */

#define TMAP_FSW_ADD_FSCORE(s, f) ((s).match_score -= f, (s).ins_score -= f, (s).del_score -= f)

int32_t tmap_fsw_sm_short[] = {
    TMAP_MAP_OPT_SCORE_MATCH*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100,
    TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_SCORE_MATCH*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, 
    TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_SCORE_MATCH*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, 
    TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_SCORE_MATCH*-100, TMAP_MAP_OPT_PEN_MM*-100, 
    TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_PEN_MM*-100, TMAP_MAP_OPT_SCORE_MATCH*-100, 
};

static inline void
tmap_fsw_set_match(tmap_fsw_dpcell_t **dpcell, tmap_fsw_dpscore_t **dpscore, 
                   int32_t i, int32_t j, 
                   int64_t score, int32_t sub_path) 
{ 
  if(dpscore[i-1][j-1].match_score >= dpscore[i-1][j-1].ins_score) {
      if(dpscore[i-1][j-1].match_score >= dpscore[i-1][j-1].del_score) {
          if(1 == sub_path) dpcell[i][j].match_from = TMAP_FSW_FROM_M | 4; 
          else dpcell[i][j].match_from = 4 + dpcell[i-1][j-1].match_from;
          dpcell[i][j].match_bc = 1 + dpcell[i-1][j-1].match_bc;
          dpscore[i][j].match_score = dpscore[i-1][j-1].match_score + score;
      } 
      else { 
          if(1 == sub_path) dpcell[i][j].match_from = TMAP_FSW_FROM_D | 4;
          else dpcell[i][j].match_from = 4 + dpcell[i-1][j-1].del_from; 
          dpcell[i][j].match_bc = 1 + dpcell[i-1][j-1].del_bc;
          dpscore[i][j].match_score = dpscore[i-1][j-1].del_score + score;
      }
  } else {
      if(dpscore[i-1][j-1].ins_score >= dpscore[i-1][j-1].del_score) {
          if(1 == sub_path) dpcell[i][j].match_from = TMAP_FSW_FROM_I | 4;
          else dpcell[i][j].match_from = 4 + dpcell[i-1][j-1].ins_from; 
          dpcell[i][j].match_bc = 1 + dpcell[i-1][j-1].ins_bc;
          dpscore[i][j].match_score = dpscore[i-1][j-1].ins_score + score;
      } 
      else {
          if(1 == sub_path) dpcell[i][j].match_from = TMAP_FSW_FROM_D | 4;
          else dpcell[i][j].match_from = 4 + dpcell[i-1][j-1].del_from; 
          dpcell[i][j].match_bc = 1 + dpcell[i-1][j-1].del_bc;
          dpscore[i][j].match_score = dpscore[i-1][j-1].del_score + score;
      }
  }
}

static inline void 
tmap_fsw_set_ins(tmap_fsw_dpcell_t **dpcell, tmap_fsw_dpscore_t **dpscore, 
                 int32_t i, int32_t j, 
                 int32_t gap_open, int32_t gap_ext, int32_t sub_path) 
{
  if(dpscore[i-1][j].match_score - gap_open > dpscore[i-1][j].ins_score) {
      if(1 == sub_path) dpcell[i][j].ins_from = TMAP_FSW_FROM_M;
      else dpcell[i][j].ins_from = dpcell[i-1][j].match_from;
      dpcell[i][j].ins_bc = 1 + dpcell[i-1][j].match_bc;
      dpscore[i][j].ins_score = dpscore[i-1][j].match_score - gap_open - gap_ext;
  } else {
      if(1 == sub_path) dpcell[i][j].ins_from = TMAP_FSW_FROM_I;
      else dpcell[i][j].ins_from = dpcell[i-1][j].ins_from;
      dpcell[i][j].ins_bc = 1 + dpcell[i-1][j].ins_bc;
      dpscore[i][j].ins_score = dpscore[i-1][j].ins_score - gap_ext;
  }
}

static inline void 
tmap_fsw_set_end_ins(tmap_fsw_dpcell_t **dpcell, tmap_fsw_dpscore_t **dpscore, 
                     int32_t i, int32_t j, 
                     int32_t gap_open, int32_t gap_ext, int32_t gap_end, int32_t sub_path) 
{
  if(gap_end >= 0) {
      tmap_fsw_set_ins(dpcell, dpscore, i, j, gap_open, gap_end, sub_path);
  }
  else {
      tmap_fsw_set_ins(dpcell, dpscore, i, j, gap_open, gap_ext, sub_path);
  }
}

static inline void 
tmap_fsw_set_del(tmap_fsw_dpcell_t **dpcell, tmap_fsw_dpscore_t **dpscore, 
                 int32_t i, int32_t j, 
                 int32_t gap_open, int32_t gap_ext, int32_t sub_path) 
{
  if(dpscore[i][j-1].match_score - gap_open > dpscore[i][j-1].del_score) {
      if(1 == sub_path) dpcell[i][j].del_from = TMAP_FSW_FROM_M | 4;
      else dpcell[i][j].del_from = 4 + dpcell[i][j-1].match_from;
      dpcell[i][j].del_bc = dpcell[i][j-1].match_bc;
      dpscore[i][j].del_score = dpscore[i][j-1].match_score - gap_open - gap_ext;
  } else {
      if(1 == sub_path) dpcell[i][j].del_from = TMAP_FSW_FROM_D | 4;
      else dpcell[i][j].del_from = 4 + dpcell[i][j-1].del_from;
      dpcell[i][j].del_bc = dpcell[i][j-1].del_bc;
      dpscore[i][j].del_score = dpscore[i][j-1].del_score - gap_ext;
  }
}

static inline void 
tmap_fsw_set_end_del(tmap_fsw_dpcell_t **dpcell, tmap_fsw_dpscore_t **dpscore, 
                     int32_t i, int32_t j, 
                     int32_t gap_open, int32_t gap_ext, int32_t gap_end, int32_t sub_path) 
{
  if(gap_end >= 0) {
      tmap_fsw_set_del(dpcell, dpscore, i, j, gap_open, gap_end, sub_path);
  }
  else {
      tmap_fsw_set_del(dpcell, dpscore, i, j, gap_open, gap_ext, sub_path);
  }
}

inline void
tmap_fsw_sub_core(uint8_t *seq, int32_t len,
                  uint8_t flow_base, uint8_t base_call, uint16_t flow_signal,
                  const tmap_fsw_param_t *ap,
                  tmap_fsw_dpcell_t **sub_dpcell,
                  tmap_fsw_dpscore_t **sub_dpscore, 
                  tmap_fsw_dpcell_t *dpcell_last,
                  tmap_fsw_dpscore_t *dpscore_last,
                  tmap_fsw_dpcell_t *dpcell_curr,
                  tmap_fsw_dpscore_t *dpscore_curr,
                  tmap_fsw_path_t *path, int32_t *path_len, int32_t best_ctype,
                  uint8_t key_bases,
                  int32_t flowseq_start_clip, int32_t flowseq_end_clip)
{
  register int32_t i, j;
  int32_t low_offset, high_offset, flow_score;
  int32_t gap_open, gap_ext, gap_end, bw;
  int32_t *mat, *score_matrix, N_MATRIX_ROW;
  uint8_t offset;
  int32_t num_bases; 
  int32_t sub_path;

  gap_open = ap->gap_open;
  gap_ext = ap->gap_ext;
  gap_end = ap->gap_end;
  bw = ap->band_width;
  score_matrix = ap->matrix;
  N_MATRIX_ROW = ap->row;
  offset = ap->offset;

  // get the homopolymer bounds
  low_offset = (base_call < offset) ? 0 : base_call - offset;
  high_offset = base_call + offset;

  //fprintf(stderr, "base_call=%d low_offset=%d high_offset=%d\n", base_call, low_offset, high_offset);

  // scoring matrix will always be the same
  mat = score_matrix + flow_base * N_MATRIX_ROW;

  sub_path = (NULL == path) ? 0 : 1; // only if we wish to recover the path  
  /*
  fprintf(stderr, "flow_base=%c base_call=%d flow_signal=%d flowseq_start_clip=%d flowseq_end_clip=%d len=%d\n",
          "ACGTN"[flow_base], base_call, flow_signal,
          flowseq_start_clip, flowseq_end_clip, len);
          */

  // copy previous row
  for(j=0;j<=len;j++) {
      sub_dpscore[0][j] = dpscore_last[j];
      TMAP_FSW_INIT_CELL(sub_dpcell[0][j]); // empty
      sub_dpcell[0][j].match_from = TMAP_FSW_FROM_M;
      sub_dpcell[0][j].ins_from = TMAP_FSW_FROM_I;
      sub_dpcell[0][j].del_from = TMAP_FSW_FROM_D;
      sub_dpcell[0][j].match_bc = sub_dpcell[0][j].ins_bc = sub_dpcell[0][j].del_bc = 0;
  }

  // fill in sub_dpcell and sub_dpscore
  for(i=1;i<=high_offset;i++) { // for each row in the sub-alignment
      // initialize the first column
      TMAP_FSW_SET_SCORE_INF(sub_dpscore[i][0]); 
      TMAP_FSW_INIT_CELL(sub_dpcell[i][0]);
      tmap_fsw_set_end_ins(sub_dpcell, sub_dpscore, i, 0, gap_open, gap_ext, gap_end, sub_path);
      // fill in the rest of the columns
      for(j=1;j<=len;j++) { // for each col
          //if(1 == i) fprintf(stderr, "j=%d flow_base=%d seq[j-1]=%d mat[seq[j-1]]=%d\n", j, flow_base, seq[j-1], mat[seq[j-1]]); 
          tmap_fsw_set_match(sub_dpcell, sub_dpscore, i, j, mat[seq[j-1]], sub_path);
          tmap_fsw_set_ins(sub_dpcell, sub_dpscore, i, j, gap_open, gap_ext, sub_path);
          tmap_fsw_set_del(sub_dpcell, sub_dpscore, i, j, gap_open, gap_ext, sub_path);
          /*
          if(1 == flowseq_start_clip && sub_dpscore[i][j].match_score < 0) {
              sub_dpscore[i][j].match_score = 0;
              sub_dpcell[i][j].match_from = TMAP_FSW_FROM_S;
          }
          */
      }
  }
  
  // add flow scores
  for(i=low_offset;i<=high_offset;i++) { // for each possible base call within +-offset
      // get flow score for "(i-low_offset)" bases
      num_bases = i + key_bases; // key bases will shift this
      flow_score = (flow_signal < 100*num_bases) ? (100*num_bases - flow_signal) : (flow_signal - 100*num_bases);
      flow_score *= ap->fscore / 100;
      // factor out extra match scores for over calls.  Is this correct?
      //if(base_call < i) flow_score += mat[flow_base] * (i - base_call);
      //fprintf(stderr, "flow_score=%d i=%d\n", flow_score, i);
      if(flow_score < 0) tmap_error("bug encountered", Exit, OutOfRange); // we will subtract it
      for(j=0;j<=len;j++) { // for each col
          TMAP_FSW_ADD_FSCORE(sub_dpscore[i][j], flow_score);
          /*
          if(TMAP_FSW_FROM_S == sub_dpcell[i][j].match_from) {
              TMAP_FSW_ADD_FSCORE(sub_dpscore[i][j], flow_score);
              sub_dpscore[i][j].match_score = 0; // reset
          }
          else {
              TMAP_FSW_ADD_FSCORE(sub_dpscore[i][j], flow_score);
          }
          */
      }
  }

  if(NULL != dpcell_curr && NULL != dpscore_curr) {
      // set the best cell to be [base_call][0,len]
      // NOTE: set this to the original base call to get consistency between
      // calling homopolymer over/under calls
      for(j=0;j<=len;j++) { // for each col
          dpcell_curr[j] = sub_dpcell[base_call][j];
          dpscore_curr[j] = sub_dpscore[base_call][j];
          if(1 == flowseq_start_clip) { // start anywhere
              if(dpscore_curr[j].match_score < 0) {
                  dpcell_curr[j].match_from = TMAP_FSW_FROM_S;
                  dpcell_curr[j].match_bc = 0;
                  dpscore_curr[j].match_score = TMAP_SW_MINOR_INF;
              }
              if(dpscore_curr[j].ins_score < 0) {
                  dpcell_curr[j].ins_from = TMAP_FSW_FROM_S;
                  dpcell_curr[j].ins_bc = 0;
                  dpscore_curr[j].ins_score = TMAP_SW_MINOR_INF;
              }
              if(dpscore_curr[j].del_score < 0) {
                  dpcell_curr[j].del_from = TMAP_FSW_FROM_S;
                  dpcell_curr[j].del_bc = 0;
                  dpscore_curr[j].del_score = TMAP_SW_MINOR_INF;
              }
          }
      }

      // get the best cells within [low_offset+1,high_offset][0,len]
      for(i=low_offset;i<=high_offset;i++) {
          if(base_call != i) { 
              // break ties by preferring hp errors, hence "<=" below
              for(j=0;j<=len;j++) { // for each col
                  // match
                  if(dpscore_curr[j].match_score <= sub_dpscore[i][j].match_score) {
                      dpcell_curr[j].match_from = sub_dpcell[i][j].match_from;
                      dpcell_curr[j].match_bc = sub_dpcell[i][j].match_bc;
                      dpscore_curr[j].match_score = sub_dpscore[i][j].match_score;
                  }
                  // ins
                  if(dpscore_curr[j].ins_score <= sub_dpscore[i][j].ins_score) {
                      dpcell_curr[j].ins_from = sub_dpcell[i][j].ins_from;
                      dpcell_curr[j].ins_bc = sub_dpcell[i][j].ins_bc;
                      dpscore_curr[j].ins_score = sub_dpscore[i][j].ins_score;
                  }
                  // del
                  if(dpscore_curr[j].del_score <= sub_dpscore[i][j].del_score) {
                      dpcell_curr[j].del_from = sub_dpcell[i][j].del_from;
                      dpcell_curr[j].del_bc = sub_dpcell[i][j].del_bc;
                      dpscore_curr[j].del_score = sub_dpscore[i][j].del_score;
                  }
              }
          }
      }
  }

  if(NULL != path) {
      int32_t ctype, ctype_next = 0;
      tmap_fsw_path_t *p;

      p = path;
      i = high_offset; j = len; 
      ctype = best_ctype;

      //fprintf(stderr, "  START sub_core i=%d j=%d ctype=%d\n", i, j, ctype);
      while(TMAP_FSW_FROM_S != ctype && 0 < i) {
          //fprintf(stderr, "  sub_core i=%d j=%d ctype=%d\n", i, j, ctype);
          if(i < 0 || j < 0) tmap_error("bug encountered", Exit, OutOfRange);
          // is there enough path
          if(*path_len < p - path) tmap_error("bug encountered", Exit, OutOfRange);

          switch(ctype) { 
            case TMAP_FSW_FROM_M: 
              ctype_next = sub_dpcell[i][j].match_from & 0x3;
              break;
            case TMAP_FSW_FROM_I: 
              ctype_next = sub_dpcell[i][j].ins_from & 0x3;
              break;
            case TMAP_FSW_FROM_D: 
              ctype_next = sub_dpcell[i][j].del_from & 0x3;
              break;
            default:
              tmap_error(NULL, Exit, OutOfRange);
          }
          p->ctype = ctype;
          p->i = i-1;
          p->j = j-1;
          p++;


          // move the row and column (as necessary)
          switch(ctype) {
            case TMAP_FSW_FROM_M: 
              --i; --j; break;
            case TMAP_FSW_FROM_I: 
              --i; break;
            case TMAP_FSW_FROM_D: 
              --j; break;
            default:
              tmap_error(NULL, Exit, OutOfRange);
          }

          // move to the next cell type
          ctype = ctype_next;
          //fprintf(stderr, "  sub_core i=%d j=%d ctype_next=%d\n", i, j, ctype_next);
      }

      (*path_len) = p - path;
  }
  else if(NULL != path_len) {
      (*path_len) = 0;
  }
}

static void
tmap_fsw_get_path(uint8_t *seq, uint8_t *flow_order, int32_t flow_order_len, uint8_t *base_calls, uint16_t *flowgram,
                  int32_t key_index, int32_t key_bases,
                  tmap_fsw_dpcell_t **dpcell, tmap_fsw_dpscore_t **dpscore,
                  tmap_fsw_dpcell_t **sub_dpcell,
                  tmap_fsw_dpscore_t **sub_dpscore, 
                  const tmap_fsw_param_t *ap,
                  int32_t best_i, int32_t best_j, uint8_t best_ctype, 
                  int32_t flowseq_start_clip, int32_t flowseq_end_clip,
                  int32_t right_j,
                  tmap_fsw_path_t *path, int32_t *path_len)
{
  register int32_t i, j;
  int32_t k, l;
  int32_t base_call = 0, col_offset = 0, base_call_diff = 0;
  uint8_t ctype, ctype_next = 0;
  tmap_fsw_path_t *p;
  tmap_fsw_path_t *sub_path = NULL;
  int32_t sub_path_len = 0, sub_path_mem = 0;

  /*
     fprintf(stderr, "PRINTING\n");
     for(i=0;i<=best_i;i++) {
     for(j=0;j<=best_j;j++) {
     fprintf(stderr, "(%d,%d M[%d,%d,%d,%d] I[%d,%d,%d,%d] D[%d,%d,%d,%d])\n",
     i, j,
     dpcell[i][j].match_bc, dpcell[i][j].match_from & 0x3, dpcell[i][j].match_from >> 2, (int)dpscore[i][j].match_score,
     dpcell[i][j].ins_bc, dpcell[i][j].ins_from & 0x3, dpcell[i][j].ins_from >> 2, (int)dpscore[i][j].ins_score,
     dpcell[i][j].del_bc, dpcell[i][j].del_from & 0x3, dpcell[i][j].del_from >> 2, (int)dpscore[i][j].del_score);
     }
     }
     */

  // get best scoring end cell
  i = best_i; j = best_j; p = path;
  ctype = best_ctype;

  /*
  fprintf(stderr, "GETTING PATH\n");
  fprintf(stderr, "best_i=%d best_j=%d ctype=%d\n",
          best_i, best_j, ctype);
          */

  while(TMAP_FSW_FROM_S != ctype && (0 < i || 0 < j)) {
      base_call = 0;
      col_offset = 0;
      base_call_diff = 0;
      //fprintf(stderr, "CORE i=%d j=%d ctype=%d\n", i, j, ctype);

      // local
      //if(i <= 0 && 1 == flowseq_start_clip) {
      if(i <= 0) { // do not add in leading deletion...
          break;
      }

      switch(ctype) { 
        case TMAP_FSW_FROM_M: 
          if(i < 0 || j < 0) tmap_error("bug encountered", Exit, OutOfRange);
          base_call = dpcell[i][j].match_bc;
          col_offset = dpcell[i][j].match_from >> 2;
          ctype_next = dpcell[i][j].match_from & 0x3;
          break;
        case TMAP_FSW_FROM_I: 
          if(i < 0 || j < 0) tmap_error("bug encountered", Exit, OutOfRange);
          base_call = dpcell[i][j].ins_bc;
          col_offset = dpcell[i][j].ins_from >> 2;
          ctype_next = dpcell[i][j].ins_from & 0x3;
          break;
        case TMAP_FSW_FROM_D: 
          if(i < 0 || j < 0) tmap_error("bug encountered", Exit, OutOfRange);
          base_call = dpcell[i][j].del_bc;
          col_offset = dpcell[i][j].del_from >> 2;
          ctype_next = dpcell[i][j].del_from & 0x3;
          break;
        default:
          tmap_error(NULL, Exit, OutOfRange);
      }

      if(0 < i && 0 < j && ctype_next == TMAP_FSW_FROM_S) {
          break;
      }

      /*
      fprintf(stderr, "i=%d j=%d base_call=%d col_offset=%d ctype=%d ctype_next=%d base_calls[%d]=%d\n",
              i, j, base_call, col_offset, ctype, ctype_next,
              i-1,
              (0 < i) ? base_calls[i-1] : 0);
              */

      if(j <= 0) {
          // if base_call_diff > 0, add insertions (more read bases than reference bases)
          // if base_call_diff < 0, add deletions (fewer read bases than reference bases)
          base_call_diff = base_call - base_calls[i-1];
          // starts with an insertion?
          for(k=0;k<base_call;k++) {
              p->i = i-1;
              p->j = j-1;
              p->ctype = TMAP_FSW_FROM_I;
              p++;
              if(*path_len <= p - path) tmap_error("bug encountered", Exit, OutOfRange);
          }
          for(k=base_call_diff;k<0;k++) {
              p->i = i-1;
              p->j = j-1;
              p->ctype = TMAP_FSW_FROM_D;
              p++;
              if(*path_len <= p - path) tmap_error("bug encountered", Exit, OutOfRange);
          }
          //fprintf(stderr, "base_call=%d base_calls[%d]=%d base_call_diff=%d\n", base_call, i-1, base_calls[i-1], base_call_diff);
          i--;
      }
      else if(i <= 0) {
          // starts with a deletion?
          p->i = i-1;
          p->j = j-1;
          p->ctype = TMAP_FSW_FROM_D;
          p++;
          if(*path_len <= p - path) tmap_error("bug encountered", Exit, OutOfRange);
          j--;
      }
      else if(0 == base_call) {
          if(0 < i && 0 < base_calls[i-1]) {
              for(k=0;k<base_calls[i-1];k++) {
                  p->j = j - 1;
                  p->i = i - 1;
                  p->ctype = TMAP_FSW_FROM_HP_MINUS;
                  p++;
                  if(*path_len <= p - path) tmap_error("bug encountered", Exit, OutOfRange);
              }
          }
          i--;
          j -= col_offset;
      }
      else {
          while(sub_path_mem < j + ap->offset + 1) { // more memory please
              sub_path_mem = (0 == sub_path_mem) ? 4 : (sub_path_mem << 1);
              sub_path = tmap_realloc(sub_path, sizeof(tmap_fsw_path_t) * sub_path_mem, "sub_path"); 
          } 
          sub_path_len = sub_path_mem;

          tmap_fsw_param_t ap_tmp;
          ap_tmp = (*ap);
          ap_tmp.offset = 0; // no offset, since we will use the previous base call

          // solve the sub-problem and get the path
          if(j - col_offset < 0) tmap_error("bug encountered", Exit, OutOfRange);
          tmap_fsw_sub_core(seq, j,
                            flow_order[(i-1) % flow_order_len], base_call, flowgram[i-1], 
                            &ap_tmp,
                            sub_dpcell, sub_dpscore,
                            dpcell[i-1], dpscore[i-1],
                            NULL, NULL, // do not update
                            sub_path, &sub_path_len, ctype, // get the path
                            ((key_index+1) == i) ? key_bases : 0,
                            0, 0);

          // if base_call_diff > 0, add insertions (more read bases than reference bases)
          // if base_call_diff < 0, add deletions (fewer read bases than reference bases)
          base_call_diff = base_call - base_calls[i-1];
          /*
          fprintf(stderr, "base_call=%d base_calls[i-1]=%d sub_path_len=%d base_call_diff=%d\n", 
                  base_call, base_calls[i-1],
                  sub_path_len, base_call_diff);
          for(l=0;l<sub_path_len;l++) {
              fprintf(stderr, "l=%d sub_path[l].i=%d sub_path[l].j=%d sub_path[l].ctype=%d\n",
                      l, sub_path[l].i, sub_path[l].j, sub_path[l].ctype);
          }
          */

          // correctly justify indels
          if(1 == right_j) {
              if(0 < base_call_diff) { // change beginning matches to insertions
                  for(l=0;l<sub_path_len && 0 < base_call_diff;l++) {
                      if(TMAP_FSW_FROM_M == sub_path[l].ctype) { // delete a "match"
                          sub_path[l].ctype = TMAP_FSW_FROM_HP_PLUS;
                          base_call_diff--;
                      }
                  }
              }
              // add deletions before storing the path
              while(base_call_diff < 0) { // there are bases left that were deleted
                  p->ctype = TMAP_FSW_FROM_HP_MINUS; // add an insertion
                  p->i = i - 1;
                  p->j = j - 1; // is this correct ? 
                  p++;
                  if(*path_len <= p - path) tmap_error("bug encountered", Exit, OutOfRange);
                  base_call_diff++;
              }
              for(l=0;l<sub_path_len;l++) {
                  (*p) = sub_path[l]; // store p->j and p->ctype
                  p->i = i - 1; // same flow index
                  p++;
                  if(*path_len <= p - path) tmap_error("bug encountered", Exit, OutOfRange);
              }
          }
          else {
              if(0 < base_call_diff) { // change end matches to insertions
                  for(l=sub_path_len-1;0<=l && 0 < base_call_diff;l--) {
                      if(TMAP_FSW_FROM_M == sub_path[l].ctype) { // delete a "match"
                          sub_path[l].ctype = TMAP_FSW_FROM_HP_PLUS;
                          base_call_diff--;
                      }
                  }
              }
              k = j - 1; // for base_call_diff < 0
              for(l=0;l<sub_path_len;l++) {
                  (*p) = sub_path[l]; // store p->j and p->ctype
                  p->i = i - 1; // same flow index
                  k = p->j; // for base_call_diff < 0
                  p++;
                  if(*path_len <= p - path) tmap_error("bug encountered", Exit, OutOfRange);
              }
              // add deletions after storing the path
              while(base_call_diff < 0) { // there are bases left that were deleted
                  p->ctype = TMAP_FSW_FROM_HP_MINUS; // add an insertion
                  p->i = i - 1;
                  p->j = k; // is this correct ? 
                  p++;
                  if(*path_len <= p - path) tmap_error("bug encountered", Exit, OutOfRange);
                  base_call_diff++;
              }
          }
          if(0 != base_call_diff) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }

          // move the row and column (as necessary)
          i--;
          j -= col_offset;
      }

      // move to the next cell type
      ctype = ctype_next;

      //fprintf(stderr, "HERE i=%d j=%d ctype=%d\n", i, j, ctype);
  }
  (*path_len) = p - path;

  /*
  for(i=0;i<(*path_len);i++) {
      fprintf(stderr, "i=%d path[i].i=%d path[i].j=%d\n",
              i, path[i].i, path[i].j);
  }
  */

  free(sub_path);
}

/*
Notes: key_index is zero-base and should be -1, 0, or (num_flows-1)
*/
int64_t
tmap_fsw_clipping_core(uint8_t *seq, int32_t len, 
                    tmap_fsw_flowseq_t *flowseq,
                    const tmap_fsw_param_t *ap,
                    int32_t flowseq_start_clip, int32_t flowseq_end_clip,
                    int32_t right_j,
                    tmap_fsw_path_t *path, int32_t *path_len)
{
  register int32_t i, j;
  int32_t max_bc = 0, bw;

  // main cells 
  tmap_fsw_dpcell_t **dpcell;
  tmap_fsw_dpscore_t **dpscore;

  // for homopolymer re-calling 
  tmap_fsw_dpcell_t **sub_dpcell;
  tmap_fsw_dpscore_t **sub_dpscore;

  int32_t gap_open, gap_ext, gap_end;
  int32_t *score_matrix, N_MATRIX_ROW;
  uint8_t offset;

  int32_t best_i=-1, best_j=-1;
  uint8_t best_ctype=0;
  int64_t best_score = TMAP_SW_MINOR_INF;

  if(0 == flowseq->num_flows || 0 == len) {
      (*path_len) = 0;
      return 0;
  }
  
  /*
  fprintf(stderr, "%s flowseq_start_clip=%d flowseq_end_clip=%d\n",
          __func__, flowseq_start_clip, flowseq_end_clip);
          */

  gap_open = ap->gap_open;
  gap_ext = ap->gap_ext;
  gap_end = ap->gap_end;
  bw = ap->band_width;
  score_matrix = ap->matrix;
  N_MATRIX_ROW = ap->row;
  offset = ap->offset;

  // maximum length base call
  for(i=0;i<flowseq->num_flows;i++) {
      if(max_bc < flowseq->base_calls[i]) {
          max_bc = flowseq->base_calls[i];
      }
  }

  /*
  fprintf(stderr, "flowseq_start_clip=%d flowseq_end_clip=%d\n", flowseq_start_clip, flowseq_end_clip);
  fprintf(stderr, "flowseq->flow_order_len=%d\n", flowseq->flow_order_len);
  for(i=0;i<flowseq->flow_order_len;i++) {
      fputc("ACGTN"[flowseq->flow_order[i]], stderr);
  }
  fputc('\n', stderr);
  fprintf(stderr, "max_bc=%d offset=%d len=%d\n", max_bc, offset, len);
  fprintf(stderr, "flowseq->num_flows=%d len=%d\n", flowseq->num_flows, len);
  fprintf(stderr, "max i=%d max j=%d\n", max_bc+offset, len);
  */

  // allocate memory for the sub-cells
  sub_dpcell = tmap_malloc(sizeof(tmap_fsw_dpcell_t*) * (max_bc + offset + 1), "sub_dpcell");
  sub_dpscore = tmap_malloc(sizeof(tmap_fsw_dpscore_t*) * (max_bc + offset + 1), "sub_dpscore");
  for(i=0;i<=max_bc+offset;i++) {
      sub_dpcell[i] = tmap_malloc(sizeof(tmap_fsw_dpcell_t) * (len + 1), "sub_dpcell");
      sub_dpscore[i] = tmap_malloc(sizeof(tmap_fsw_dpscore_t) * (len + 1), "sub_dpscore");
  }

  // allocate memory for the main cells
  dpcell = tmap_malloc(sizeof(tmap_fsw_dpcell_t*) * (flowseq->num_flows + 1), "dpcell");
  dpscore = tmap_malloc(sizeof(tmap_fsw_dpscore_t*) * (flowseq->num_flows + 1), "dpscore");
  for(i=0;i<=flowseq->num_flows;i++) {
      dpcell[i] = tmap_malloc(sizeof(tmap_fsw_dpcell_t) * (len + 1), "dpcell");
      dpscore[i] = tmap_malloc(sizeof(tmap_fsw_dpscore_t) * (len + 1), "dpscore");
  }

  // set first row
  TMAP_FSW_SET_SCORE_INF(dpscore[0][0]); 
  TMAP_FSW_INIT_CELL(dpcell[0][0]);
  dpscore[0][0].match_score = 0;
  if(1 == flowseq_start_clip) { // start anywhere in flowseq
      for(j=1;j<=len;j++) { // for each col (flow in the reference)
          TMAP_FSW_SET_SCORE_INF(dpscore[0][j]);
          TMAP_FSW_INIT_CELL(dpcell[0][j]);
          // the alignment can start anywhere within seq and anywhere within
          // flowseq
          for(i=0;i<=flowseq->num_flows;i++) {
              dpscore[i][j].match_score = 0; 
              dpcell[i][j].match_from = TMAP_FSW_FROM_S; 
          }
      }
  }
  else { // start at the first flow in seq2
      for(j=1;j<=len;j++) { // for each col
          TMAP_FSW_SET_SCORE_INF(dpscore[0][j]);
          TMAP_FSW_INIT_CELL(dpcell[0][j]);
          tmap_fsw_set_end_del(dpcell, dpscore, 0, j, gap_open, gap_ext, gap_end, 0);
          // the alignment can start anywhere within seq 
          dpscore[0][j].match_score = 0; 
          dpcell[0][j].match_from = TMAP_FSW_FROM_S; 
      }
  }

  // core loop
  for(i=1;i<=flowseq->num_flows;i++) { // for each row (flow in the read)
      // fill in the columns
      /*
      fprintf(stderr, "i=%d num_flows=%d flow_base=%c base_calls=%d flowgram=%d\n", 
              i, flowseq->num_flows,
              "ACGTN"[flowseq->flow_order[(i-1) % flowseq->flow_order_len]], 
              flowseq->base_calls[i-1], 
              flowseq->flowgram[i-1]);
              */
      tmap_fsw_sub_core(seq, len,
                        flowseq->flow_order[(i-1) % flowseq->flow_order_len], 
                        flowseq->base_calls[i-1], 
                        flowseq->flowgram[i-1], 
                        ap,
                        sub_dpcell, sub_dpscore,
                        dpcell[i-1], dpscore[i-1],
                        dpcell[i], dpscore[i],
                        NULL, NULL, 0,
                        ((flowseq->key_index+1) == i) ? flowseq->key_bases : 0,
                        flowseq_start_clip, flowseq_end_clip);

      // deal with start clipping
      if(1 == flowseq_start_clip) {
          for(j=0;j<=len;j++) {
              if(dpscore[i][j].match_score < 0) {
                  //fprintf(stderr, "%s HERE 1 i=%d j=%d base_calls[i-1]=%d\n", __func__, i, j, flowseq->base_calls[i-1]);
                  dpcell[i][j].match_from = TMAP_FSW_FROM_S;
                  dpscore[i][j].match_score = 0; 
              }
          }
      }

      // Update best
      if(1 == flowseq_end_clip // end anywhere in flowseq
         || i == flowseq->num_flows) {
          for(j=1;j<=len;j++) {
              /*
              fprintf(stderr, "i=%d j=%d scores=[%d,%d,%d] from=[%d,%d,%d]\n",
                      i, j, 
                      dpscore[i][j].match_score,
                      dpscore[i][j].ins_score,
                      dpscore[i][j].del_score,
                      dpcell[i][j].match_from,
                      dpcell[i][j].ins_from,
                      dpcell[i][j].del_from);
                      */
              if(best_score < dpscore[i][j].match_score) {
                  best_score = dpscore[i][j].match_score;
                  best_ctype = TMAP_FSW_FROM_M;
                  best_i = i; best_j = j;
              }
              if(best_score < dpscore[i][j].ins_score) {
                  best_score = dpscore[i][j].ins_score;
                  best_ctype = TMAP_FSW_FROM_I;
                  best_i = i; best_j = j;
              }
              if(best_score < dpscore[i][j].del_score) {
                  best_score = dpscore[i][j].del_score;
                  best_ctype = TMAP_FSW_FROM_D;
                  best_i = i; best_j = j;
              }
          }
      }
  }
  //fprintf(stderr, "%s best_score=%d best_i=%d best_j=%d\n", __func__, (int)best_score, best_i, best_j);

  if(best_i < 0 || best_j < 0) { // was not updated
      (*path_len) = 0;
      // go to the end and free
  }
  else if(NULL == path) { 
      // do nothing
  }
  else if(NULL == path_len) {
      path[0].i = best_i; path[0].j = best_j;
  }
  else {
      // recover the path
      tmap_fsw_get_path(seq, flowseq->flow_order, flowseq->flow_order_len, flowseq->base_calls, flowseq->flowgram,
                        flowseq->key_index, flowseq->key_bases,
                        dpcell, dpscore, 
                        sub_dpcell, sub_dpscore, 
                        ap, 
                        best_i, best_j, best_ctype, 
                        flowseq_start_clip, flowseq_end_clip,
                        right_j,
                        path, path_len);
  }

  // free memory for the sub-cells
  for(i=0;i<=max_bc+offset;i++) {
      free(sub_dpcell[i]);
      free(sub_dpscore[i]);
  }
  free(sub_dpcell);
  free(sub_dpscore);

  // free memory for the main cells
  for(i=0;i<=flowseq->num_flows;i++) {
      free(dpcell[i]);
      free(dpscore[i]);
  }
  free(dpcell);
  free(dpscore);

  return best_score;
}

tmap_fsw_flowseq_t *
tmap_fsw_flowseq_init(uint8_t *flow_order, int32_t flow_order_len, uint8_t *base_calls, uint16_t *flowgram,
                                          int32_t num_flows, int32_t key_index, int32_t key_bases)
{
  tmap_fsw_flowseq_t *flowseq;

  flowseq = tmap_calloc(1, sizeof(tmap_fsw_flowseq_t), "flowseq");

  // shallow copy
  flowseq->flow_order = flow_order;
  flowseq->flow_order_len = flow_order_len;
  flowseq->base_calls = base_calls;
  flowseq->flowgram = flowgram;
  flowseq->num_flows = num_flows;
  flowseq->key_index = key_index;
  flowseq->key_bases = key_bases;
  flowseq->mem = num_flows;

  return flowseq;
}

void
tmap_fsw_flowseq_print(tmap_file_t *fp, tmap_fsw_flowseq_t *flowseq)
{
  int32_t i, j;

  for(i=0;i<flowseq->flow_order_len;i++) {
      tmap_file_fprintf(fp, "%c", "ACGTN"[flowseq->flow_order[i]]);
  }
  tmap_file_fprintf(fp, "\n");
  for(i=0;i<flowseq->num_flows;i++) {
      for(j=0;j<flowseq->base_calls[i];j++) {
          tmap_file_fprintf(fp, "%c", "ACGTN"[flowseq->flow_order[i % flowseq->flow_order_len]]);
      }
  }
  tmap_file_fprintf(fp, "\n");
  for(i=0;i<flowseq->num_flows;i++) {
      if(0 < i) tmap_file_fprintf(fp, ",");
      tmap_file_fprintf(fp, "%d", flowseq->flowgram[i]);
  }
  tmap_file_fprintf(fp, "\n");
  tmap_file_fprintf(fp, "key_index=%d\tkey_bases=%d\n", flowseq->key_index, flowseq->key_bases);
}

void
tmap_fsw_flowseq_destroy_shallow(tmap_fsw_flowseq_t *flowseq)
{
  free(flowseq);
}

void
tmap_fsw_flowseq_destroy(tmap_fsw_flowseq_t *flowseq)
{
  free(flowseq->base_calls);
  free(flowseq->flowgram);
  free(flowseq->flow_order);
  free(flowseq);
}

int64_t
tmap_fsw_extend_core(uint8_t *seq, int32_t len, 
                     tmap_fsw_flowseq_t *flowseq,
                     const tmap_fsw_param_t *ap,
                     int32_t right_j,
                     tmap_fsw_path_t *path, int32_t *path_len, int32_t prev_score)
{
  return prev_score + tmap_fsw_clipping_core(seq, len, flowseq,
                             ap, 0, 1, right_j, path, path_len);
}

int64_t
tmap_fsw_extend_fitting_core(uint8_t *seq, int32_t len, 
                             tmap_fsw_flowseq_t *flowseq,
                             const tmap_fsw_param_t *ap,
                             int32_t right_j,
                             tmap_fsw_path_t *path, int32_t *path_len, int32_t prev_score)
{
  return prev_score + tmap_fsw_clipping_core(seq, len, flowseq,
                             ap, 0, 0, right_j, path, path_len);
}

static inline
uint32_t tmap_fsw_path2cigar_get_type(uint32_t ctype)
{
  switch(ctype) {
    case TMAP_FSW_FROM_HP_PLUS: // deletion
      return TMAP_FSW_FROM_D;
      break;
    case TMAP_FSW_FROM_HP_MINUS: // insertion
      return TMAP_FSW_FROM_I;
      break;
    default:
      break;
  }
  return ctype;
}

uint32_t *
tmap_fsw_path2cigar(const tmap_fsw_path_t *path, int32_t path_len, int32_t *n_cigar, int32_t rm_hp)
{
  int32_t i, n;
  uint32_t *cigar;
  uint8_t dpscore_last_type, dpscore_cur_type;

  // Note: we could just use the function 'tmap_sw_path2cigar'

  if (path_len == 0 || path == 0) {
      *n_cigar = 0;
      return 0;
  }

  if(0 == rm_hp) {
      dpscore_last_type = path->ctype;
      for (i = n = 1; i < path_len; ++i) {
          if (dpscore_last_type != path[i].ctype) ++n;
          dpscore_last_type = path[i].ctype;
      }
      *n_cigar = n;
      cigar = tmap_malloc(*n_cigar * 4, "cigar");

      TMAP_SW_CIGAR_STORE(cigar[0], path[path_len-1].ctype, 1u);
      dpscore_last_type = path[path_len-1].ctype;
      for (i = path_len - 2, n = 0; i >= 0; --i) {
          if (path[i].ctype == dpscore_last_type) TMAP_SW_CIGAR_ADD_LENGTH(cigar[n], 1u);
          else {
              TMAP_SW_CIGAR_STORE(cigar[++n], path[i].ctype, 1u);
              dpscore_last_type = path[i].ctype;
          }
      }
  }
  else {
      // get the # of cigar operators
      dpscore_last_type = tmap_fsw_path2cigar_get_type(path[path_len-1].ctype);
      for(i = path_len - 2, n = 1;0 <= i; i--) {
          // get the current type
          dpscore_cur_type = tmap_fsw_path2cigar_get_type(path[i].ctype);
          if (dpscore_last_type != dpscore_cur_type) ++n;
          dpscore_last_type = dpscore_cur_type;
      }
      *n_cigar = n;
      cigar = tmap_malloc(*n_cigar * 4, "cigar");
          
      // get the last type
      dpscore_last_type = tmap_fsw_path2cigar_get_type(path[path_len-1].ctype);
      cigar[0] = 1u << 4 | dpscore_last_type;
      TMAP_SW_CIGAR_STORE(cigar[0], dpscore_last_type, 1u);
      for(i = path_len - 2, n = 0;0 <= i; i--) {
          // get the current type
          dpscore_cur_type = tmap_fsw_path2cigar_get_type(path[i].ctype);
          if (dpscore_cur_type == dpscore_last_type) TMAP_SW_CIGAR_ADD_LENGTH(cigar[n], 1u);
          else {
              TMAP_SW_CIGAR_STORE(cigar[++n], dpscore_cur_type, 1u);
              dpscore_last_type = dpscore_cur_type;
          }
      }
      if(n+1 != (*n_cigar)) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }
  }

  return cigar;
}

static int 
tmap_fsw_left_justify(char *ref, char *read, char *aln, int32_t len)
{
  char c;
  int32_t i, prev_del, prev_ins, start_ins, start_del, end_ins, end_del;
  int32_t justified = 0;

  prev_del = prev_ins = 0;
  start_del = start_ins = end_del = end_ins = -1;

  // reverse
  for(i=0;i<(len>>1);i++) {
      c = ref[i]; ref[i] = ref[len-i-1]; ref[len-i-1] = c;
      c = read[i]; read[i] = read[len-i-1]; read[len-i-1] = c;
      c = aln[i]; aln[i] = aln[len-i-1]; aln[len-i-1] = c;
  }

  for(i=0;i<len;) {
      if('-' == read[i]) { // deletion
          if(0 == prev_del) {
              start_del = i;
          }
          prev_del = 1;
          end_del = i;
          prev_ins = 0;
          start_ins = end_ins = -1;
          i++;
      }
      else if('-' == ref[i]) { // insert
          if(0 == prev_ins) {
              start_ins = i;
          }
          prev_ins = 1;
          end_ins = i;
          prev_del = 0;
          start_del = -1;
          end_del = -1;
          i++;
      }
      else {
          if(1 == prev_del) { // previous was an deletion
              start_del--;
              while(0 <= start_del && // bases remaining to examine 
                    read[start_del] != '-' && // hit another deletion 
                    ref[start_del] != '-' && // hit an insertion 
                    ref[start_del] == ref[end_del]) { // src ref base matches dest ref base 
                  // swap end_del and start_del for the read and aln
                  c = read[end_del]; read[end_del] = read[start_del]; read[start_del] = c;
                  c = aln[end_del]; aln[end_del] = aln[start_del]; aln[start_del] = c;
                  start_del--;
                  end_del--;
                  justified = 1;
              }
              end_del++; // we decremented when we exited the loop 
              i = end_del;
          }
          else if(1 == prev_ins) { // previous was an insertion
              start_ins--;
              while(0 <= start_ins && // bases remaining to examine
                    read[start_ins] != '-' && // hit another deletion 
                    ref[start_ins] != '-' && // hit an insertion 
                    read[start_ins] == read[end_ins]) { // src read base matches dest read base 
                  // swap end_ins and start_ins for the ref and aln
                  c = ref[end_ins]; ref[end_ins] = ref[start_ins]; ref[start_ins] = c;
                  c = aln[end_ins]; aln[end_ins] = aln[start_ins]; aln[start_ins] = c;
                  start_ins--;
                  end_ins--;
                  justified = 1;
              }
              end_ins++; // e decremented when we exited the loop 
              i = end_ins;
          }
          else {
              i++;
          }
          // reset
          prev_del = prev_ins = 0;
          start_del = start_ins = end_del = end_ins = -1;
      }
  }

  // reverse
  for(i=0;i<(len>>1);i++) {
      c = ref[i]; ref[i] = ref[len-i-1]; ref[len-i-1] = c;
      c = read[i]; read[i] = read[len-i-1]; read[len-i-1] = c;
      c = aln[i]; aln[i] = aln[len-i-1]; aln[len-i-1] = c;
  }

  return justified;
}

static int 
tmap_fsw_right_justify(char *ref, char *read, char *aln, int32_t len)
{
  int32_t i;
  int32_t justified = 0;
  // reverse
  for(i=0;i<(len>>1);i++) {
      char tmp;
      tmp = ref[i]; ref[i] = ref[len-i-1]; ref[len-i-1] = tmp;
      tmp = read[i]; read[i] = read[len-i-1]; read[len-i-1] = tmp;
      tmp = aln[i]; aln[i] = aln[len-i-1]; aln[len-i-1] = tmp;
  }
  // left-justify
  justified = tmap_fsw_left_justify(ref, read, aln, len);
  // reverse
  for(i=0;i<(len>>1);i++) {
      char tmp;
      tmp = ref[i]; ref[i] = ref[len-i-1]; ref[len-i-1] = tmp;
      tmp = read[i]; read[i] = read[len-i-1]; read[len-i-1] = tmp;
      tmp = aln[i]; aln[i] = aln[len-i-1]; aln[len-i-1] = tmp;
  }
  return justified;
}

void
tmap_fsw_get_aln(tmap_fsw_path_t *path, int32_t path_len,
                 uint8_t *flow_order, int32_t flow_order_len, uint8_t *target, uint8_t strand,
                 char **ref, char **read, char **aln, int32_t j_type)
{
  int32_t i, j;

  (*ref) = tmap_malloc(sizeof(char) * (1 + path_len), "(*ref)");
  (*read) = tmap_malloc(sizeof(char) * (1 + path_len), "(*read)");
  (*aln) = tmap_malloc(sizeof(char) * (1 + path_len), "(*aln)");

  // ref
  for(i=path_len-1,j=0;0<=i;i--,j++) {
      if(TMAP_FSW_FROM_M == path[i].ctype 
         || TMAP_FSW_FROM_D == path[i].ctype
         || TMAP_FSW_FROM_HP_PLUS == path[i].ctype) {
          (*ref)[j] = "ACGTN"[target[path[i].j]];
      }
      else {
          (*ref)[j] = '-';
      }
  }
  (*ref)[path_len] = '\0';

  // alignment string
  for(i=path_len-1,j=0;0<=i;i--,j++) {
      switch(path[i].ctype) {
        case TMAP_FSW_FROM_M:
          if(flow_order[path[i].i % flow_order_len] == target[path[i].j]) {
              (*aln)[j] = '|';
          }
          else { 
              (*aln)[j] = ' ';
          }
          break;
        case TMAP_FSW_FROM_I:
          (*aln)[j] = '+'; break;
        case TMAP_FSW_FROM_D:
          (*aln)[j] = '-'; break;
        case TMAP_FSW_FROM_HP_PLUS: // overcall
          (*aln)[j] = 'h'; break;
        case TMAP_FSW_FROM_HP_MINUS: // undercall
          (*aln)[j] = 'H'; break;
        default:
          (*aln)[j] = ' '; break;
      }
  }
  (*aln)[path_len] = '\0';

  // read
  for(i=path_len-1,j=0;0<=i;i--,j++) {
      if(TMAP_FSW_FROM_M == path[i].ctype 
         || TMAP_FSW_FROM_I == path[i].ctype
         || TMAP_FSW_FROM_HP_MINUS == path[i].ctype) {
          (*read)[j] = "ACGTN"[flow_order[path[i].i % flow_order_len]];
      }
      else {
          (*read)[j] = '-';
      }
  }
  (*read)[path_len] = '\0';

  // the input is assumed to be in read-order
  switch(j_type) {
    case TMAP_FSW_NO_JUSTIFY:
      // do nothing
      break;
    case TMAP_FSW_JUSTIFY_LEFT_REF:
      if(0 == strand) { 
          tmap_fsw_left_justify((*ref), (*read), (*aln), path_len);
      }
      else {
          tmap_fsw_right_justify((*ref), (*read), (*aln), path_len);
      }
      break;
    case TMAP_FSW_JUSTIFY_LEFT_READ:
      tmap_fsw_left_justify((*ref), (*read), (*aln), path_len);
      break;
    default:
      break;
  }
}

void 
tmap_fsw_print_aln(tmap_file_t *fp, int64_t score, tmap_fsw_path_t *path, int32_t path_len,
                   uint8_t *flow_order, int32_t flow_order_len, uint8_t *target, uint8_t strand, int32_t j_type, char sep)
{
  char *ref=NULL, *read=NULL, *aln=NULL;

  tmap_fsw_get_aln(path, path_len, flow_order, flow_order_len, target, strand, &ref, &read, &aln, j_type);

  tmap_file_fprintf(fp, "%lld%c%s%c%s%c%s",
                    (long long int)score, sep, read, sep, aln, sep, ref);

  free(ref);
  free(read);
  free(aln);
}

tmap_fsw_flowseq_t *
tmap_fsw_flowseq_from_seq(tmap_fsw_flowseq_t *fs, tmap_seq_t *seq, uint8_t *flow_order, int32_t flow_order_len, uint8_t *key_seq, uint8_t key_seq_len, int32_t use_flowgram)
{

  int32_t i, j;
  tmap_string_t *bases = NULL;
  int32_t last_key_flow;
  int32_t to_fill = 0, was_int;
  int32_t num_flows = 0;

  if(NULL == fs) {
      fs = tmap_calloc(1, sizeof(tmap_fsw_flowseq_t), "fs");
  }
  
  // base sequence
  bases = tmap_seq_get_bases(seq);
  if(0 == bases->l) tmap_error("bug encountered", Exit, OutOfRange);

  // convert bases to integers
  was_int = tmap_seq_is_int(seq);
  if(0 == tmap_seq_is_int(seq)) {
      tmap_seq_to_int(seq);
  }

  /*
  tmap_print_debug_int((int)bases->l);
  for(i=0;i<bases->l;i++) {
      fputc("ACGTN"[(int)bases->s[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<flow_order_len;i++) {
      fputc("ACGTN"[(int)flow_order[i]], stderr);
  }
  fputc('\n', stderr);
  if(0 < key_seq_len) {
      for(i=0;i<key_seq_len;i++) {
          fputc("ACGTN"[(int)key_seq[i]], stderr);
      }
      fputc('\n', stderr);
  }
  */
  
  // key bases
  if(0 < key_seq_len && key_seq[key_seq_len-1] == bases->s[0]) {
      fs->key_index = 0;
      fs->key_bases = 0;
      for(i=key_seq_len-1;0<=i && key_seq[i] == bases->s[0];i--) {
          fs->key_bases++;
      }
  }
  else {
      fs->key_index = -1;
      fs->key_bases = 0;
  }
      
  // get the flow for the last base in the key
  i = 0;
  if(0 < key_seq_len) {
      last_key_flow = -1;
      while(i < key_seq_len) {
          last_key_flow++;
          while(key_seq[i] != flow_order[last_key_flow % flow_order_len]) {
              last_key_flow++;
          }
          while(i < key_seq_len && key_seq[i] == flow_order[last_key_flow]) {
              i++;
          }
      }
  }
  else {
      last_key_flow = 0;
  }

  // flow order
  // NB: shift based on the key sequence
  to_fill = 0;
  if(flow_order_len == fs->flow_order_len) {
      j = last_key_flow;
      if(0 < key_seq_len && 0 == fs->key_bases) j++; 
      for(i=0;i<flow_order_len;i++,j++) {
          if(fs->flow_order[i] != flow_order[j % flow_order_len]) {
              to_fill = 1;
              break;
          }
      }
  }
  else {
      to_fill = 1;
  }
  if(1 == to_fill) {
      fs->flow_order = tmap_realloc(fs->flow_order, sizeof(uint8_t) * flow_order_len, "fs->flow_order");
      j = last_key_flow;
      if(0 < key_seq_len && 0 == fs->key_bases) j++; 
      for(i=0;i<flow_order_len;i++,j++) {
          fs->flow_order[i] = flow_order[j % flow_order_len];
      }
      fs->flow_order_len = flow_order_len;
  }

  // find the flowgram length
  i = num_flows = 0;
  while(i < bases->l) {
      while(bases->s[i] != fs->flow_order[num_flows % flow_order_len]) {
          num_flows++;
      }
      while(i < bases->l && bases->s[i] == fs->flow_order[num_flows % flow_order_len]) {
          i++;
      }
      num_flows++;
  }

  // get flowgram
  // NB: the returned flowgram should always include the key sequence
  to_fill = 0;
  if(1 == use_flowgram) {
      fs->num_flows = tmap_seq_get_flowgram(seq, &fs->flowgram, fs->mem);
  }
  else {
      fs->num_flows = 0;
  }
  if(0 == fs->num_flows) { // unsuccessful
      // will generate one from the base sequence
      to_fill = 1;
      
      // flowgram
      fs->num_flows = num_flows; 
      if(fs->mem < fs->num_flows) {
          fs->flowgram = tmap_calloc(fs->num_flows, sizeof(uint16_t), "flowgram");
      }
  }
  else if(0 < key_seq_len) {
      // remove the key from the flowgram
      j = last_key_flow + 1;
      if(0 < fs->key_bases) {
          j--;
      }
      for(i=0;i<fs->num_flows-j;i++) {
          fs->flowgram[i] = fs->flowgram[i+j];
      }
      fs->num_flows -= j;
      // NB: we must remove ending flows that are not part of the base sequence
      if(num_flows < fs->num_flows) {
          fs->num_flows = num_flows;
          fs->flowgram = tmap_realloc(fs->flowgram, sizeof(uint16_t) * fs->num_flows, "fs->flowgram");
      }
      else if(fs->num_flows < num_flows) {
          tmap_print_debug_int(fs->num_flows);
          tmap_print_debug_int(num_flows);
          tmap_print_debug_int(j);
          tmap_print_debug_string(tmap_seq_get_name(seq)->s);
          tmap_print_debug_int(((int)tmap_seq_get_bases(seq)->l));
          tmap_error("bug encountered", Exit, OutOfRange);
      }
      // NB: do not include the key sequence in the key bases
  }

  // base calls and memory
  if(fs->mem < fs->num_flows) {
      fs->base_calls = tmap_realloc(fs->base_calls, sizeof(uint8_t) * fs->num_flows, "base_calls");
      fs->mem = fs->num_flows;
  }

  // flowgram and bases
  i = 0;
  j = 0;
  // last key base and first base
  while(i < bases->l) {
      while(bases->s[i] != fs->flow_order[j % flow_order_len]) {
          if(fs->num_flows <= j) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }
          if(1 == to_fill) {
              fs->flowgram[j] = 0;
          }
          fs->base_calls[j] = 0;
          j++;
      }
      if(1 == to_fill) {
          fs->flowgram[j] = 0;
      }
      fs->base_calls[j] = 0;
      while(i < bases->l && bases->s[i] == fs->flow_order[j % flow_order_len]) {
          if(fs->num_flows <= j) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }
          if(1 == to_fill) {
              fs->flowgram[j] += 100;
          }
          fs->base_calls[j]++;
          i++;
      }
      j++;
  }
  if(j != fs->num_flows) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }
  // key sequence in the flowgram
  if(1 == to_fill) {
      fs->flowgram[0] += 100 * fs->key_bases;
  }
  // NB: do not include the key sequence in the key bases

  if(0 == was_int) {
      tmap_seq_to_char(seq);
  }

  return fs;
}

#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
typedef struct {
    char *flow_order;
    char *base_calls;
    char *flowgram;
    int32_t num_flows;
    char *target;
    int32_t target_length;
    tmap_fsw_param_t param;
    int32_t j_type;
} tmap_fsw_main_opt_t;

static tmap_fsw_main_opt_t *
tmap_fsw_main_opt_init()
{
  tmap_fsw_main_opt_t *opt = NULL;
  opt = tmap_calloc(1, sizeof(tmap_fsw_main_opt_t), "opt");

  opt->flow_order = tmap_strdup("TAGC");
  opt->base_calls = NULL;
  opt->flowgram = NULL;
  opt->target = NULL;
  opt->target_length = 0;
  opt->j_type = TMAP_FSW_NO_JUSTIFY;

  // param
  opt->param.gap_open = TMAP_MAP_OPT_PEN_GAPO*100;
  opt->param.gap_ext = TMAP_MAP_OPT_PEN_GAPE*100;
  opt->param.gap_end = TMAP_MAP_OPT_PEN_GAPE*100;
  opt->param.matrix = tmap_fsw_sm_short; // 11*100 - match, -19*100 - mismatch
  opt->param.fscore = TMAP_MAP_OPT_FSCORE*100; // set this to score_match + gap_open + gap_ext
  opt->param.offset = 0;
  opt->param.row = 5;
  opt->param.band_width = 50;

  return opt;
}

static void
tmap_fsw_main_opt_destroy(tmap_fsw_main_opt_t *opt)
{
  free(opt->flow_order);
  free(opt->base_calls);
  free(opt->flowgram);
  free(opt->target);
  free(opt);
}

static uint16_t*
tmap_fsw_main_get_flowgram(char *flowgram, int32_t *num_flows)
{
  int32_t i, j, m;
  uint16_t *f = NULL;
  int32_t state = 0;
  int32_t v;

  // state  |  todo
  // 0      |  read in before decimal
  // 1      |  read in decimal or comma
  // 2      |  read in after decimal
  // 3      |  read in comma

  m = 4;
  f = tmap_calloc(m, sizeof(uint16_t), "f");

  i=j=0;
  while('\0' != flowgram[i]) {
      //fprintf(stderr, "state=%d flowgram[%d]=\"%s\" j=%d\n", state, i, flowgram+i, j);
      if('.' == flowgram[i]) {
          if(1 != state) {
              tmap_error("unexpected decimal point", Exit, OutOfRange);
          }
          i++; // skip over the decimal
          state = 2;
      }
      else if(',' == flowgram[i]) {
          if(1 != state && 3 != state) {
              tmap_error("unexpected comma", Exit, OutOfRange);
          }
          i++; // skip over the comma
          j++; // next flowgram
          state = 0;
      }
      else if(0 != isgraph(flowgram[i])) {
          if(0 != state && 2 != state) {
              tmap_error("unexpected digits", Exit, OutOfRange);
          }
          if(m <= j) { // more memory
              m <<= 1;
              f = tmap_realloc(f, m*sizeof(uint16_t), "f");
          }
          if(0 == state) {
              v = atoi(flowgram + i);
              f[j] = 100*v;
          }
          else { // state 2
              int32_t x = 10;
              // keep only the 2 left-most digits
              for(v=0;v<2 && '0'<=flowgram[i+v] && flowgram[i+v] <= '9';v++) {
                  f[j] += x * ((int)(flowgram[i+v]) - (int)'0');
                  x /= 10;
              }
          }
          while('\0' != flowgram[i] && '0' <= flowgram[i] && flowgram[i] <= '9') { // skip over
              i++;
          }
          state++;
      }
      else {
          tmap_error("could not parse the flowgram", Exit, OutOfRange);
      }
  }
  j++;
  if(1 != state && 3 != state) {
      tmap_error("ended with a comma or a decimal point", Exit, OutOfRange);
  } 

  (*num_flows) = j;
  f = tmap_realloc(f, (*num_flows)*sizeof(uint16_t), "f");
  return f;
}

uint8_t *
tmap_fsw_main_get_base_calls(char *optarg, int32_t num_flows, char *flow_order)
{
  int32_t i, j, k;
  uint8_t *base_calls = NULL;
  int32_t flow_order_len = strlen(flow_order);

  base_calls = tmap_calloc(num_flows, sizeof(uint8_t), "base_calls");

  for(i=j=0;i<strlen(optarg);i++) {
      k=0;
      while(toupper(flow_order[j % flow_order_len]) != toupper(optarg[i])) {
          j++; k++;
          if(3 < k) tmap_error("could not understand the base", Exit, OutOfRange);
      }
      base_calls[j]++;
  }
  j++;
  while(j < num_flows) {
      base_calls[j]=0;
      j++;
  }

  return base_calls;
}

static int32_t
usage(tmap_fsw_main_opt_t *opt)
{
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Usage: %s fsw [options]", PACKAGE);
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Options (required):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -b STRING   base calls ([ACGT]+) [%s]\n",
                    opt->base_calls);
  tmap_file_fprintf(tmap_file_stderr, "         -f STRING   flowgram float values (comma separated) [%s]\n",
                    opt->flowgram);
  tmap_file_fprintf(tmap_file_stderr, "         -t STRING   target sequence ([ACGT]+) [%s]\n",
                    opt->target);
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Options (optional):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -x STRING   the flow order ([ACGT]+) [%s]\n",
                    opt->flow_order);
  tmap_file_fprintf(tmap_file_stderr, "         -X INT      flow penalty [%d]\n", opt->param.fscore);
  tmap_file_fprintf(tmap_file_stderr, "         -o INT      search for homopolymer errors +- offset [%d]\n", opt->param.offset);
  tmap_file_fprintf(tmap_file_stderr, "         -l INT      indel justification type: 0 - none, 1 - 5' strand of the reference, 2 - 5' strand of the read [%d]\n", opt->j_type);
  tmap_file_fprintf(tmap_file_stderr, "         -v          print verbose progress information\n");
  tmap_file_fprintf(tmap_file_stderr, "         -h          print this message\n");
  tmap_file_fprintf(tmap_file_stderr, "\n");

  return 1;
}

int tmap_fsw_main(int argc, char *argv[])
{
  int32_t i, c;
  tmap_fsw_path_t *path=NULL;
  int32_t path_len = 0;
  tmap_fsw_main_opt_t *opt;
  int32_t num_flows = 0;
  uint16_t *flowgram = NULL;
  uint8_t *base_calls = NULL;
  uint8_t *flow_order=NULL;
  int32_t flow_order_len;
  tmap_fsw_flowseq_t *flowseq = NULL;

  opt = tmap_fsw_main_opt_init();

  while((c = getopt(argc, argv, "b:f:t:x:X:o:l:vh")) >= 0) {
      switch(c) {
        case 'b':
          opt->base_calls = tmap_strdup(optarg); break;
        case 'f':
          opt->flowgram = tmap_strdup(optarg); break;
        case 't':
          opt->target = tmap_strdup(optarg); opt->target_length = strlen(opt->target); break;
        case 'x':
          free(opt->flow_order);
          opt->flow_order = tmap_strdup(optarg); break;
        case 'X':
          opt->param.fscore = 100*atoi(optarg); break;
        case 'o':
          opt->param.offset = atoi(optarg); break;
        case 'l':
          opt->j_type = atoi(optarg); break;
        case 'v':
          tmap_progress_set_verbosity(1); break;
        case 'h':
        default:
          return usage(opt);
      }
  }

  if(argc != optind || 1 == argc) {
      return usage(opt);
  }
  else { // check the command line options
      if(NULL == opt->base_calls) {
          tmap_error("option -b must be specified", Exit, CommandLineArgument);
      }
      if(NULL == opt->flowgram) {
          tmap_error("option -f must be specified", Exit, CommandLineArgument);
      }
      if(NULL == opt->target) {
          tmap_error("option -t must be specified", Exit, CommandLineArgument);
      }
      tmap_validate_flow_order(opt->flow_order);
      tmap_error_cmd_check_int(opt->param.fscore, 0, INT32_MAX, "-X");
      tmap_error_cmd_check_int(opt->param.offset, 0, INT32_MAX, "-o");
      tmap_error_cmd_check_int(opt->param.offset, TMAP_FSW_NO_JUSTIFY, TMAP_FSW_NO_JUSTIFY, "-l");
  }

  // get the flowgram
  flowgram = tmap_fsw_main_get_flowgram(opt->flowgram, &num_flows); 

  // get the flow base calls
  base_calls = tmap_fsw_main_get_base_calls(opt->base_calls, num_flows, opt->flow_order);

  // convert the target to 2-bit format 
  for(i=0;i<opt->target_length;i++) {
      opt->target[i] = tmap_nt_char_to_int[(int)opt->target[i]];
  }

  flow_order_len = strlen(opt->flow_order);
  flow_order = tmap_malloc(sizeof(uint8_t) * flow_order_len, "flow_order"); 
  for(i=0;i<flow_order_len;i++) {
      flow_order[i] = tmap_nt_char_to_int[(int)opt->flow_order[i]];
  }

  path = tmap_malloc(sizeof(tmap_fsw_path_t)*TMAP_FSW_MAX_PATH_LENGTH(opt->target_length, num_flows, opt->param.offset), "path");

  flowseq = tmap_fsw_flowseq_init(flow_order, flow_order_len, base_calls, flowgram, num_flows, -1, -1);

  // align
  int64_t score = tmap_fsw_clipping_core((uint8_t*)opt->target, opt->target_length,
                                       flowseq,
                                       &opt->param,
                                       0, 0, 0,
                                       path, &path_len);

  tmap_file_stdout = tmap_file_fdopen(fileno(stdout), "wb", TMAP_FILE_NO_COMPRESSION);

  // print
  tmap_fsw_print_aln(tmap_file_stdout, score, path, path_len,
                     flow_order, flow_order_len, (uint8_t*)opt->target, 0, opt->j_type, '\t');
  tmap_file_fprintf(tmap_file_stdout, "\n");

  tmap_file_fclose(tmap_file_stdout);

  // destroy memory
  tmap_fsw_main_opt_destroy(opt);
  free(flowgram);
  free(base_calls);
  free(path);
  free(flow_order);
  tmap_fsw_flowseq_destroy_shallow(flowseq);

  return 0;
}
#endif
