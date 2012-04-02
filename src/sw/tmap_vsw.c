/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/* The MIT License

   Copyright (c) 2011 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   */

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <unistd.h>
#include <stdio.h>
#include <config.h>
#include "../util/tmap_alloc.h"
#include "../util/tmap_error.h"
#include "../util/tmap_definitions.h"
#include "tmap_sw.h"
#include "tmap_vsw_definitions.h"
#include "tmap_vsw_s0.h"
#include "tmap_vsw.h"

//#define TMAP_VSW_DEBUG

tmap_vsw_t*
tmap_vsw_init(const uint8_t *query, int32_t qlen,
                    int32_t query_start_clip, int32_t query_end_clip,
                    int32_t type,
                    tmap_vsw_opt_t *opt)
{
  tmap_vsw_t *vsw = NULL;
  vsw = tmap_calloc(1, sizeof(tmap_vsw_t), "vsw");
  vsw->type = type;
  vsw->query_start_clip = query_start_clip;
  vsw->query_end_clip = query_end_clip;
  vsw->opt = opt;
  switch(type) {
    case TMAP_VSW_TYPE_S0:
      vsw->data.s0 = tmap_vsw_data_init_s0(query, qlen, query_start_clip, query_end_clip, opt);
      break;
    case TMAP_VSW_TYPE_S1:
    default:
      break;
  }
  return vsw;
}

void
tmap_vsw_destroy(tmap_vsw_t *vsw)
{
  if(NULL == vsw) return;
  switch(vsw->type) {
    case TMAP_VSW_TYPE_S0:
      tmap_vsw_data_destroy_s0(vsw->data.s0);
      vsw->data.s0 = NULL;
      break;
    case TMAP_VSW_TYPE_S1:
    default:
      break;
  }
  free(vsw);
}

tmap_vsw_t*
tmap_vsw_update(tmap_vsw_t *vsw, const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen)
{
  switch(vsw->type) {
    case TMAP_VSW_TYPE_S0:
      tmap_vsw_data_update_s0(vsw->data.s0, query, qlen, target, tlen);
      break;
    case TMAP_VSW_TYPE_S1:
    default:
      break;
  }
  return vsw;
}

int32_t
tmap_vsw_process(tmap_vsw_t *vsw,
              const uint8_t *query, int32_t qlen,
              uint8_t *target, int32_t tlen, 
              tmap_vsw_result_t *result,
              int32_t *overflow, int32_t score_thr, int32_t is_rev)
{
  int32_t found_forward = 1, query_end, target_end, n_best, score;
#ifdef TMAP_VSW_DEBUG
  int32_t i;
#endif
  // TODO: check potential overflow
  // TODO: check that gap penalties will not result in an overflow
  // TODO: check that the max/min alignment score do not result in an overflow

#ifdef TMAP_VSW_DEBUG
  fprintf(stderr, "in %s is_rev=%d\n", __func__, is_rev);
  fprintf(stderr, "query_start_clip=%d\n", vsw->query_start_clip);
  fprintf(stderr, "query_end_clip=%d\n", vsw->query_end_clip);
  fprintf(stderr, "qlen=%d tlen=%d\n", qlen, tlen);
  for(i=0;i<qlen;i++) {
      fputc("ACGTN"[query[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<tlen;i++) {
      fputc("ACGTN"[target[i]], stderr);
  }
  fputc('\n', stderr);
#endif
      
  // update
  tmap_vsw_update(vsw, query, qlen, target, tlen);

  query_end = target_end = n_best = 0;
      
  // perform the alignment
  switch(vsw->type) {
    case TMAP_VSW_TYPE_S0:
      score = tmap_vsw_process_s0(vsw->data.s0, 
                          query, qlen, target, tlen, 
                          vsw->query_start_clip, vsw->query_end_clip, vsw->opt,
                          is_rev, score_thr, 
                          &query_end, &target_end, &n_best, overflow);
      break;
    case TMAP_VSW_TYPE_S1:
    default:
      score = INT32_MIN;
      break;
  }
  
  if(0 == is_rev) {

      result->query_end = query_end;
      result->target_end = target_end;
      result->n_best = n_best;
      result->score_fwd = score;

      // check forward results
      if(NULL != overflow && 1 == (*overflow)) {
          found_forward = 0;
      }
      else if(result->score_fwd < score_thr) {
          found_forward = 0;
      }
      else if((result->query_end == result->query_start || result->target_end == result->target_start)
              && result->score_fwd <= 0) {
          found_forward = 0;
      }
      else if(-1 == result->query_end) {
          tmap_bug();
      }

      // return if we found no legal/good forward results
      if(0 == found_forward) {
          result->query_end = result->query_start = 0;
          result->target_end = result->target_start = 0;
          result->score_fwd = result->score_rev = INT16_MIN;
          result->n_best = 0;
          return INT32_MIN;
      }

      result->query_start = result->target_start = 0;
      result->score_rev = INT16_MIN;

#ifdef TMAP_VSW_DEBUG
      fprintf(stderr, "result->score_fwd=%d result->score_rev=%d\n",
              result->score_fwd, result->score_rev);
      fprintf(stderr, "{?-%d] {?-%d}\n",
              result->query_end,
              result->target_end);
#endif

      return result->score_fwd;
  }
  else {

      result->query_start = qlen - query_end - 1;
      result->target_start = tlen - target_end - 1;
      result->n_best = n_best;
      result->score_rev = score;

#ifdef TMAP_VSW_DEBUG
      fprintf(stderr, "is_rev=%d result->score_fwd=%d result->score_rev=%d\n",
              is_rev, result->score_fwd, result->score_rev);
#endif

      // check reverse results
      if(1 == (*overflow)) {
          result->query_end = result->query_start = 0;
          result->target_end = result->target_start = 0;
          result->score_fwd = result->score_rev = INT16_MIN;
          result->n_best = 0;
          return INT32_MIN; 
      }
      else if(result->score_fwd != result->score_rev) {
          fprintf(stderr, "{%d-%d} {%d-%d}\n",
                  result->query_start, result->query_end,
                  result->target_start, result->target_end);
          fprintf(stderr, "result->score_fwd=%d result->score_rev=%d\n",
                  result->score_fwd, result->score_rev);
          tmap_bug();
      }

      // return
      return result->score_fwd;
  }
}
