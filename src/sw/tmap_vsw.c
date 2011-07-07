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
#include "tmap_vsw16.h"
#include "tmap_vsw.h"

tmap_vsw_opt_t*
tmap_vsw_opt_init(int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape, int32_t score_thr)
{
  tmap_vsw_opt_t *opt;

  opt = tmap_malloc(sizeof(tmap_vsw_opt_t), "opt");
  opt->score_match = score_match; 
  opt->pen_mm = pen_mm;
  opt->pen_gapo = pen_gapo;
  opt->pen_gape = pen_gape;
  opt->score_thres = score_thr;

  return opt;
}

void
tmap_vsw_opt_destroy(tmap_vsw_opt_t *opt)
{
  free(opt);
}

tmap_vsw_result_t *
tmap_vsw_result_init()
{
  return tmap_calloc(1, sizeof(tmap_vsw_result_t), "return");
}

void
tmap_vsw_result_destroy(tmap_vsw_result_t *result)
{
  if(NULL == result) return;
  free(result);
}

tmap_vsw_query_t*
tmap_vsw_query_init(const uint8_t *query_fwd, int32_t qlen_fwd,
                    const uint8_t *query_rev, int32_t qlen_rev,
                    int32_t tlen,
                    int32_t query_start_clip, int32_t query_end_clip,
                    tmap_vsw_opt_t *opt)
{
  tmap_vsw_query_t *query;
  query = tmap_calloc(1, sizeof(tmap_vsw_query_t), "query");
  query->query16_fwd = tmap_vsw16_query_init(query->query16_fwd, query_fwd, qlen_fwd, tlen, query_start_clip, query_end_clip, opt);
  query->query16_rev = tmap_vsw16_query_init(query->query16_rev, query_rev, qlen_rev, tlen, query_end_clip, query_start_clip, opt);
  return query;
}

void
tmap_vsw_query_destroy(tmap_vsw_query_t *query)
{
  if(NULL == query) return;
  if(NULL != query->query16_fwd) tmap_vsw16_query_destroy(query->query16_fwd);
  if(NULL != query->query16_rev) tmap_vsw16_query_destroy(query->query16_rev);
  free(query);
}

int32_t
tmap_vsw_sse2(tmap_vsw_query_t *query,
              const uint8_t *query_fwd, int32_t qlen_fwd,
              const uint8_t *query_rev, int32_t qlen_rev,
              uint8_t *target, int32_t tlen, 
              int32_t query_start_clip, int32_t query_end_clip,
              tmap_vsw_opt_t *opt, tmap_vsw_result_t *result,
              int32_t *overflow, int32_t score_thr)
{
  uint8_t *tmp_target;
  int32_t tmp_tlen, found_forward = 1;
#ifdef TMAP_VSW_DEBUG
  int32_t i;
#endif
  // TODO: check potential overflow
  // TODO: check that gap penalties will not result in an overflow
  // TODO: check that the max/min alignment score do not result in an overflow
  

#ifdef TMAP_VSW_DEBUG
  fprintf(stderr, "in %s\n", __func__);
  fprintf(stderr, "query_start_clip=%d\n", query_start_clip);
  fprintf(stderr, "query_end_clip=%d\n", query_end_clip);
  for(i=0;i<qlen_rev;i++) {
      fputc("ACGTN"[query_fwd[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<tlen;i++) {
      fputc("ACGTN"[target[i]], stderr);
  }
  fputc('\n', stderr);
#endif

  // init forward
  query->query16_fwd = tmap_vsw16_query_init(query->query16_fwd, query_fwd, qlen_fwd, tlen, query_start_clip, query_end_clip, opt);

  // run forward
  result->score_fwd = tmap_vsw16_sse2_forward(query->query16_fwd, target, tlen,
                                              query_start_clip, query_end_clip,
                                              opt, &result->query_end, &result->target_end,
                                              0, overflow, score_thr);

  // check forward results
  if(NULL != overflow && 1 == *overflow) {
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
      tmap_error("bug encountered", Exit, OutOfRange);
  }

  // return if we found no legal/good forward results
  if(0 == found_forward) {
      result->query_end = result->query_start = 0;
      result->target_end = result->target_start = 0;
      result->score_fwd = result->score_rev = INT32_MIN;
      return INT32_MIN;
  }

  tmp_target = target;
  tmp_tlen = tlen;
  // reverse compliment
  tmap_reverse_compliment_int(tmp_target, tmp_tlen);

  // adjust based on the forward
  query_rev += qlen_rev - result->query_end - 1;
  qlen_rev = result->query_end + 1;
  target += tlen - result->target_end - 1;
  tlen = result->target_end + 1;

#ifdef TMAP_VSW_DEBUG
  fprintf(stderr, "{?-%d] {?-%d}\n",
          result->query_end,
          result->target_end);
  fprintf(stderr, "qlen_rev=%d tlen=%d\n", qlen_rev, tlen);
  for(i=0;i<qlen_rev;i++) {
      fputc("ACGTN"[query_rev[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<tlen;i++) {
      fputc("ACGTN"[target[i]], stderr);
  }
  fputc('\n', stderr);
#endif

  // init reverse
  // NB: reverse the clipping
  query->query16_rev = tmap_vsw16_query_init(query->query16_rev, query_rev, qlen_rev, tlen, query_end_clip, query_start_clip, opt);

  // run the reverse
  // NB: we do not start clip since we know here the alignment should begin on
  // the reverse
  result->score_rev = tmap_vsw16_sse2_forward(query->query16_rev, target, tlen,
                                              0, query_start_clip,
                                              opt, &result->query_start, &result->target_start, 
                                              1, overflow, score_thr);

#ifdef TMAP_VSW_DEBUG
  fprintf(stderr, "result->score_fwd=%d result->score_rev=%d\n",
          result->score_fwd, result->score_rev);
#endif

  // check reverse results
  if(1 == *overflow) {
      result->query_end = result->query_start = 0;
      result->target_end = result->target_start = 0;
      result->score_fwd = result->score_rev = INT32_MIN;
      return INT32_MIN; 
  }
  else if(result->score_fwd != result->score_rev) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }

#ifdef TMAP_VSW_DEBUG
  fprintf(stderr, "{%d-%d] {%d-%d}\n",
          result->query_start, result->query_end,
          result->target_start, result->target_end);
#endif

  // adjust query_start and query_end
  result->query_start = qlen_rev - result->query_start - 1;
  result->target_start = tlen - result->target_start - 1;
  
#ifdef TMAP_VSW_DEBUG
  fprintf(stderr, "{%d-%d] {%d-%d}\n",
          result->query_start, result->query_end,
          result->target_start, result->target_end);
#endif
  // reverse compliment (is this necessary?)
  tmap_reverse_compliment_int(tmp_target, tmp_tlen);

  // return
  return result->score_fwd;
}
