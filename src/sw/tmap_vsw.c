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
#include <smmintrin.h>
#include <unistd.h>
#include <stdio.h>
#include "../util/tmap_alloc.h"
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
tmap_vsw_query_init()
{
  tmap_vsw_query_t *query;
  query = tmap_calloc(1, sizeof(tmap_vsw_query_t), "query");
  return query;
}

void
tmap_vsw_query_destroy(tmap_vsw_query_t *query)
{
  if(NULL == query) return;
  if(NULL != query->query16) tmap_vsw16_query_destroy(query->query16);
  free(query);
}

int32_t
tmap_vsw_sse2(tmap_vsw_query_t *q,
              const uint8_t *query, int32_t qlen, 
              const uint8_t *target, int32_t tlen, 
              tmap_vsw_opt_t *opt, tmap_vsw_result_t *result,
              int32_t *overflow)
{
  // TODO: check potential overflow
  // TODO: check that gap penalties will not result in an overflow
  // TODO: check that the max/min alignment score do not result in an overflow

  q->query16 = tmap_vsw16_query_init_short(q->query16, query, qlen, 0, 0, opt); 
  tmap_vsw16_sse2(q->query16, target, tlen, opt, result, overflow);
  if(1 == *overflow) {
      return INT32_MIN; 
  }
  return result->score_fwd;
}

void
tmap_vsw_sse2_get_path(const uint8_t *query, int32_t qlen,
                        const uint8_t *target, int32_t tlen,
                        tmap_vsw_query_t *q,
                        tmap_vsw_result_t *result,
                        tmap_sw_path_t *path,
                        int32_t *path_len,
                        tmap_vsw_opt_t *opt)
{
  q->query16 = tmap_vsw16_query_init_full(q->query16, query, qlen, tlen, 0, 0, opt); 
  tmap_vsw16_sse2_get_path(query, qlen, target, tlen,
                           q->query16, result, 
                           path, path_len,
                           opt);
}

