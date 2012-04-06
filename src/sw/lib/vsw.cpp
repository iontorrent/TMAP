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
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include "vsw16.h"
#include "vsw.h"

vsw_opt_t*
vsw_opt_init(int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape, int32_t score_thr)
{
  vsw_opt_t *opt;

  opt = (vsw_opt_t*)malloc(sizeof(vsw_opt_t));
  opt->score_match = score_match; 
  opt->pen_mm = pen_mm;
  opt->pen_gapo = pen_gapo;
  opt->pen_gape = pen_gape;
  opt->score_thres = score_thr;

  return opt;
}

void
vsw_opt_destroy(vsw_opt_t *opt)
{
  free(opt);
}

vsw_query_t*
vsw_query_init(const uint8_t *query, int32_t qlen,
                    int32_t tlen,
                    int32_t query_start_clip, int32_t query_end_clip,
                    vsw_opt_t *opt)
{
  vsw_query_t *vsw_query;
  vsw_query = (vsw_query_t*)calloc(1, sizeof(vsw_query_t));
  vsw_query->query16 = vsw16_query_init(vsw_query->query16, query, qlen, query_start_clip, query_end_clip, opt);
  // TMP
  vsw_query->query = query;
  vsw_query->qlen = qlen;
  return vsw_query;
}

void
vsw_query_destroy(vsw_query_t *query)
{
  if(NULL == query) return;
  if(NULL != query->query16) vsw16_query_destroy(query->query16);
  free(query);
}
