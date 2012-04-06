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

#ifndef VSW16_H
#define VSW16_H


#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <unistd.h>
#include <stdio.h>
#include "vsw.h"

/* For branch prediction */
#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif


/*!
  @param  prev              the previous vectorized query sequence, NULL otherwise
  @param  query             the query sequence
  @param  qlen              the query sequence length
  @param  query_start_clip  1 if we are to clip the start of the query, 0 otherwise
  @param  query_end_clip    1 if we are to clip the end of the query, 0 otherwise
  @param  opt               the previous alignment parameters, NULL if none exist
  @return                   the query sequence in vectorized form
  */
vsw16_query_t *
vsw16_query_init(vsw16_query_t *prev, const uint8_t *query, int32_t qlen, 
                              int32_t query_start_clip, int32_t query_end_clip,
                              vsw_opt_t *opt);

/*!
  @param  vsw  the query sequence in vectorized form
 */
void
vsw16_query_destroy(vsw16_query_t *vsw);

/*!
  @param  query             the query sequence
  @param  target            the target sequence
  @param  tlen              the target sequence length
  @param  query_start_clip  1 if we are to clip the start of the query, 0 otherwise
  @param  query_end_clip    1 if we are to clip the end of the query, 0 otherwise
  @param  opt               the alignment parameters
  @param  query_end         the query end position in the alignment (0-based) 
  @param  target_end        the target end position in the alignment (0-based) 
  @param  direction         1 if we are performing forward alignment, 0 otherwise
  @param  overflow          returns 1 if overflow occurs, 0 otherwise
  @param  n_best            the number of bset scoring alignments found
  @param  score_thr         the minimum scoring threshold (inclusive)
  @return                   the alignment score
  */
int32_t
vsw16_sse2_forward(vsw16_query_t *query, const uint8_t *target, int32_t tlen,
                        int32_t query_start_clip, int32_t query_end_clip,
                        vsw_opt_t *opt, int16_t *query_end, int16_t *target_end,
                        int32_t direction, int32_t *overflow, int32_t *n_best, int32_t score_thr);
#endif
