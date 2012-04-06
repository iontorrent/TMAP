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

#ifndef TMAP_VSW_S3_H
#define TMAP_VSW_S3_H

// TODO: document


#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <unistd.h>
#include <stdio.h>
#include "tmap_vsw_definitions.h"

// TODO
typedef union { int16_t s[8]; __m128i m; } m128si16;

/*!
  The parameter and memory for the vectorized alignment.
 */
typedef struct {
    int16_t *abuf; //[512];
    int16_t *B; //[1024+16];
    int16_t *MV; //[1024+16];
    m128si16 *X; //[1024+16];
    int32_t mem_qlen;
    int32_t mem_tlen;
    int32_t max_qlen;
    int32_t max_tlen;
    int32_t query_start_clip;
    int32_t query_end_clip;
    tmap_vsw_opt_t *opt; // TODO
} tmap_vsw_data_s3_t;

/*!
  @param  query             the query sequence
  @param  qlen              the query sequence length
  @param  query_start_clip  1 if we are to clip the start of the query, 0 otherwise
  @param  query_end_clip    1 if we are to clip the end of the query, 0 otherwise
  @param  opt               the previous alignment parameters, NULL if none exist
  @return                   the query sequence in vectorized form
  */
tmap_vsw_data_s3_t*
tmap_vsw_data_init_s3(const uint8_t *query, int32_t qlen, int32_t query_start_clip, int32_t query_end_clip, tmap_vsw_opt_t *opt);

// TODO
tmap_vsw_data_s3_t*
tmap_vsw_data_update_s3(tmap_vsw_data_s3_t *vsw, const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen);

/*!
  @param  vsw  the query sequence in vectorized form
 */
void
tmap_vsw_data_destroy_s3(tmap_vsw_data_s3_t *vsw);

/*!
  @param  vsw               TODO
  @param  query             the query sequence
  @param  qlen              the query sequence length
  @param  target            the target sequence
  @param  tlen              the target sequence length
  @param  query_start_clip  1 if we are to clip the start of the query, 0 otherwise
  @param  query_end_clip    1 if we are to clip the end of the query, 0 otherwise
  @param  opt               the alignment parameters
  @param  direction         1 if we are performing forward alignment, 0 otherwise
  @param  score_thr         the minimum scoring threshold (inclusive)
  @param  result            TODO
  @param  overflow           returns 1 if overflow occurs, 0 otherwise
  @return                   the alignment score
  */
int32_t
tmap_vsw_process_s3(tmap_vsw_data_s3_t *vsw,
                    const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen, 
                    int32_t query_start_clip, int32_t query_end_clip, tmap_vsw_opt_t *opt, 
                    int32_t direction, int32_t score_thr, 
                    int32_t *query_end, int32_t *target_end, 
                    int32_t *n_best, int32_t *overflow);
#endif
