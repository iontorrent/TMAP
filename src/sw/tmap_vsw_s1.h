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

#ifndef TMAP_VSW_S1_H
#define TMAP_VSW_S1_H

// TODO: document


#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <unistd.h>
#include <stdio.h>
#include "tmap_vsw_definitions.h"
    
typedef struct {
    uint32_t hash;
    short opt;
    short n_best;
    int res_min_pos;
    int res_max_pos;
} tmap_vsw_s1_qres_t;

/*!
  The parameter and memory for the vectorized alignment.
 */
typedef struct {
    tmap_vsw_s1_qres_t *htdata; //htdata[HMSIZE];
    int16_t *buffer __attribute__((aligned(64))); //BUFFER[MAX_DIMB * 9] __attribute__((aligned(64)));
    /*
    int32_t n, m;
    int32_t segNo;
    int32_t len;
    __m128i* XP[4];
    __m128i* M0;
    __m128i* M1;
    __m128i* V;
    int16_t* POS;
    int16_t INVALID_POS[16];
    int32_t INVALID_POS_NO;
    int32_t lastMax;
    */

    int32_t mem_qlen; // TODO
    int32_t mem_tlen; // TODO
    uint8_t query_start_clip:2; /*!< 1 if we are to clip the start of the query, 0 otherwise */
    uint8_t query_end_clip:2;  /*!< 1 if we are to clip the end of the query, 0 otherwise */
    tmap_vsw_opt_t *opt; // TODO
    int32_t max_qlen; // TODO
    int32_t max_tlen; // TODO
} tmap_vsw_data_s1_t;

/*!
  @param  query             the query sequence
  @param  qlen              the query sequence length
  @param  query_start_clip  1 if we are to clip the start of the query, 0 otherwise
  @param  query_end_clip    1 if we are to clip the end of the query, 0 otherwise
  @param  opt               the previous alignment parameters, NULL if none exist
  @return                   the query sequence in vectorized form
  */
tmap_vsw_data_s1_t*
tmap_vsw_data_init_s1(const uint8_t *query, int32_t qlen, int32_t query_start_clip, int32_t query_end_clip, tmap_vsw_opt_t *opt);

// TODO
tmap_vsw_data_s1_t*
tmap_vsw_data_update_s1(tmap_vsw_data_s1_t *vsw, const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen);

/*!
  @param  vsw  the query sequence in vectorized form
 */
void
tmap_vsw_data_destroy_s1(tmap_vsw_data_s1_t *vsw);

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
tmap_vsw_process_s1(tmap_vsw_data_s1_t *vsw,
                    const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen, 
                    int32_t query_start_clip, int32_t query_end_clip, tmap_vsw_opt_t *opt, 
                    int32_t direction, int32_t score_thr, 
                    int32_t *query_end, int32_t *target_end, 
                    int32_t *n_best, int32_t *overflow);
#endif
