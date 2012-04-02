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

#ifndef TMAP_VSW_S0_H
#define TMAP_VSW_S0_H

// TODO: document


#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <unistd.h>
#include <stdio.h>
#include "tmap_vsw_definitions.h"

/*!
  The parameter and memory for the vectorized alignment.
 */
typedef struct {
    uint8_t query_start_clip:2; /*!< 1 if we are to clip the start of the query, 0 otherwise */
    uint8_t query_end_clip:2;  /*!< 1 if we are to clip the end of the query, 0 otherwise */
    uint16_t qlen; /*!< the query length */
    uint16_t tlen; /*!< the target length */
    int32_t slen; /*!< the number of stripes needed */
    int32_t qlen_mem; /*!<  get the amount of memory the stripes will occupy */
    tmap_vsw16_int_t min_edit_score; /*!< the minimum edit score */
    tmap_vsw16_int_t max_edit_score; /*!< the minimum edit score */
    tmap_vsw16_int_t zero_aln_score; /*!< the zero-scoring starting alignment value  */
    tmap_vsw16_int_t min_aln_score; /*!< the minimum possible alignment score */
    tmap_vsw16_int_t max_aln_score; /*!< the maximum global alignment score */
    __m128i *query_profile; /*!< the query scoring profile */
    __m128i *H0; /*!< H/H' */
    __m128i *H1; /*!< H/H' */
    __m128i *E; /*!< deletion from the query, insertion into the target */
    int32_t qlen_max; /*!< the maximum query size */
    tmap_vsw_opt_t *opt; // TODO
    int32_t max_qlen; // TODO
    int32_t max_tlen; // TODO
} tmap_vsw_data_s0_t;

// returns cell j in the scoring matrix row
#define __tmap_vsw16_get_matrix_cell1(_Xs, _j, _slen) \
  (((tmap_vsw16_int_t*)(_Xs + ((_j) % _slen)))[(_j) / _slen])
// returns cell (i,j) in the scoring matrix
#define __tmap_vsw16_get_matrix_cell(_Xs, _slen, _i, _j) \
  __tmap_vsw16_get_matrix_cell1(_Xs + ((_i) * _slen), _j, _slen)
// returns the score in the query profile
#define __tmap_vsw16_get_query_profile_value(_query, _base, _j) \
  (((tmap_vsw16_int_t*)((_query)->query_profile + (_base * (_query)->slen) + ((_j) % (_query)->slen)))[(_j) / (_query)->slen])

/** Auxiliary functions **/
// returns the amount of memory the stripes will occupy  (align the memory by 16 bytes)
#define __tmap_vsw_calc_qlen_mem(_slen) __tmap_vsw_16((_slen) << 4) // divide by 16, since there 16 bytes per 128
// zero-based query index to zero-based stripe #
#define __tmap_vsw_query_index_to_stripe_number(_q, _slen) ((_q) & _slen) 
// zero-based query index to zero-based bit #
#define __tmap_vsw_query_index_to_byte_number(_q, _slen) (((int32_t)(_q)) / _slen)

/*!
  @param  query             the query sequence
  @param  qlen              the query sequence length
  @param  query_start_clip  1 if we are to clip the start of the query, 0 otherwise
  @param  query_end_clip    1 if we are to clip the end of the query, 0 otherwise
  @param  opt               the previous alignment parameters, NULL if none exist
  @return                   the query sequence in vectorized form
  */
tmap_vsw_data_s0_t*
tmap_vsw_data_init_s0(const uint8_t *query, int32_t qlen, int32_t query_start_clip, int32_t query_end_clip, tmap_vsw_opt_t *opt);

// TODO
tmap_vsw_data_s0_t*
tmap_vsw_data_update_s0(tmap_vsw_data_s0_t *vsw, const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen);

/*!
  @param  vsw  the query sequence in vectorized form
 */
void
tmap_vsw_data_destroy_s0(tmap_vsw_data_s0_t *vsw);

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
tmap_vsw_process_s0(tmap_vsw_data_s0_t *vsw,
                    const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen, 
                    int32_t query_start_clip, int32_t query_end_clip, tmap_vsw_opt_t *opt, 
                    int32_t direction, int32_t score_thr, 
                    int32_t *query_end, int32_t *target_end, 
                    int32_t *n_best, int32_t *overflow);
#endif
