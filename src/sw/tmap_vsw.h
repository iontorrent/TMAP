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

#ifndef TMAP_VSW_H
#define TMAP_VSW_H

//#define TMAP_VSW_DEBUG_CMP
//#define TMAP_VSW_DEBUG

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <unistd.h>
#include "tmap_vsw_definitions.h"
#include "lib/AffineSWOptimizationWrapper.h"
  
/*!
  Used to run the underlying vectorized smith waterman (VSW) types
  */
typedef struct {
    int32_t type; /*!< the vectorized smith waterman type */
    tmap_vsw_wrapper_t *algorithm; /*!< the main VSW algorithm */
    tmap_vsw_wrapper_t *algorithm_default; /*!< the VSW algorithm to test if the main algorithm fails */
    int32_t query_start_clip; /*!< 1 if we are to clip the start of the query, 0 otherwise */
    int32_t query_end_clip; /*!< 1 if we are to clip the end of the query, 0 otherwise */
    tmap_vsw_opt_t *opt; /*!< the alignment parameters */
} tmap_vsw_t;

/*!
  @param  query             the query sequence
  @param  qlen              the query sequence length
  @param  query_start_clip  1 if we are to clip the start of the query, 0 otherwise
  @param  query_end_clip    1 if we are to clip the end of the query, 0 otherwise
  @param  type              the VSW type
  @param  opt               the previous alignment parameters, NULL if none exist
  @return                   the query sequence in vectorized form
  */
tmap_vsw_t*
tmap_vsw_init(const uint8_t *query, int32_t qlen,
              int32_t query_start_clip, int32_t query_end_clip,
              int32_t type,
              tmap_vsw_opt_t *opt);

/*!
  @param  vsw  the structure to destroy
  */
void
tmap_vsw_destroy(tmap_vsw_t *vsw);

/*!
  Performs alignment in the sequencing direction.  This will update query_end and target_end
  in the results.
  @param  vsw               the query in its vectorized form
  @param  query             the query sequence
  @param  qlen              the query sequence length
  @param  target            the target sequence
  @param  tlen              the target sequence length
  @param  result            the structure in which to store the results
  @param  overflow          returns 1 if overflow occurs, 0 otherwise
  @param  score_thr         the minimum scoring threshold (inclusive)
  @param  direction         how to break ties
  @return                   the alignment score
  @details direction explains how to break ties. If direction = 0, then query_end needs to be as 
  small as possible.  If direction = 1, then query_end needs to be as large as possible.  In both
  cases, if there are still several possibilities, choose the one with the largest value of 
  target_end among them). 
  */
int32_t
tmap_vsw_process_fwd(tmap_vsw_t *vsw,
                 const uint8_t *query, int32_t qlen,
                 uint8_t *target, int32_t tlen, 
                 tmap_vsw_result_t *result,
                 int32_t *overflow, int32_t score_thr, int32_t direction);

/*!
  Performs alignment in the reverse of the sequencing direction.  This will update query_start 
  and target_start in the results.
  @param  vsw               the query in its vectorized form
  @param  query             the query sequence
  @param  qlen              the query sequence length
  @param  target            the target sequence
  @param  tlen              the target sequence length
  @param  result            the structure in which to store the results
  @param  overflow          returns 1 if overflow occurs, 0 otherwise
  @param  score_thr         the minimum scoring threshold (inclusive)
  @param  direction         how to break ties
  @return                   the alignment score
  @details direction explains how to break ties. If direction = 0, then query_end needs to be as 
  small as possible.  If direction = 1, then query_end needs to be as large as possible.  In both
  cases, if there are still several possibilities, choose the one with the largest value of 
  target_end among them). 
  */
int32_t
tmap_vsw_process_rev(tmap_vsw_t *vsw,
                 const uint8_t *query, int32_t qlen,
                 uint8_t *target, int32_t tlen, 
                 tmap_vsw_result_t *result,
                 int32_t *overflow, int32_t score_thr, int32_t direction);

#endif
