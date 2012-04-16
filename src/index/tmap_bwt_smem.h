/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_BWT_SMEM_H
#define TMAP_BWT_SMEM_H

/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

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

#include <stdint.h>
#include "../util/tmap_definitions.h"

/*!
  Bi-directional occurrence intervals
 */
typedef struct {
    tmap_bwt_int_t x[2]; /*!< the forward and backward suffix intervals */
    tmap_bwt_int_t size; /*!< the interval size */
    uint64_t info; /*!< lower 32 bits store the query index for the forward search, the upper 32 bits store the query index for the reverse search */
    int8_t flag; /*!< 1 if this is a representative repetitive hit, 0 otherwise */
} tmap_bwt_smem_intv_t;

/*!
  A vector of intervals
 */
typedef struct { 
    size_t n; /*!< the number of intervals */
    size_t m; /*!< the amount of memory allocated for the intervals */
    tmap_bwt_smem_intv_t *a; /*!< the intervals */
} tmap_bwt_smem_intv_vec_t;

/*!
  @param  bwt  the bwt
  @param  c    the query base
  @param  ik   the bi-directional occurrence interval
  @details  resets the interval structure based on the given query base
 */
#define tmap_bwt_smem_set_intv(bwt, c, ik) ((ik).x[0] = (bwt)->L2[(int32_t)(c)]+1, (ik).size = (bwt)->L2[(int32_t)(c)+1]-(bwt)->L2[(int32_t)(c)], (ik).x[1] = (bwt)->L2[3-(c)]+1, (ik).info = 0)

/*!
  @param  bwt     the bwt
  @param  len     the query length
  @param  q       the query
  @param  x       the start index into the query for the search (0-based)
  @param  mem     where the matches should be stored
  @param  tmpvec  temporary memory
  @return         the new index into the query (0-based)
 */
int32_t
tmap_bwt_smem1(const tmap_bwt_t *bwt, int32_t len, const uint8_t *q, int32_t x, tmap_bwt_smem_intv_vec_t *mem, tmap_bwt_smem_intv_vec_t *tmpvec[2]);

#endif
