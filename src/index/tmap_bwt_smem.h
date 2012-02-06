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

// TODO
typedef struct {
    tmap_bwt_int_t x[3];
    tmap_bwt_int_t info;
} tmap_bwt_smem_intv_t;

// TODO
typedef struct { 
    size_t n;
    size_t m;
    tmap_bwt_smem_intv_t *a; 
} tmap_bwt_smem_intv_vec_t;

// TODO
#define tmap_bwt_smem_set_intv(bwt, c, ik) ((ik).x[0] = (bwt)->L2[(int32_t)(c)]+1, (ik).x[2] = (bwt)->L2[(int32_t)(c)+1]-(bwt)->L2[(int32_t)(c)], (ik).x[1] = (bwt)->L2[3-(c)]+1, (ik).info = 0)

// TODO
tmap_bwt_int_t
tmap_bwt_smem1(const tmap_bwt_t *bwt, int32_t len, const uint8_t *q, int32_t x, tmap_bwt_smem_intv_vec_t *mem, tmap_bwt_smem_intv_vec_t *tmpvec[2]);

#endif
