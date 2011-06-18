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

#ifndef TMAP_VSW8_H
#define TMAP_VSW8_H

// TODO: document


#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <unistd.h>
#include <stdio.h>
#include "tmap_vsw.h"

tmap_vsw8_query_t *
tmap_vsw8_query_init_full(tmap_vsw8_query_t *prev, const uint8_t *query, int32_t qlen, int32_t tlen, tmap_vsw_opt_t *opt);

tmap_vsw8_query_t *
tmap_vsw8_query_init_short(tmap_vsw8_query_t *prev, const uint8_t *query, int32_t qlen, tmap_vsw_opt_t *opt);

void
tmap_vsw8_query_reinit(tmap_vsw8_query_t *query, tmap_vsw_result_t *result);

void
tmap_vsw8_query_destroy(tmap_vsw8_query_t *vsw);

void
tmap_vsw8_sse2_forward(tmap_vsw8_query_t *query, const uint8_t *target, int32_t tlen, tmap_vsw_opt_t *opt, tmap_vsw_result_t *result, int32_t *overflow);

void
tmap_vsw8_sse2_reverse(tmap_vsw8_query_t *query, const uint8_t *target, int32_t tlen, tmap_vsw_opt_t *opt, tmap_vsw_result_t *result, int32_t *overflow);

void
tmap_vsw8_sse2(tmap_vsw8_query_t *query, const uint8_t *target, int32_t tlen, tmap_vsw_opt_t *opt, tmap_vsw_result_t *result, int32_t *overflow);

void
tmap_vsw8_sse2_trackback_test(tmap_vsw8_query_t *query, const uint8_t *target, int32_t tlen, tmap_vsw_opt_t *opt, tmap_vsw_result_t *result);

void
tmap_vsw8_sse2_get_path(const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen,
                        tmap_vsw8_query_t *q, tmap_vsw_result_t *result, tmap_sw_path_t *path,
                        int32_t *path_len, tmap_vsw_opt_t *opt);

#endif
