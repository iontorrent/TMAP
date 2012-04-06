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

#ifndef TMAP_VSW_DEFINITIONS_H
#define TMAP_VSW_DEFINITIONS_H

// TODO: document

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <unistd.h>

// Gives the # of bytes to align the memory in 16-byte increments
#define __tmap_vsw_16(_val) (((_val) + 15) >> 4 << 4)

/** Global SIMD macros **/
#define __tmap_vsw_mm_store_si128(_dest, _val) _mm_store_si128(_dest, _val)
#define __tmap_vsw_mm_load_si128(_src) _mm_load_si128(_src)
#define __tmap_vsw_mm_slli_si128(_src, _imm) _mm_slli_si128(_src, _imm)
#define __tmap_vsw_mm_srli_si128(_src, _imm) _mm_srli_si128(_src, _imm)

/** VSW MM 16-byte mode (8 values per 128)  **/
#define __tmap_vsw16_mm_store_si128(_dest, _val) _mm_store_si128(_dest, _val)
#define __tmap_vsw16_mm_load_si128(_src) _mm_load_si128(_src)
#define __tmap_vsw16_mm_slli_si128(_src, _imm) _mm_slli_si128(_src, _imm)
#define __tmap_vsw16_mm_srli_si128(_src, _imm) _mm_srli_si128(_src, _imm)
#define __tmap_vsw16_mm_set1_epi16(_val) _mm_set1_epi16(_val)
#define __tmap_vsw16_mm_set_epi16(_r0, _r1, _r2, _r3, _r4, _r5, _r6, _r7) _mm_set_epi16(_r7, _r6, _r5, _r4, _r3, _r2, _r1, _r0)
#define __tmap_vsw16_mm_adds_epi16(_a, _b) _mm_adds_epi16(_a, _b)
#define __tmap_vsw16_mm_subs_epi16(_a, _b) _mm_subs_epi16(_a, _b)
#define __tmap_vsw16_mm_max_epi16(_a, _b) _mm_max_epi16(_a, _b)
#define __tmap_vsw16_mm_min_epi16(_a, _b) _mm_min_epi16(_a, _b)
#define __tmap_vsw16_mm_movemask_epi16(_a) _mm_movemask_epi8(_a) // NB: this should be 8-bit...
#define __tmap_vsw16_mm_cmpeq_epi16(_a, _b) _mm_cmpeq_epi16(_a, _b)
#define __tmap_vsw16_mm_cmplt_epi16(_a, _b) _mm_cmplt_epi16(_a, _b)
#define __tmap_vsw16_mm_cmpgt_epi16(_a, _b) _mm_cmpgt_epi16(_a, _b)
#define __tmap_vsw16_mm_insert_epi16(_a, _b, _imm) _mm_insert_epi16(_a, _b, _imm)
#define __tmap_vsw16_mm_extract_epi16(_a, _imm) _mm_extract_epi16(_a, _imm)
#define tmap_vsw16_values_per_128_bits 8 
#define tmap_vsw16_values_per_128_bits_log2 3 // should be log2(values_per_128_bits)
#define tmap_vsw16_max_value INT16_MAX
#define tmap_vsw16_mid_value ((tmap_vsw16_max_value >> 1)-1)
#define tmap_vsw16_min_value INT16_MIN 
#define tmap_vsw16_shift_bytes 2 
//typedef uint16_t tmap_vsw16_uint_t;
typedef int16_t tmap_vsw16_int_t;
// returns the number of stripes required given the query length
#define __tmap_vsw16_calc_slen(_qlen) (((_qlen) + tmap_vsw16_values_per_128_bits - 1) >> tmap_vsw16_values_per_128_bits_log2)
// returns the maximum value in stored in the vector
#define __tmap_vsw16_max(ret, xx) do { \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    (ret) = (tmap_vsw16_int_t)(__tmap_vsw16_mm_extract_epi16((xx), 0) & 0xffff); \
} while (0) 
// returns the minimum value in stored in the vector
#define __tmap_vsw16_min(ret, xx) do { \
    (xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 8)); \
    (xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 4)); \
    (xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 2)); \
    (ret) = (tmap_vsw16_int_t)(__tmap_vsw16_mm_extract_epi16((xx), 0) & 0xffff); \
} while (0) 
// Inserts a value into the given location in the vector.  This overcomes the immediate problem when compiling.
#define __tmap_vsw16_mm_insert(_a, _val, _k) \
  (0 == _k ? __tmap_vsw16_mm_insert_epi16(_a, _val, 0) : \
   (1 == _k ? __tmap_vsw16_mm_insert_epi16(_a, _val, 1) : \
    (2 == _k ? __tmap_vsw16_mm_insert_epi16(_a, _val, 2) : \
     (3 == _k ? __tmap_vsw16_mm_insert_epi16(_a, _val, 3) : \
      (4 == _k ? __tmap_vsw16_mm_insert_epi16(_a, _val, 4) : \
       (5 == _k ? __tmap_vsw16_mm_insert_epi16(_a, _val, 5) : \
        (6 == _k ? __tmap_vsw16_mm_insert_epi16(_a, _val, 6) : \
         (7 == _k ? __tmap_vsw16_mm_insert_epi16(_a, _val, 7) : \
          _a))))))))
// ACGTN
#define TMAP_VSW_ALPHABET_SIZE 5

/*!
  The alignment scoring parameters.
  */
typedef struct {
    int32_t score_match; /*!< the match score */
    int32_t pen_mm; /*!< the mismatch penalty */
    int32_t pen_gapo; /*!< the gap open penalty */
    int32_t pen_gape; /*!< the gap extension penalty */
    int32_t score_thres; /*!< the minimum scoring threshold (inclusive) */
} tmap_vsw_opt_t;

typedef struct {
    int32_t score_fwd; /*!< the alignment score for the forward algorithm */
    int32_t score_rev; /*!< the alignment score for the reverse algorithm */
    int32_t query_start; /*!< the start alignment position in the query (0-based) */
    int32_t query_end; /*!< the end alignment position in the query (0-based) */
    int32_t target_start; /*!< the start alignment position in the target (0-based) */
    int32_t target_end; /*!< the end alignment position in the target (0-based) */
    int32_t n_best; // TODO
} tmap_vsw_result_t;

/*!
  @param  score_match  the match score
  @param  pen_mm       the mismatch penalty
  @param  pen_gapo     the gap open penalty
  @param  pen_gape     the gap extension penalty
  @param  score_thr    the minimum scoring threshold (inclusive)
  @return              the initialized scoring parameters
 */
tmap_vsw_opt_t*
tmap_vsw_opt_init(int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape, int32_t score_thr);

/*!
  @param  opt  the parameters to destroy 
 */
void
tmap_vsw_opt_destroy(tmap_vsw_opt_t *opt);

#endif
