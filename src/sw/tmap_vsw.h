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

// TODO: document

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <unistd.h>

/* For branch prediction */
#ifdef __GNUC__
#define TMAP_VSW_LIKELY(x) __builtin_expect((x),1)
#define TMAP_VSW_UNLIKELY(x) __builtin_expect((x),0)
#else
#define TMAP_VSW_LIKELY(x) (x)
#define TMAP_VSW_UNLIKELY(x) (x)
#endif

// Gives the # of bytes to align the memory in 16-byte increments
#define __tmap_vsw_16(_val) (((_val) + 15) >> 4 << 4)

/** Global SIMD macros **/
#define __tmap_vsw_mm_store_si128(_dest, _val) _mm_store_si128(_dest, _val)
#define __tmap_vsw_mm_load_si128(_src) _mm_load_si128(_src)
#define __tmap_vsw_mm_slli_si128(_src, _imm) _mm_slli_si128(_src, _imm)
#define __tmap_vsw_mm_srli_si128(_src, _imm) _mm_srli_si128(_src, _imm)

/** VSW MM 8-byte mode (16 values per 128)  **/
#define __tmap_vsw8_mm_set_epi(_val) _mm_set1_epi8(_val)
#define __tmap_vsw8_mm_adds_epu(_a, _b) _mm_adds_epu8(_a, _b)
#define __tmap_vsw8_mm_subs_epu(_a, _b) _mm_subs_epu8(_a, _b)
#define __tmap_vsw8_mm_max_epu(_a, _b) _mm_max_epu8(_a, _b)
#define __tmap_vsw8_mm_movemask_epi(_a) _mm_movemask_epi8(_a)
#define tmap_vsw8_values_per_128_bits 16 
#define tmap_vsw8_values_per_128_bits_log2 4 // should be log2(values_per_128_bits)
#define tmap_vsw8_max_value 256 
#define tmap_vsw8_mid_value ((tmap_vsw8_max_value >> 1)-1)
#define tmap_vsw8_shift_bits 1
typedef uint8_t tmap_vsw8_uint_t;
typedef int8_t tmap_vsw8_int_t;
// returns the number of stripes required 
#define __tmap_vsw8_calc_slen(_qlen) (((_qlen) + tmap_vsw8_values_per_128_bits - 1) >> tmap_vsw8_values_per_128_bits_log2)
// ret is max in xx.
#define __tmap_vsw8_max(ret, xx) do { \
    (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
    (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
    (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
    (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
    (ret) = _mm_extract_epi16((xx), 0) & 0x00ff; \
} while (0)
// returns cell (i,j) in the scoring matrix
#define __tmap_vsw8_get_matrix_cell(_query, _i, _j) \
  (((tmap_vsw8_uint_t*)((_query)->H + ((_i) * ((_query)->qlen_mem >> 4)) + ((_j) % (_query)->slen)))[(_j) / (_query)->slen])
// returns the score in the query profile
#define __tmap_vsw8_get_query_profile_value(_query, _base, _j) \
  (((tmap_vsw8_uint_t*)((_query)->query_profile + (_base * (_query)->slen) + ((_j) % (_query)->slen)))[(_j) / (_query)->slen])

/** VSW MM 16-byte mode (8 values per 128)  **/
#define __tmap_vsw16_mm_store_si128(_dest, _val) _mm_store_si128(_dest, _val)
#define __tmap_vsw16_mm_load_si128(_src) _mm_load_si128(_src)
#define __tmap_vsw16_mm_slli_si128(_src, _imm) _mm_slli_si128(_src, _imm)
#define __tmap_vsw16_mm_srli_si128(_src, _imm) _mm_srli_si128(_src, _imm)
#define __tmap_vsw16_mm_set_epi(_val) _mm_set1_epi16(_val)
#define __tmap_vsw16_mm_adds_epu(_a, _b) _mm_adds_epu16(_a, _b)
#define __tmap_vsw16_mm_subs_epu(_a, _b) _mm_subs_epu16(_a, _b)
#define __tmap_vsw16_mm_max_epu(_a, _b) _mm_max_epu16(_a, _b)
#define __tmap_vsw16_mm_movemask_epi(_a) _mm_movemask_epi8(_a)
#define tmap_vsw16_values_per_128_bits 8 
#define tmap_vsw16_values_per_128_bits_log2 3 // should be log2(values_per_128_bits)
#define tmap_vsw16_max_value 65536 
#define tmap_vsw16_mid_value ((tmap_vsw16_max_value >> 1)-1)
#define tmap_vsw16_shift_bits 2
typedef uint16_t tmap_vsw16_uint_t;
typedef int16_t tmap_vsw16_int_t;
// returns the number of stripes required 
#define __tmap_vsw16_calc_slen(_qlen) (((_qlen) + tmap_vsw16_values_per_128_bits - 1) >> tmap_vsw16_values_per_128_bits_log2)
// ret is max in xx.
#define __tmap_vsw16_max(ret, xx) do { \
    (xx) = _mm_max_epu16((xx), _mm_srli_si128((xx), 8)); \
    (xx) = _mm_max_epu16((xx), _mm_srli_si128((xx), 4)); \
    (xx) = _mm_max_epu16((xx), _mm_srli_si128((xx), 2)); \
    (ret) = _mm_extract_epi16((xx), 0) & 0xffff; \
} while (0)
// returns cell (i,j) in the scoring matrix
#define __tmap_vsw16_get_matrix_cell(_query, _i, _j) \
  (((tmap_vsw16_uint_t*)((_query)->H + ((_i) * ((_query)->qlen_mem >> 4)) + ((_j) % (_query)->slen)))[(_j) / (_query)->slen])
// returns the score in the query profile
#define __tmap_vsw16_get_query_profile_value(_query, _base, _j) \
  (((tmap_vsw16_uint_t*)((_query)->query_profile + (_base * (_query)->slen) + ((_j) % (_query)->slen)))[(_j) / (_query)->slen])

/** Auxiliary functions **/
// returns the amount of memory the stripes will occupy  (align the memory by 16 bytes)
#define __tmap_vsw_calc_qlen_mem(_slen) __tmap_vsw_16((_slen) << 4) // divide by 16, since there 16 bytes per 128
// zero-based query index to zero-based stripe #
#define __tmap_vsw_query_index_to_stripe_number(_q, _slen) ((_q) & _slen) 
// zero-based query index to zero-based bit #
#define __tmap_vsw_query_index_to_byte_number(_q, _slen) (((int32_t)(_q)) / _slen)

typedef struct {
    int32_t type; // 0 = score only, 1 = full  
    int32_t qlen; // the query length
    int32_t tlen; // the target length
    int32_t slen; // the number of stripes needed
    int32_t qlen_mem; // get the amount of memory the stripes will occupy
    tmap_vsw8_uint_t min_score; // the minimum score
    tmap_vsw8_uint_t range_score; // the difference between the minimum and maximum score
    __m128i *query_profile; // the query scoring profile
    /* -- For the score only VSW -- */
    __m128i *H0; // H/H'
    __m128i *H1; // H/H'
    __m128i *E; // deletion from the query, insert into the target
    int32_t qlen_max; // the maximum query size
    /* -- End score VSW -- */
    /* -- For the full VSW --*/
    int32_t hlen; // the # of cells in the matrix (type == 1)
    int32_t hlen_mem; // the memory allocated for the matrix
    int32_t hlen_max; // the maximum matrix size
    __m128i *H; // H
    /* -- End full VSW -- */
} tmap_vsw8_query_t;

typedef struct {
    int32_t type; // 0 = score only, 1 = full  
    int32_t qlen; // the query length
    int32_t tlen; // the target length
    int32_t slen; // the number of stripes needed
    int32_t qlen_mem; // get the amount of memory the stripes will occupy
    tmap_vsw16_uint_t min_score; // the minimum score
    tmap_vsw16_uint_t range_score; // the difference between the minimum and maximum score
    __m128i *query_profile; // the query scoring profile
    /* -- For the score only VSW -- */
    __m128i *H0; // H/H'
    __m128i *H1; // H/H'
    __m128i *E; // deletion from the query, insert into the target
    int32_t qlen_max; // the maximum query size
    /* -- End score VSW -- */
    /* -- For the full VSW --*/
    int32_t hlen; // the # of cells in the matrix (type == 1)
    int32_t hlen_mem; // the memory allocated for the matrix
    int32_t hlen_max; // the maximum matrix size
    __m128i *H; // H
    /* -- End full VSW -- */
} tmap_vsw16_query_t;

typedef struct {
    tmap_vsw8_query_t *query8;
    tmap_vsw16_query_t *query16;
} tmap_vsw_query_t;

// ACGTN
#define TMAP_VSW_ALPHABET_SIZE 5

typedef struct {
    // scoring parameters
    int32_t score_match;
    int32_t pen_mm;
    uint32_t pen_gapo;
    uint32_t pen_gape;
    int32_t score_thres;
} tmap_vsw_opt_t;

typedef struct {
    int32_t byte_type;
    int32_t score_fwd;
    int32_t score_rev;
    int32_t query_start;
    int32_t query_end;
    int32_t target_start;
    int32_t target_end;
} tmap_vsw_result_t;

tmap_vsw_opt_t*
tmap_vsw_opt_init(int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape, int32_t score_thr);

void
tmap_vsw_opt_destroy(tmap_vsw_opt_t *opt);

tmap_vsw_result_t *
tmap_vsw_result_init();

void
tmap_vsw_result_destroy(tmap_vsw_result_t *result);

tmap_vsw_query_t*
tmap_vsw_query_init();

void
tmap_vsw_query_destroy(tmap_vsw_query_t *query);

int32_t
tmap_vsw_sse2(tmap_vsw_query_t *q,
              const uint8_t *query, int32_t qlen,
              const uint8_t *target, int32_t tlen,
              tmap_vsw_opt_t *opt, tmap_vsw_result_t *result);

void
tmap_vsw_sse2_get_path(const uint8_t *query, int32_t qlen,
                       const uint8_t *target, int32_t tlen,
                       tmap_vsw_query_t *q,
                       tmap_vsw_result_t *result,
                       tmap_sw_path_t *path,
                       int32_t *path_len,
                       tmap_vsw_opt_t *opt);

#endif
