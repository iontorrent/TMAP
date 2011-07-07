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
#include <unistd.h>

//#define TMAP_VSW_DEBUG 1

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
// returns the number of stripes required 
#define __tmap_vsw16_calc_slen(_qlen) (((_qlen) + tmap_vsw16_values_per_128_bits - 1) >> tmap_vsw16_values_per_128_bits_log2)
// ret is max in xx.
#define __tmap_vsw16_max(ret, xx) do { \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    (ret) = (tmap_vsw16_int_t)(__tmap_vsw16_mm_extract_epi16((xx), 0) & 0xffff); \
} while (0) 
// ret is min in xx.
#define __tmap_vsw16_min(ret, xx) do { \
    (xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 8)); \
    (xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 4)); \
    (xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 2)); \
    (ret) = (tmap_vsw16_int_t)(__tmap_vsw16_mm_extract_epi16((xx), 0) & 0xffff); \
} while (0) 
// overcomes the immediate problem
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
// returns cell (i,j) in the scoring matrix
#define __tmap_vsw16_get_matrix_cell1(_Xs, _j, _slen) \
  (((tmap_vsw16_int_t*)(_Xs + ((_j) % _slen)))[(_j) / _slen])
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

typedef struct {
    uint8_t query_start_clip:2;
    uint8_t query_end_clip:2; 
    uint16_t qlen; // the query length
    uint16_t tlen; // the target length
    int32_t slen; // the number of stripes needed
    int32_t qlen_mem; // get the amount of memory the stripes will occupy
    tmap_vsw16_int_t min_edit_score; // the minimum edit score
    tmap_vsw16_int_t max_edit_score; // the minimum edit score
    tmap_vsw16_int_t zero_aln_score; // the zero-scoring starting alignment value 
    tmap_vsw16_int_t min_aln_score; // the minimum possible alignment score
    tmap_vsw16_int_t max_aln_score; // the maximum global alignment score
    __m128i *query_profile; // the query scoring profile
    __m128i *H0; // H/H'
    __m128i *H1; // H/H'
    __m128i *E; // deletion from the query, insertion into the target
    int32_t qlen_max; // the maximum query size
} tmap_vsw16_query_t;

typedef struct {
    tmap_vsw16_query_t *query16_fwd;
    tmap_vsw16_query_t *query16_rev;
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
tmap_vsw_query_init(const uint8_t *query_fwd, int32_t qlen_fwd,
                    const uint8_t *query_rev, int32_t qlen_rev,
                    int32_t tlen,
                    int32_t query_start_clip, int32_t query_end_clip,
                    tmap_vsw_opt_t *opt);

void
tmap_vsw_query_destroy(tmap_vsw_query_t *query);

int32_t
tmap_vsw_sse2(tmap_vsw_query_t *query,
              const uint8_t *query_fwd, int32_t qlen_fwd,
              const uint8_t *query_rev, int32_t qlen_rev,
              uint8_t *target, int32_t tlen,
              int32_t query_start_clip, int32_t query_end_clip,
              tmap_vsw_opt_t *opt, tmap_vsw_result_t *result,
              int32_t *overflow);

#endif
