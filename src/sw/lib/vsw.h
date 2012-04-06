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

#ifndef VSW_H
#define VSW_H

// TODO: document

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <unistd.h>

//#define VSW_DEBUG 1

// Gives the # of bytes to align the memory in 16-byte increments
#define __vsw_16(_val) (((_val) + 15) >> 4 << 4)

/** Global SIMD macros **/
#define __vsw_mm_store_si128(_dest, _val) _mm_store_si128(_dest, _val)
#define __vsw_mm_load_si128(_src) _mm_load_si128(_src)
#define __vsw_mm_slli_si128(_src, _imm) _mm_slli_si128(_src, _imm)
#define __vsw_mm_srli_si128(_src, _imm) _mm_srli_si128(_src, _imm)

/** VSW MM 16-byte mode (8 values per 128)  **/
#define __vsw16_mm_store_si128(_dest, _val) _mm_store_si128(_dest, _val)
#define __vsw16_mm_load_si128(_src) _mm_load_si128(_src)
#define __vsw16_mm_slli_si128(_src, _imm) _mm_slli_si128(_src, _imm)
#define __vsw16_mm_srli_si128(_src, _imm) _mm_srli_si128(_src, _imm)
#define __vsw16_mm_set1_epi16(_val) _mm_set1_epi16(_val)
#define __vsw16_mm_set_epi16(_r0, _r1, _r2, _r3, _r4, _r5, _r6, _r7) _mm_set_epi16(_r7, _r6, _r5, _r4, _r3, _r2, _r1, _r0)
#define __vsw16_mm_adds_epi16(_a, _b) _mm_adds_epi16(_a, _b)
#define __vsw16_mm_subs_epi16(_a, _b) _mm_subs_epi16(_a, _b)
#define __vsw16_mm_max_epi16(_a, _b) _mm_max_epi16(_a, _b)
#define __vsw16_mm_min_epi16(_a, _b) _mm_min_epi16(_a, _b)
#define __vsw16_mm_movemask_epi16(_a) _mm_movemask_epi8(_a) // NB: this should be 8-bit...
#define __vsw16_mm_cmpeq_epi16(_a, _b) _mm_cmpeq_epi16(_a, _b)
#define __vsw16_mm_cmplt_epi16(_a, _b) _mm_cmplt_epi16(_a, _b)
#define __vsw16_mm_cmpgt_epi16(_a, _b) _mm_cmpgt_epi16(_a, _b)
#define __vsw16_mm_insert_epi16(_a, _b, _imm) _mm_insert_epi16(_a, _b, _imm)
#define __vsw16_mm_extract_epi16(_a, _imm) _mm_extract_epi16(_a, _imm)
#define vsw16_values_per_128_bits 8 
#define vsw16_values_per_128_bits_log2 3 // should be log2(values_per_128_bits)
#define vsw16_max_value INT16_MAX
#define vsw16_mid_value ((vsw16_max_value >> 1)-1)
#define vsw16_min_value INT16_MIN 
#define vsw16_shift_bytes 2 
//typedef uint16_t vsw16_uint_t;
typedef int16_t vsw16_int_t;
// returns the number of stripes required given the query length
#define __vsw16_calc_slen(_qlen) (((_qlen) + vsw16_values_per_128_bits - 1) >> vsw16_values_per_128_bits_log2)
// returns the maximum value in stored in the vector
#define __vsw16_max(ret, xx) do { \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    (ret) = (vsw16_int_t)(__vsw16_mm_extract_epi16((xx), 0) & 0xffff); \
} while (0) 
// returns the minimum value in stored in the vector
#define __vsw16_min(ret, xx) do { \
    (xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 8)); \
    (xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 4)); \
    (xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 2)); \
    (ret) = (vsw16_int_t)(__vsw16_mm_extract_epi16((xx), 0) & 0xffff); \
} while (0) 
// Inserts a value into the given location in the vector.  This overcomes the immediate problem when compiling.
#define __vsw16_mm_insert(_a, _val, _k) \
  (0 == _k ? __vsw16_mm_insert_epi16(_a, _val, 0) : \
   (1 == _k ? __vsw16_mm_insert_epi16(_a, _val, 1) : \
    (2 == _k ? __vsw16_mm_insert_epi16(_a, _val, 2) : \
     (3 == _k ? __vsw16_mm_insert_epi16(_a, _val, 3) : \
      (4 == _k ? __vsw16_mm_insert_epi16(_a, _val, 4) : \
       (5 == _k ? __vsw16_mm_insert_epi16(_a, _val, 5) : \
        (6 == _k ? __vsw16_mm_insert_epi16(_a, _val, 6) : \
         (7 == _k ? __vsw16_mm_insert_epi16(_a, _val, 7) : \
          _a))))))))
// returns cell j in the scoring matrix row
#define __vsw16_get_matrix_cell1(_Xs, _j, _slen) \
  (((vsw16_int_t*)(_Xs + ((_j) % _slen)))[(_j) / _slen])
// returns cell (i,j) in the scoring matrix
#define __vsw16_get_matrix_cell(_Xs, _slen, _i, _j) \
  __vsw16_get_matrix_cell1(_Xs + ((_i) * _slen), _j, _slen)
// returns the score in the query profile
#define __vsw16_get_query_profile_value(_query, _base, _j) \
  (((vsw16_int_t*)((_query)->query_profile + (_base * (_query)->slen) + ((_j) % (_query)->slen)))[(_j) / (_query)->slen])

/** Auxiliary functions **/
// returns the amount of memory the stripes will occupy  (align the memory by 16 bytes)
#define __vsw_calc_qlen_mem(_slen) __vsw_16((_slen) << 4) // divide by 16, since there 16 bytes per 128
// zero-based query index to zero-based stripe #
#define __vsw_query_index_to_stripe_number(_q, _slen) ((_q) & _slen) 
// zero-based query index to zero-based bit #
#define __vsw_query_index_to_byte_number(_q, _slen) (((int32_t)(_q)) / _slen)

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
    vsw16_int_t min_edit_score; /*!< the minimum edit score */
    vsw16_int_t max_edit_score; /*!< the minimum edit score */
    vsw16_int_t zero_aln_score; /*!< the zero-scoring starting alignment value  */
    vsw16_int_t min_aln_score; /*!< the minimum possible alignment score */
    vsw16_int_t max_aln_score; /*!< the maximum global alignment score */
    __m128i *query_profile; /*!< the query scoring profile */
    __m128i *H0; /*!< H/H' */
    __m128i *H1; /*!< H/H' */
    __m128i *E; /*!< deletion from the query, insertion into the target */
    int32_t qlen_max; /*!< the maximum query size */
} vsw16_query_t;

/*!
  Wrapper for the query memory.
 */
typedef struct {
    vsw16_query_t *query16; /*!< the query memory */
    const uint8_t *query; // TMP
    uint32_t qlen; // TMP
} vsw_query_t;

// ACGTN
#define VSW_ALPHABET_SIZE 5

/*!
  The alignment scoring parameters.
  */
typedef struct {
    int32_t score_match; /*!< the match score */
    int32_t pen_mm; /*!< the mismatch penalty */
    int32_t pen_gapo; /*!< the gap open penalty */
    int32_t pen_gape; /*!< the gap extension penalty */
    int32_t score_thres; /*!< the minimum scoring threshold (inclusive) */
} vsw_opt_t;

typedef struct {
    int32_t score_fwd; /*!< the alignment score for the forward algorithm */
    int32_t score_rev; /*!< the alignment score for the reverse algorithm */
    int32_t query_start; /*!< the start alignment position in the query (0-based) */
    int32_t query_end; /*!< the end alignment position in the query (0-based) */
    int32_t target_start; /*!< the start alignment position in the target (0-based) */
    int32_t target_end; /*!< the end alignment position in the target (0-based) */
} vsw_result_t;

/*!
  @param  score_match  the match score
  @param  pen_mm       the mismatch penalty
  @param  pen_gapo     the gap open penalty
  @param  pen_gape     the gap extension penalty
  @param  score_thr    the minimum scoring threshold (inclusive)
  @return              the initialized scoring parameters
 */
vsw_opt_t*
vsw_opt_init(int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape, int32_t score_thr);

/*!
  @param  opt  the parameters to destroy 
 */
void
vsw_opt_destroy(vsw_opt_t *opt);

/*!
  @param  query             the query sequence
  @param  qlen              the query sequence length
  @param  tlen              the target length
  @param  query_start_clip  1 if we are to clip the start of the query, 0 otherwise
  @param  query_end_clip    1 if we are to clip the end of the query, 0 otherwise
  @param  opt               the previous alignment parameters, NULL if none exist
  @return                   the query sequence in vectorized form
  */
vsw_query_t*
vsw_query_init(const uint8_t *query, int32_t qlen,
                    int32_t tlen,
                    int32_t query_start_clip, int32_t query_end_clip,
                    vsw_opt_t *opt);

/*!
  @param  query  the query sequence
 */
void
vsw_query_destroy(vsw_query_t *query);

/*!
  @param  vsw_query         the query in its vectorized form
  @param  query             the query sequence
  @param  qlen              the query sequence length
  @param  target            the target sequence
  @param  tlen              the target sequence length
  @param  query_start_clip  1 if we are to clip the start of the query, 0 otherwise
  @param  query_end_clip    1 if we are to clip the end of the query, 0 otherwise
  @param  opt               the alignment parameters
  @param  score_fwd         the alignment score for the forward smith waterman
  @param  score_rev         the alignment score for the reverse smith waterman
  @param  query_start       the query start position in the alignment (0-based) 
  @param  query_end         the query end position in the alignment (0-based) 
  @param  target_start      the target start position in the alignment (0-based) 
  @param  target_end        the target end position in the alignment (0-based) 
  @param  overflow          returns 1 if overflow occurs, 0 otherwise
  @param  n_best            the number of bset scoring alignments found
  @param  score_thr         the minimum scoring threshold (inclusive)
  @param  is_rev            1 if the reverse alignment is being performed, 0 for the forward
  @return                   the alignment score
  */
int32_t
vsw_sse2(vsw_query_t *vsw_query,
              const uint8_t *query, int32_t qlen,
              uint8_t *target, int32_t tlen,
              int32_t query_start_clip, int32_t query_end_clip,
              vsw_opt_t *opt, 
              int16_t *score_fwd, int16_t *score_rev,
              int16_t *query_start, int16_t *query_end,
              int16_t *target_start, int16_t *target_end,
              int32_t *overflow, int32_t *n_best, int32_t score_thr, int32_t is_rev);

#endif
