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

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <unistd.h>
#include <stdio.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "tmap_sw.h"
#include "tmap_vsw.h"
#include "tmap_vsw8.h"

static tmap_vsw8_query_t *
tmap_vsw8_query_init_by_type(tmap_vsw8_query_t *prev, const uint8_t *query, int32_t qlen, int32_t tlen, int32_t type, tmap_vsw_opt_t *opt)
{
  int32_t qlen_mem, slen, a;
  tmap_vsw8_uint_t *t;
  int32_t hlen = 0, hlen_mem = 0;

  // get the number of stripes
  slen = __tmap_vsw8_calc_slen(qlen); 
  if(1 == type) {
      hlen = qlen * tlen;
      hlen_mem = __tmap_vsw_calc_qlen_mem(hlen); // same calculation
      if(NULL != prev && hlen == prev->hlen) return prev; // no need to carry on
  }
  else if(NULL != prev && qlen == prev->qlen) return prev; // no need to carry on

  if(NULL == prev 
     || (0 == type && prev->qlen_max < qlen) 
     || (1 == type && prev->hlen_mem < hlen_mem) 
     || type != prev->type) { // recompute
      // free previous
      if(NULL != prev) {
          free(prev); prev = NULL;
      }
      // get the memory needed to hold the stripes for the query 
      qlen_mem = __tmap_vsw_calc_qlen_mem(slen);
      if(0 == type) {
          prev = tmap_malloc(sizeof(tmap_vsw8_query_t) + 15 + qlen_mem * (TMAP_VSW_ALPHABET_SIZE + 3), "prev"); // add three for H0, H1, and E
          // update the type
          prev->type = 0;
          prev->qlen_max = qlen; 
      }
      else if(1 == type) {
          // allocate a single block of memory (add 15 for 16-byte aligning the prev
          // pointer): struct + (16-byte alignment) + query profile + scoring matrix 
          prev = tmap_malloc(sizeof(tmap_vsw8_query_t) + 15 + (qlen_mem * TMAP_VSW_ALPHABET_SIZE + 3) + (hlen_mem), "prev"); 
          // update the type
          prev->type = 1;
          prev->hlen_mem = hlen_mem;
          prev->hlen_max = hlen; 
      }
      else {
          tmap_error("bug encountered", Exit, OutOfRange);
      }
      // update the memory
      prev->qlen_mem = qlen_mem;
      // compute the min/max score and range
      prev->min_score = (opt->pen_gapo+opt->pen_gape < opt->pen_mm) ? -opt->pen_mm : -(opt->pen_gapo+opt->pen_gape); // minimum score
      prev->min_score = tmap_vsw8_max_value - prev->min_score; // NB: prev->min_score is tmap_vsw8_uint_t
      prev->range_score = tmap_vsw8_max_value - opt->score_match;
      prev->range_score += prev->min_score; // this is the difference between the min and max scores
  }

  prev->qlen = qlen; // update the query length
  prev->hlen = hlen; // the scoring matrix size 
  prev->slen = slen; // update the number of stripes

  // NB: align all the memory from one block
  prev->query_profile = (__m128i*)__tmap_vsw_16((size_t)prev + sizeof(tmap_vsw8_query_t)); // skip over the struct variables, align memory 
  prev->H0 = prev->query_profile + (prev->slen * TMAP_VSW_ALPHABET_SIZE); // skip over the query profile
  prev->H1 = prev->H0 + prev->slen; // skip over H0
  prev->E = prev->H1 + prev->slen; // skip over H1
  if(1 == type) {
      prev->H = prev->E + prev->slen; // skip over E
  }
  
  // create the query profile
  t = (tmap_vsw8_uint_t*)prev->query_profile;
  tmap_vsw8_uint_t u = (__tmap_vsw_16((size_t)prev + sizeof(tmap_vsw8_query_t))-(size_t)prev);
  for(a = 0; a < TMAP_VSW_ALPHABET_SIZE; a++) {
      int32_t i, k;
      for(i = 0; i < prev->slen; ++i) { // for each stripe
          // fill in this stripe
          for(k = i; k < prev->slen << tmap_vsw8_values_per_128_bits_log2; k += prev->slen) { //  for q_{i+1}, q_{2i+1}, ..., q_{2s+1}
              // NB: pad with zeros
              //fprintf(stderr, "On byte %lu\n", u);
              *t++ = ((k >= qlen) ? 0 : ((a == query[k]) ? opt->score_match : -opt->pen_mm)) + prev->min_score;
              u += sizeof(tmap_vsw8_uint_t);
          }
      }
  }
  return prev;
}

tmap_vsw8_query_t *
tmap_vsw8_query_init_full(tmap_vsw8_query_t *prev, const uint8_t *query, int32_t qlen, int32_t tlen, tmap_vsw_opt_t *opt)
{
  return tmap_vsw8_query_init_by_type(prev, query, qlen, tlen, 1, opt);
}

tmap_vsw8_query_t *
tmap_vsw8_query_init_short(tmap_vsw8_query_t *prev, const uint8_t *query, int32_t qlen, tmap_vsw_opt_t *opt)
{
  return tmap_vsw8_query_init_by_type(prev, query, qlen, 0, 0, opt);
}
  
/*
static void
tmap_vsw8_query_print_query_profile(tmap_vsw8_query_t *query)
{
  int32_t a, i, j;

  for(a=0;a<TMAP_VSW_ALPHABET_SIZE;a++) {
      fprintf(stderr, "a=%d ", a);
      __m128i *S = query->query_profile + a * query->slen;
      for(i = 0; i < tmap_vsw8_values_per_128_bits; i++) {
          for(j = 0; j < query->slen; j++) {
              fprintf(stderr, " %d", ((tmap_vsw8_uint_t*)(S+j))[i] - query->min_score);
          }
      }
      fprintf(stderr, "\n");
  }
}
*/

void
tmap_vsw8_query_reinit(tmap_vsw8_query_t *query, tmap_vsw_result_t *result)
{
  tmap_vsw8_int_t *new, *old;
  int32_t a, i, k;
  int32_t slen;

  if(0 <= result->query_end && result->query_end+1 < query->qlen) {
      // update the query profile
      slen = __tmap_vsw8_calc_slen(query->qlen);
      new = old = (tmap_vsw8_int_t*)query->query_profile;
      for(a = 0; a < TMAP_VSW_ALPHABET_SIZE; a++) {
          for(i = 0; i < slen; ++i) { // for each new stripe
              for(k = i; k < slen << tmap_vsw8_values_per_128_bits_log2; k += slen) { //  for q_{i+1}, q_{2i+1}, ..., q_{2s+1}
                  // NB: pad with zeros
                  if(query->qlen <= k) {
                      *new = query->min_score;
                  }
                  else {
                      // use the old query profile
                      *new = *old;
                  }
                  new++;
                  old++;
              }
          }
          // skip over old stripes
          for(i = slen; i < query->slen; i++) { // for each old stripe (to skip)
              for(k = i; k < query->slen << tmap_vsw8_values_per_128_bits_log2; k += query->slen) { //  for q_{i+1}, q_{2i+1}, ..., q_{2s+1}
                  old++;
              }
          }
      }
      // update variables
      query->qlen = result->query_end+1;
      query->slen = slen;
      // update H0, H1, and E pointers
      query->H0 = query->query_profile + (query->slen * TMAP_VSW_ALPHABET_SIZE); // skip over the query profile
      query->H1 = query->H0 + query->slen; // skip over H0
      query->E = query->H1 + query->slen; // skip over H1
      if(1 == query->type) {
          // update H
          query->H = query->E + query->slen; // skip over E
      }
  }
}

void
tmap_vsw8_query_destroy(tmap_vsw8_query_t *vsw)
{
  // all memory was allocated as one block
  free(vsw);
}

// H -> diagonal
// E -> deletion (query)
// F -> insertion (query)
// H' -> max{H, E, F}
void
tmap_vsw8_sse2_forward(tmap_vsw8_query_t *query, const uint8_t *target, int32_t tlen, tmap_vsw_opt_t *opt, tmap_vsw_result_t *result, int32_t *overflow)
{
  int32_t slen, i, sum, m_b, n_b, _thres = tmap_vsw8_max_value - query->range_score - 2;
  int32_t gmax = 0;
  __m128i zero, pen_gapoe, pen_gape, min_score, reduce, thres, *H0, *H1, *E, *H;
      
  // initialization
  zero = _mm_set1_epi32(0); // mmz zeros
  m_b = n_b = 0; 
  pen_gapoe = __tmap_vsw8_mm_set_epi(opt->pen_gapo + opt->pen_gape); // gap open penalty
  pen_gape = __tmap_vsw8_mm_set_epi(opt->pen_gape); // gap extend penalty
  min_score = __tmap_vsw8_mm_set_epi(query->min_score); // the minimum score
  thres = __tmap_vsw8_mm_set_epi(_thres); // when max score exceeds this, min_score all scores by "reduce" below
  reduce = __tmap_vsw8_mm_set_epi(tmap_vsw8_mid_value);
  H0 = query->H0; 
  H1 = query->H1; 
  E = query->E; 
  H = query->H;
  slen = query->slen;
  if(NULL != overflow) *overflow = 0;

  // HERE
  //tmap_vsw8_query_print_query_profile(query);
  // select query end only
  // check stripe #: __tmap_vsw8_query_index_to_stripe_number(query->qlen, slen) 
  // check byte # __tmap_vsw8_query_index_to_byte_number(query->qlen, slen)

  // initialize all zeros
  for(i = 0; i < slen; ++i) {
      __tmap_vsw_mm_store_si128(E + i, zero);
      __tmap_vsw_mm_store_si128(H0 + i, zero);
  }

  // the core loop
  for(i = 0, sum = 0; i < tlen; i++) { // for each base in the target
      int j, k, cmp, imax = 0;
      // NB: f and max are initially zero
      __m128i e, h, f = zero, max = zero, *S = query->query_profile + target[i] * slen; // s is the 1st score vector
      h = __tmap_vsw_mm_load_si128(H0 + slen - 1); // the last stripe
      h = __tmap_vsw_mm_slli_si128(h, tmap_vsw8_shift_bits); // h=H(i-1,-1); << instead of >> because x64 is little-endian
      for(j = 0; TMAP_VSW_LIKELY(j < slen); ++j) { // for each stripe in the query
          /* SW cells are computed in the following order:
           *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
           *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
           *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
           */
          // compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
          h = __tmap_vsw8_mm_adds_epu(h, __tmap_vsw_mm_load_si128(S + j));
          h = __tmap_vsw8_mm_subs_epu(h, min_score); // h=H'(i-1,j-1)+S(i,j)
          e = __tmap_vsw_mm_load_si128(E + j); // e=E'(i,j)
          h = __tmap_vsw8_mm_max_epu(h, e); // h=H'(i,j) = max{E'(i,j), H'(i-1,j-1)+S(i,j)}
          h = __tmap_vsw8_mm_max_epu(h, f); // h=H'(i,j) = max{max{E'(i,j), H'(i-1,j-1)+S(i,j)}, F}
          max = __tmap_vsw8_mm_max_epu(max, h); // save the max values in this stripe versus the last
          __tmap_vsw_mm_store_si128(H1 + j, h); // save h to H'(i,j)
          // next, compute E'(i+1,j)
          h = __tmap_vsw8_mm_subs_epu(h, pen_gapoe); // h=H'(i,j)-pen_gapo
          e = __tmap_vsw8_mm_subs_epu(e, pen_gape); // e=E'(i,j)-pen_gape
          e = __tmap_vsw8_mm_max_epu(e, h); // e=E'(i+1,j) = max{E'(i,j)-pen_gape, H'(i,j)-pen_gapo}
          __tmap_vsw_mm_store_si128(E + j, e); // save e to E'(i+1,j)
          // now compute F'(i,j+1)
          //h = __tmap_vsw8_mm_subs_epu(h, pen_gapoe); // h=H'(i,j)-pen_gapo
          f = __tmap_vsw8_mm_subs_epu(f, pen_gape); // f=F'(i,j)-pen_gape
          f = __tmap_vsw8_mm_max_epu(f, h); // f=F'(i,j+1) = max{F'(i,j)-pen_gape, H'(i,j)-pen_gapo}
          // get H'(i-1,j) and prepare for the next j
          h = __tmap_vsw_mm_load_si128(H0 + j); // h=H'(i-1,j)
      }
      // NB: we do not need to set E(i,j) as we disallow adjacent insertion and then deletion
      for(k = 0; TMAP_VSW_LIKELY(k < tmap_vsw8_values_per_128_bits); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
          f = __tmap_vsw_mm_slli_si128(f, tmap_vsw8_shift_bits);
          for(j = 0; TMAP_VSW_LIKELY(j < slen); ++j) {
              h = __tmap_vsw_mm_load_si128(H1 + j);
              h = __tmap_vsw8_mm_max_epu(h, f); // h=H'(i,j)
              __tmap_vsw_mm_store_si128(H1 + j, h);
              h = __tmap_vsw8_mm_subs_epu(h, pen_gapoe);
              f = __tmap_vsw8_mm_subs_epu(f, pen_gape);
              cmp = __tmap_vsw8_mm_movemask_epi(_mm_cmpeq_epi8(__tmap_vsw8_mm_subs_epu(f, h), zero));
              if (TMAP_VSW_UNLIKELY(cmp == 0xffff)) goto end_loop; // ACK: goto statement
          }
      }
end_loop:
      if(1 == query->type) { 
          // save the matrix
          for(j = 0; TMAP_VSW_LIKELY(j < slen); ++j) { // for each stripe in the query
              __tmap_vsw_mm_store_si128(H, __tmap_vsw_mm_load_si128(H1 + j));
              H++;
          }
      }
      /*
      fprintf(stderr, "i=%d", i);
      for(j = 0; j < tmap_vsw8_values_per_128_bits; j++) { // for each start position in the stripe
          for(k = 0; k < query->slen; k++) {
              fprintf(stderr, " %d", ((tmap_vsw8_uint_t*)(H1+k))[j]);
          }
      }
      fprintf(stderr, "\n");
      for(k=0;k<tmap_vsw8_values_per_128_bits;++k) {
          fprintf(stderr, "%d ", ((tmap_vsw8_uint_t*)&max)[k]);
      }
      fprintf(stderr, "\n");
      */
      __tmap_vsw8_max(imax, max); // imax is the maximum number in max
      if(imax > gmax) {
          int32_t m;
          gmax = imax; // global maximum score 
          result->target_end = i;
          tmap_vsw8_uint_t *t = (tmap_vsw8_uint_t*)H1;
          for(j = 0, result->query_end = m = -1; TMAP_VSW_LIKELY(j < slen); ++j) { // for each stripe in the query
              for(k = 0; k < tmap_vsw8_values_per_128_bits; k++, t++) { // for each cell in the stripe
                  if((int32_t)*t > m) { // find the first instance
                      m = (int32_t)*t;
                      result->query_end = j + ((k & (tmap_vsw8_values_per_128_bits-1)) * slen);
                  }
              }
          }
      }
      if (gmax >= _thres) { // reduce
          if(NULL != overflow) {
              *overflow = 1;
              return;
          }
          // When overflow is going to happen, subtract tmap_vsw8_mid_value from all scores. This heuristic
          // may miss the best alignment, but in practice this should happen very rarely.
          sum += tmap_vsw8_mid_value; gmax -= tmap_vsw8_mid_value;
          for(j = 0; TMAP_VSW_LIKELY(j < slen); ++j) {
              h = __tmap_vsw8_mm_subs_epu(__tmap_vsw_mm_load_si128(H1 + j), reduce);
              __tmap_vsw_mm_store_si128(H1 + j, h);
              e = __tmap_vsw8_mm_subs_epu(__tmap_vsw_mm_load_si128(E + j), reduce);
              __tmap_vsw_mm_store_si128(E + j, e);
          }
      }
      S = H1; H1 = H0; H0 = S; // swap H0 and H1
  }
  result->score_fwd = gmax + sum;
}

// H -> diagonal
// E -> deletion (query)
// F -> insertion (query)
// H' -> max{H, E, F}
void
tmap_vsw8_sse2_reverse(tmap_vsw8_query_t *query, const uint8_t *target, int32_t tlen, tmap_vsw_opt_t *opt, tmap_vsw_result_t *result, int32_t *overflow)
{
  int32_t slen, i, j, k, sum, m_b, n_b, _thres = tmap_vsw8_max_value - query->range_score - 2;
  int32_t gmax = 0;
  __m128i zero, pen_gapoe, pen_gape, min_score, reduce, thres, *H0, *H1, *E;

  // initialization
  zero = _mm_set1_epi32(0); // mmz zeros
  m_b = n_b = 0; 
  pen_gapoe = __tmap_vsw8_mm_set_epi(opt->pen_gapo + opt->pen_gape); // gap open penalty
  pen_gape = __tmap_vsw8_mm_set_epi(opt->pen_gape); // gap extend penalty
  min_score = __tmap_vsw8_mm_set_epi(query->min_score); // the minimum score
  thres = __tmap_vsw8_mm_set_epi(_thres); // when max score exceeds this, min_score all scores by "reduce" below
  reduce = __tmap_vsw8_mm_set_epi(tmap_vsw8_mid_value);
  H0 = query->H0; 
  H1 = query->H1; 
  E = query->E; 
  // re-initialize the query
  tmap_vsw8_query_reinit(query, result);
  // get the new ranges
  slen = query->slen;
  tlen = (0 <= result->target_end) ? (result->target_end+1) : tlen;
  if(NULL != overflow) *overflow = 0;

  // initialize all zeros
  for(i = 0; i < slen; ++i) {
      __tmap_vsw_mm_store_si128(E + i, zero);
      __tmap_vsw_mm_store_si128(H0 + i, zero);
  }
  
  // HERE
  //fprintf(stderr, "\n");
  //tmap_vsw8_query_print_query_profile(query);
  // select query end only
  // check stripe #: __tmap_vsw8_query_index_to_stripe_number(0, slen) 
  // check byte # __tmap_vsw8_query_index_to_byte_number(0, slen)
  
  for(i = tlen-1, sum = 0; 0 <= i; i--) { // for each base in the target
      int cmp, imax = 0;
      // NB: f and max are initially zero
      __m128i e, h, f = zero, max = zero, *S = query->query_profile + target[i] * slen; // s is the 1st score vector
      h = __tmap_vsw_mm_load_si128(H0); // the first stripe
      h = __tmap_vsw_mm_srli_si128(h, tmap_vsw8_shift_bits); // h=H(i+1,+1); >> instead of << because x64 is little-endian
      for(j = slen-1; TMAP_VSW_LIKELY(0 <= j); j--) { // for each stripe in the query
          /* SW cells are computed in the following order:
           *   H(i,j)   = max{H(i+1,j+1)+S(i,j), E(i,j), F(i,j)}
           *   E(i-1,j) = max{H(i,j)-q, E(i,j)-r}
           *   F(i,j-1) = max{H(i,j)-q, F(i,j)-r}
           */
          // compute H'(i,j); note that at the beginning, h=H'(i+1,j+1)
          h = __tmap_vsw8_mm_adds_epu(h, __tmap_vsw_mm_load_si128(S + j));
          h = __tmap_vsw8_mm_subs_epu(h, min_score); // h=H'(i+1,j+1)+S(i,j)
          e = __tmap_vsw_mm_load_si128(E + j); // e=E'(i,j)
          h = __tmap_vsw8_mm_max_epu(h, e); // h=H'(i,j) = max{E'(i,j), H'(i+1,j+1)+S(i,j)}
          h = __tmap_vsw8_mm_max_epu(h, f); // h=H'(i,j) = max{max{E'(i,j), H'(i+1,j+1)+S(i,j)}, F}
          max = __tmap_vsw8_mm_max_epu(max, h); // save the max values in this stripe versus the last
          __tmap_vsw_mm_store_si128(H1 + j, h); // save h to H'(i,j)
          // next, compute E'(i-1,j)
          h = __tmap_vsw8_mm_subs_epu(h, pen_gapoe); // h=H'(i,j)-pen_gapoe
          e = __tmap_vsw8_mm_subs_epu(e, pen_gape); // e=E'(i,j)-pen_gape
          e = __tmap_vsw8_mm_max_epu(e, h); // e=E'(i-1,j) = max{E'(i,j)-pen_gape, H'(i,j)-pen_gapo}
          __tmap_vsw_mm_store_si128(E + j, e); // save e to E'(i-1,j)
          // now compute F'(i,j-1)
          //h = __tmap_vsw8_mm_subs_epu(h, pen_gapoe); // h=H'(i,j)-pen_gapoe
          f = __tmap_vsw8_mm_subs_epu(f, pen_gape); // f=F'(i,j)-pen_gape
          f = __tmap_vsw8_mm_max_epu(f, h); // f=F'(i,j-1) = max{F'(i,j)-pen_gape, H'(i,j)-pen_gapo}
          // get H'(i+1,j) and prepare for the next j
          h = __tmap_vsw_mm_load_si128(H0 + j); // h=H'(i+1,j)
      }
      // NB: we do not need to set E(i,j) as we disallow adjacent insertion and then deletion
      for(k=tmap_vsw8_values_per_128_bits-1; TMAP_VSW_LIKELY(0 <= k); k--) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
          f = __tmap_vsw_mm_srli_si128(f, tmap_vsw8_shift_bits);
          for(j = slen-1; TMAP_VSW_LIKELY(0 <= j); j--) {
              h = __tmap_vsw_mm_load_si128(H1 + j);
              h = __tmap_vsw8_mm_max_epu(h, f); // h=H'(i,j)
              __tmap_vsw_mm_store_si128(H1 + j, h);
              h = __tmap_vsw8_mm_subs_epu(h, pen_gapoe);
              f = __tmap_vsw8_mm_subs_epu(f, pen_gape);
              cmp = __tmap_vsw8_mm_movemask_epi(_mm_cmpeq_epi8(__tmap_vsw8_mm_subs_epu(f, h), zero));
              if (TMAP_VSW_UNLIKELY(cmp == 0xffff)) goto end_loop; // ACK: goto statement
          }
      }
end_loop:
      /*
      fprintf(stderr, "i=%3d b=%1d", i, target[i]);
      for(j = tmap_vsw8_values_per_128_bits-1, l = 0; 0 <= j; j--) { // for each start position in the stripe
          for(k = query->slen-1; 0 <= k; k--, l++) {
              fprintf(stderr, " %2d", ((tmap_vsw8_uint_t*)(H1+k))[j]);
          }
      }
      fprintf(stderr, "\n");
      */
      /*
      for(j = tmap_vsw8_values_per_128_bits-1, l = 0; 0 <= j; j--) { // for each start position in the stripe
          for(k = query->slen-1; 0 <= k; k--, l++) {
              fprintf(stderr, "i=%3d j=%3d tb=%1d qs=%1d %2d\n", 
                      i, ((query->qlen + 15) >> 4 << 4)-l-1, target[i],
                      ((tmap_vsw8_uint_t*)(query->query_profile+k))[j] - query->min_score,
                      ((tmap_vsw8_uint_t*)(H1+k))[j]);
          }
      }
      */
      /*
      for(k=0;k<tmap_vsw8_values_per_128_bits;++k) {
          fprintf(stderr, "%d ", ((tmap_vsw8_uint_t*)&max)[k]);
      }
      fprintf(stderr, "\n");
      */
      __tmap_vsw8_max(imax, max); // imax is the maximum number in max
      if(imax > gmax) {
          int32_t m;
          gmax = imax; // global maximum score 
          result->target_start = i;
          tmap_vsw8_uint_t *t = (tmap_vsw8_uint_t*)H1;
          for(j = 0, result->query_start = m = -1; TMAP_VSW_LIKELY(j < slen); ++j) { // for each stripe in the query
              for(k = 0; k < tmap_vsw8_values_per_128_bits; k++, t++) { // for each cell in the stripe
                  if((int32_t)*t > m) { // find the first instance
  // select query end only
                      m = (int32_t)*t;
                      result->query_start = j + ((k & (tmap_vsw8_values_per_128_bits-1)) * slen);
                  }
              }
          }
      }
      if (gmax >= _thres) { // reduce
          if(NULL != overflow) {
              *overflow = 1;
              return;
          }
          // When overflow is going to happen, subtract tmap_vsw8_mid_value from all scores. This heuristic
          // may miss the best alignment, but in practice this should happen very rarely.
          sum += tmap_vsw8_mid_value; gmax -= tmap_vsw8_mid_value;
          for(j = 0; TMAP_VSW_LIKELY(j < slen); ++j) {
              h = __tmap_vsw8_mm_subs_epu(__tmap_vsw_mm_load_si128(H1 + j), reduce);
              __tmap_vsw_mm_store_si128(H1 + j, h);
              e = __tmap_vsw8_mm_subs_epu(__tmap_vsw_mm_load_si128(E + j), reduce);
              __tmap_vsw_mm_store_si128(E + j, e);
          }
      }
      S = H1; H1 = H0; H0 = S; // swap H0 and H1
  }
  
  result->score_rev = gmax + sum;
}

void
tmap_vsw8_sse2(tmap_vsw8_query_t *query, const uint8_t *target, int32_t tlen, tmap_vsw_opt_t *opt, tmap_vsw_result_t *result, int32_t *overflow)
{
  // forward
  tmap_vsw8_sse2_forward(query, target, tlen, opt, result, overflow);
  /*
  fprintf(stderr, "qlen=%d tlen=%d\n", query->qlen, tlen);
  fprintf(stderr, "result = {%d-%d} {%d-%d} score=%d\n",
          result->query_start, result->query_end,
          result->target_start, result->target_end, result->score_fwd);
          */
  if(NULL != overflow && 1 == *overflow) return;
  // adjust lengths
  query->qlen = result->query_end+1;
  tlen = result->target_end+1;
  // reverse
  tmap_vsw8_sse2_reverse(query, target, tlen, opt, result, overflow);
  /*
  fprintf(stderr, "qlen=%d tlen=%d\n", query->qlen, tlen);
  fprintf(stderr, "result = {%d-%d} {%d-%d} score=%d\n",
          result->query_start, result->query_end,
          result->target_start, result->target_end, result->score_rev);
  */
  // check that they agree
  if(result->score_fwd != result->score_rev) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }
}

void
tmap_vsw8_sse2_get_path(const uint8_t *query, int32_t qlen, 
                              const uint8_t *target, int32_t tlen, 
                              tmap_vsw8_query_t *q,
                              tmap_vsw_result_t *result, 
                              tmap_sw_path_t *path,
                              int32_t *path_len,
                              tmap_vsw_opt_t *opt)
{
  int32_t ti, qi;
  int32_t pen_gapoe, pen_gape;
  int32_t score, overflow = 0;
  tmap_sw_path_t *p;
  // Debugging
  //int32_t slen, i, j, k, l;
  //__m128i *H;
  
  if(1 != q->type) { // ignore
      tmap_error("bug encountered", Exit, OutOfRange);
  }

  if(0 == qlen || 0 == tlen) {
      *path_len = 0;
      return; // smallest value
  }

  // store here
  pen_gapoe = opt->pen_gapo + opt->pen_gape; // gap open penalty
  pen_gape = opt->pen_gape; // gap extend penalty

  // initialize query
  q = tmap_vsw8_query_init_full(q, query, qlen, tlen, opt);
      
  /*
  fprintf(stderr, "qlen=%d tlen=%d\n", qlen, tlen);
  fprintf(stderr, "result = {%d-%d} {%d-%d}\n",
          result->query_start, result->query_end,
          result->target_start, result->target_end);
          */
  
  // run the forward VSW
  score = result->score_fwd;
  tmap_vsw8_sse2_forward(q, target, tlen, opt, result, &overflow);
  if(1 == overflow) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }
  if(score != result->score_fwd) {
      fprintf(stderr, "result = {%d-%d} {%d-%d}\n",
              result->query_start, result->query_end,
              result->target_start, result->target_end);
      fprintf(stderr, "score=%d result->score_fwd=%d result->score_rev=%d\n", score, result->score_fwd, result->score_rev);
      tmap_error("bug encountered", Exit, OutOfRange);
  }

  // trackback
  ti = tlen-1; qi = q->qlen-1;
  p = path;
  while(0 <= ti || 0 <= qi) {
      tmap_vsw8_uint_t h_cur, h_prev, e_prev, f_prev;

      // get the cells to compare
      h_cur = __tmap_vsw8_get_matrix_cell(q, (ti < 0) ? 0 : ti, (qi < 0) ? 0 : qi);
      h_prev = __tmap_vsw8_get_matrix_cell(q, (ti-1 < 0) ? 0 : ti-1, (qi-1 < 0) ? 0 : qi-1);
      e_prev = __tmap_vsw8_get_matrix_cell(q, (ti-1 < 0) ? 0 : ti-1, (qi < 0) ? 0 : qi);
      f_prev = __tmap_vsw8_get_matrix_cell(q, (ti < 0) ? 0 : ti, (qi-1 < 0) ? 0 : qi-1);

      // get the match score
      score = __tmap_vsw8_get_query_profile_value(q, target[(ti < 0) ? 0 : ti], (qi < 0) ? 0 : qi) - q->min_score;
      
      /*
      fprintf(stderr, "h_cur=%u h_prev=%u e_prev=%u f_prev=%u score=%d\n",
              h_cur, h_prev, e_prev, f_prev, score);
              */
      
      // one-based
      p->i = ti+1; p->j = qi+1; 

      // compare
      if(h_cur == h_prev + score) {
          //fprintf(stderr, "M ti=%d qi=%d\n", ti, qi);
          p->ctype = TMAP_SW_FROM_M;
          ti--;
          qi--;
      }
      else if(h_cur == e_prev - pen_gapoe || h_cur == e_prev - pen_gape) {
          //fprintf(stderr, "D ti=%d qi=%d\n", ti, qi);
          p->ctype = TMAP_SW_FROM_D;
          ti--;
      }
      else if(h_cur == f_prev - pen_gapoe || h_cur == f_prev - pen_gape) {
          //fprintf(stderr, "I ti=%d qi=%d\n", ti, qi);
          p->ctype = TMAP_SW_FROM_M;
          qi--;
      }
      else if(0 == ti && 0 == qi) { // TODO: wrong for different initializations...
          p->ctype = TMAP_SW_FROM_M;
          ti--;
          qi--;
      }
      else {
          // TODO: gaps etc.
          fprintf(stderr, "ti=%d qi=%d\n", ti, qi);
          fprintf(stderr, "pen_gapoe=%d pen_gape=%d\n", pen_gapoe, pen_gape);
          tmap_error("bug encountered", Exit, OutOfRange);
      }

      // move to the next path etc.
      p++;
  }
  *path_len = p - path;
}
