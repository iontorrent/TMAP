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
#include <unistd.h>
#include <stdio.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "tmap_sw.h"
#include "tmap_vsw.h"
#include "tmap_vsw16.h"

// TODO: remove query profile so it does not need to be filled in each time

/*
static void
tmap_vsw16_query_print_query_profile(tmap_vsw16_query_t *query)
{
  int32_t a, i, j, k;

  fprintf(stderr, "min_edit_score=%d\n", query->min_edit_score);
  for(a=0;a<TMAP_VSW_ALPHABET_SIZE;a++) {
      fprintf(stderr, "a=%d ", a);
      __m128i *S = query->query_profile + a * query->slen;
      for(i = k = 0; i < tmap_vsw16_values_per_128_bits; i++) {
          for(j = 0; j < query->slen; j++, k++) {
              fprintf(stderr, " %d:%d", k, ((tmap_vsw16_int_t*)(S+j))[i]);
          }
      }
      fprintf(stderr, "\n");
  }
}
*/

tmap_vsw16_query_t *
tmap_vsw16_query_init(tmap_vsw16_query_t *prev, const uint8_t *query, int32_t qlen, int32_t tlen, 
                              int32_t query_start_clip, int32_t query_end_clip,
                              tmap_vsw_opt_t *opt)
{
  int32_t qlen_mem, slen, a;
      int32_t min_gap_score, min_mm_score;
  tmap_vsw16_int_t *t;

  // get the number of stripes
  slen = __tmap_vsw16_calc_slen(qlen); 

  // check if we need to re-size the memory
  if(NULL == prev 
     || prev->qlen_max < qlen
     || query_start_clip != prev->query_start_clip
     || query_end_clip != prev->query_end_clip) { // recompute
      // free previous
      if(NULL != prev) {
          free(prev); prev = NULL;
      }
      // get the memory needed to hold the stripes for the query 
      qlen_mem = __tmap_vsw_calc_qlen_mem(slen);
      //prev = tmap_memalign(16, sizeof(tmap_vsw16_query_t) + 15 + qlen_mem * (TMAP_VSW_ALPHABET_SIZE + 3), "prev"); // add three for H0, H1, and E
      prev = tmap_malloc(sizeof(tmap_vsw16_query_t) + 15 + qlen_mem * (TMAP_VSW_ALPHABET_SIZE + 3), "prev"); // add three for H0, H1, and E
      prev->qlen_max = qlen; 
      // update the memory
      prev->qlen_mem = qlen_mem;
      // update clipping
      prev->query_start_clip = query_start_clip;
      prev->query_end_clip = query_end_clip;
  }

  // compute min/max edit scores
  prev->min_edit_score = (opt->pen_gapo+opt->pen_gape < opt->pen_mm) ? -opt->pen_mm : -(opt->pen_gapo+opt->pen_gape); // minimum single-edit score
  prev->min_edit_score--; // for N mismatches
  prev->max_edit_score = opt->score_match;
  // compute the zero alignment scores (assume no clipping in case of re-use)
  min_mm_score = qlen * -opt->pen_mm; // all mismatches
  min_gap_score = -opt->pen_gapo + (-opt->pen_gape * qlen); // all insertions
  prev->zero_aln_score = (min_mm_score < min_gap_score) ? min_mm_score : min_gap_score;
  // the minimum alignment score
  prev->min_aln_score = prev->zero_aln_score << 1; // double it so that we have room in case of starting at negative infinity
  prev->min_aln_score += prev->min_edit_score; // so it doesn't underflow
  // max aln score
  prev->max_aln_score = tmap_vsw16_max_value - prev->max_edit_score; // so it doesn't overflow
  // normalize
  prev->zero_aln_score = tmap_vsw16_min_value - prev->min_aln_score- prev->zero_aln_score; 
  prev->min_aln_score = tmap_vsw16_min_value - prev->min_aln_score; 
  // normalize with the minimum alignment score
  /*
  fprintf(stderr, "min_edit_score=%d\n", prev->min_edit_score);
  fprintf(stderr, "max_edit_score=%d\n", prev->max_edit_score);
  fprintf(stderr, "zero_aln_score=%d\n", prev->zero_aln_score);
  fprintf(stderr, "min_aln_score=%d\n", prev->min_aln_score);
  fprintf(stderr, "max_aln_score=%d\n", prev->max_aln_score);
  */

  prev->qlen = qlen; // update the query length
  prev->slen = slen; // update the number of stripes

  // NB: align all the memory from one block
  prev->query_profile = (__m128i*)__tmap_vsw_16((size_t)prev + sizeof(tmap_vsw16_query_t)); // skip over the struct variables, align memory 
  prev->H0 = prev->query_profile + (prev->slen * TMAP_VSW_ALPHABET_SIZE); // skip over the query profile
  prev->H1 = prev->H0 + prev->slen; // skip over H0
  prev->E = prev->H1 + prev->slen; // skip over H1
  
  // create the query profile
  t = (tmap_vsw16_int_t*)prev->query_profile;
  for(a = 0; a < TMAP_VSW_ALPHABET_SIZE; a++) {
      int32_t i, k;
      for(i = 0; i < prev->slen; ++i) { // for each stripe
          // fill in this stripe
          for(k = i; k < prev->slen << tmap_vsw16_values_per_128_bits_log2; k += prev->slen) { //  for q_{i+1}, q_{2i+1}, ..., q_{2s+1}
              // NB: pad with zeros
              *t++ = ((k >= qlen) ? prev->min_edit_score : ((a == query[k]) ? opt->score_match : -opt->pen_mm)); 
          }
      }
  }
  //tmap_vsw16_query_print_query_profile(prev);
  return prev;
}

void
tmap_vsw16_query_destroy(tmap_vsw16_query_t *vsw)
{
  // all memory was allocated as one block
  free(vsw);
}

int32_t
tmap_vsw16_sse2_forward(tmap_vsw16_query_t *query, const uint8_t *target, int32_t tlen, 
                        int32_t query_start_clip, int32_t query_end_clip,
                        tmap_vsw_opt_t *opt, int16_t *query_end, int16_t *target_end,
                        int32_t direction, int32_t *overflow, int32_t score_thr)
{
  int32_t slen, i, j, k, sum = 0;
#ifdef TMAP_VSW_DEBUG
  int32_t l;
#endif
  uint16_t cmp;
  tmap_vsw16_int_t gmax, best;
  tmap_vsw16_int_t zero, imin = 0, imax = 0;
  __m128i zero_mm, zero_start_mm, negative_infinity_mm, positive_infinity_mm, reduce_mm;
  __m128i pen_gapoe, pen_gape, *H0, *H1, *E;

  // initialization
  // normalize these
  zero = query->zero_aln_score; // where the normalized zero alignment score occurs
  negative_infinity_mm = __tmap_vsw16_mm_set1_epi16(query->min_aln_score); // the minimum possible value
  positive_infinity_mm = __tmap_vsw16_mm_set1_epi16(query->max_aln_score); // the minimum possible value
  score_thr += zero; // for the scoring threshold
  // these are not normalized
  pen_gapoe = __tmap_vsw16_mm_set1_epi16(opt->pen_gapo + opt->pen_gape); // gap open penalty
  pen_gape = __tmap_vsw16_mm_set1_epi16(opt->pen_gape); // gap extend penalty
  zero_mm = __tmap_vsw16_mm_set1_epi16(zero); // where the normalized zero alignment score occurs
  zero_start_mm = __tmap_vsw16_mm_set1_epi16(zero); // where the normalized zero alignment score occurs
  gmax = tmap_vsw16_min_value;
  best = tmap_vsw16_min_value;
  reduce_mm = __tmap_vsw16_mm_set1_epi16(tmap_vsw16_mid_value);
  // vectors
  H0 = query->H0; 
  H1 = query->H1; 
  E = query->E; 
  slen = query->slen;
  if(NULL != overflow) *overflow = 0;
#ifdef TMAP_VSW_DEBUG
  fprintf(stderr, "qlen=%d tlen=%d\n", query->qlen, tlen);
  fprintf(stderr, "query_start_clip=%d query_end_clip=%d\n",
          query_start_clip, query_end_clip);
#endif

  // select query end only
  // check stripe #: __tmap_vsw16_query_index_to_stripe_number(query->qlen, slen) 
  // check byte # __tmap_vsw16_query_index_to_byte_number(query->qlen, slen)

  if(0 == query_start_clip) {
      __m128i f;
      // Set start
      for(j = 0; j < slen; j++) {
          // initialize E to negative infinity 
          __tmap_vsw_mm_store_si128(E + j, negative_infinity_mm); // E(0,j)
      }
      // NB: setting the start == 0 will be done later
      // set the leading insertions 
      f = __tmap_vsw16_mm_set1_epi16(query->min_aln_score);
      f = __tmap_vsw16_mm_insert_epi16(f, zero, 0); 
      for(j = 0; j < tmap_vsw16_values_per_128_bits; ++j) {
          f = __tmap_vsw16_mm_insert(f , -opt->pen_gapo + -opt->pen_gape + (-opt->pen_gape * slen * j) + zero, j);
      }
      for(j = 0; j < slen; j++) {
          f = __tmap_vsw16_mm_subs_epi16(f, pen_gape); // f=F(i,j)-pen_gape
          f = __tmap_vsw16_mm_max_epi16(f, negative_infinity_mm); // bound with -inf
          __tmap_vsw_mm_store_si128(H0 + ((slen + j - 1) % slen), f); // H(0,j)
      }
  }
  else {
      // Initialize all zeros
      for(i = 0; i < slen; ++i) {
          __tmap_vsw_mm_store_si128(E + i, zero_mm);
          __tmap_vsw_mm_store_si128(H0 + i, zero_mm);
      }
  }

  // the core loop
  for(i = 0, sum = 0; i < tlen; i++) { // for each base in the target
      __m128i e, h, f, g, max, min, *S;

      max = __tmap_vsw16_mm_set1_epi16(query->min_aln_score); // max is negative infinity
      min = __tmap_vsw16_mm_set1_epi16(query->max_aln_score); // min is positive infinity
      S = query->query_profile + target[i] * slen; // s is the 1st score vector
  
      // load H(i-1,-1)
      h = __tmap_vsw_mm_load_si128(H0 + slen - 1); // the last stripe, which holds the j-1 
      if(TMAP_VSW_UNLIKELY(0 < i)) { // only if we have previous results
          h = __tmap_vsw_mm_slli_si128(h, tmap_vsw16_shift_bytes); // shift left since x64 is little endian
      }
      h = __tmap_vsw16_mm_insert_epi16(h, zero, 0);
      
      // e does not need to be set

      // set F to -inf
      f = __tmap_vsw16_mm_set1_epi16(query->min_aln_score);
      // leading insertions
      g = __tmap_vsw16_mm_set1_epi16(query->min_aln_score);
      if(0 == query_start_clip) { 
          g = __tmap_vsw16_mm_set1_epi16(query->min_aln_score);
          g = __tmap_vsw16_mm_insert_epi16(g, -opt->pen_gapo + opt->pen_gape + zero, 0);
          for(j = 1; TMAP_VSW_LIKELY(j < tmap_vsw16_values_per_128_bits); j++) {
              // NB: do not include the gap extend, that will be added below
              // BEFORE considering H
              g = __tmap_vsw16_mm_insert(g, -opt->pen_gapo + (-opt->pen_gape * (j * slen - 1)) + zero, j);
          }
      }
              
      for(j = 0; TMAP_VSW_LIKELY(j < slen); ++j) { // for each stripe in the query
          // NB: at the beginning, 
          // h=H(i-1,j-1)
          // e=E(i,j)
          // f=F(i,j)
          /* SW cells are computed in the following order:
           *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
           *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
           *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
           */
          if(1 == query_start_clip) {
              // start anywhere within the query, though not before the start
              h = __tmap_vsw16_mm_max_epi16(h, zero_start_mm);  
          }
          else {
              g = __tmap_vsw16_mm_subs_epi16(g, pen_gape); // leading insertion
#ifdef TMAP_VSW_DEBUG 
              /*
              fprintf(stderr, "h i=%d j=%d", i, j);
              for(k = 0; k < tmap_vsw16_values_per_128_bits; k++) { // for each start position in the stripe
                  fprintf(stderr, " (%d,%d)", 
                          ((tmap_vsw16_int_t*)(&h))[k] - zero,
                          ((tmap_vsw16_int_t*)(&g))[k] - zero);
              }
              fprintf(stderr, "\n");
              */
#endif
              h = __tmap_vsw16_mm_max_epi16(h, g);
          }
          // compute H(i,j); 
#ifdef TMAP_VSW_DEBUG 
          /*
          __m128i s = __tmap_vsw_mm_load_si128(S + j);
          fprintf(stderr, "s i=%d j=%d", i, j);
          for(k = 0; k < tmap_vsw16_values_per_128_bits; k++) { // for each start position in the stripe
              fprintf(stderr, " %4d", ((tmap_vsw16_int_t*)(&s))[k]);
          }
          fprintf(stderr, "\n");
          fprintf(stderr, "h i=%d j=%d", i, j);
          for(k = 0; k < tmap_vsw16_values_per_128_bits; k++) { // for each start position in the stripe
              fprintf(stderr, " %4d", ((tmap_vsw16_int_t*)(&h))[k] - zero);
          }
          fprintf(stderr, "\n");
          */
#endif
          h = __tmap_vsw16_mm_adds_epi16(h, __tmap_vsw_mm_load_si128(S + j)); // h=H(i-1,j-1)+S(i,j)
          e = __tmap_vsw_mm_load_si128(E + j); // e=E(i,j)
          h = __tmap_vsw16_mm_max_epi16(h, e); // h=H(i,j) = max{E(i,j), H(i-1,j-1)+S(i,j)}
          h = __tmap_vsw16_mm_max_epi16(h, f); // h=H(i,j) = max{max{E(i,j), H(i-1,j-1)+S(i,j)}, F(i,j)}
          h = __tmap_vsw16_mm_max_epi16(h, negative_infinity_mm); // bound with -inf
          max = __tmap_vsw16_mm_max_epi16(max, h); // save the max values in this stripe versus the last
          min = __tmap_vsw16_mm_min_epi16(min, h); // save the max values in this stripe versus the last
          __tmap_vsw_mm_store_si128(H1 + j, h); // save h to H(i,j)
          // next, compute E(i+1,j)
          h = __tmap_vsw16_mm_subs_epi16(h, pen_gapoe); // h=H(i,j)-pen_gapoe
          e = __tmap_vsw16_mm_subs_epi16(e, pen_gape); // e=E(i,j)-pen_gape
          e = __tmap_vsw16_mm_max_epi16(e, h); // e=E(i+1,j) = max{E(i,j)-pen_gape, H(i,j)-pen_gapoe}
          e = __tmap_vsw16_mm_max_epi16(e, negative_infinity_mm); // bound with -inf
          __tmap_vsw_mm_store_si128(E + j, e); // save e to E(i+1,j)
          // now compute F(i,j+1)
          //h = __tmap_vsw16_mm_subs_epi16(h, pen_gapoe); // h=H(i,j)-pen_gapoe
          f = __tmap_vsw16_mm_subs_epi16(f, pen_gape); // f=F(i,j)-pen_gape
          f = __tmap_vsw16_mm_max_epi16(f, h); // f=F(i,j+1) = max{F(i,j)-pen_gape, H(i,j)-pen_gapoe}
          f = __tmap_vsw16_mm_max_epi16(f, negative_infinity_mm); // bound with -inf
          // get H(i-1,j) and prepare for the next j
          h = __tmap_vsw_mm_load_si128(H0 + j); // h=H(i-1,j)
      }
      
      // Farrar's Lazy F-loop evaluation
      f = __tmap_vsw16_mm_set1_epi16(query->min_aln_score);
      f = __tmap_vsw_mm_slli_si128(f, tmap_vsw16_shift_bytes); // since x86 is little endian
      f = __tmap_vsw16_mm_insert_epi16(f, query->min_aln_score, 0); // set F(i-1,-1)[0] as negative infinity (normalized)
      k = 0;
      while(TMAP_VSW_UNLIKELY(k < tmap_vsw16_values_per_128_bits)) { // do not need so many iterations
          // check if h will be changed
          h = __tmap_vsw_mm_load_si128(H1 + j); // h=H(i,j)
          h = __tmap_vsw16_mm_subs_epi16(h, pen_gapoe); // h=H(i,j)-gapo
          cmp = __tmap_vsw16_mm_movemask_epi16(__tmap_vsw16_mm_cmplt_epi16(f, h));
          if (TMAP_VSW_UNLIKELY(cmp == 0xffff)) break; // h is not changed
          // changed, so save
          __tmap_vsw_mm_store_si128(H1 + j, f); // save h to H(i,j) 
          f = __tmap_vsw16_mm_subs_epi16(f, pen_gape); // f=F(i,j)-pen_gape
          f = __tmap_vsw16_mm_max_epi16(f, h); // f=F(i,j+1) = max{F(i,j)-pen_gape, H(i,j)-pen_gapoe}
          f = __tmap_vsw16_mm_max_epi16(f, negative_infinity_mm); // bound with -inf
          j++;
          if(TMAP_VSW_UNLIKELY(slen == j)) { // entire segment has been processed
              f = __tmap_vsw_mm_slli_si128(f, tmap_vsw16_shift_bytes); // since x86 is little endian
              f = __tmap_vsw16_mm_insert_epi16(f, query->min_aln_score, 0); // set F(i-1,-1)[0] as negative infinity (normalized)
              j = 0;
              k++;
          }
      }
      /*
      // NB: we do not need to set E(i,j) as we disallow adjacent insertion and then deletion
      // iterate through each value stored in a stripe
      f = __tmap_vsw16_mm_set1_epi16(query->min_aln_score);
      // we require at most 'tmap_vsw16_values_per_128_bits' iterations to guarantee F propagation
      for(k = 0; TMAP_VSW_LIKELY(k < tmap_vsw16_values_per_128_bits); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
          f = __tmap_vsw_mm_slli_si128(f, tmap_vsw16_shift_bytes); // since x86 is little endian
          f = __tmap_vsw16_mm_insert_epi16(f, query->min_aln_score, 0); // set F(i-1,-1)[0] as negative infinity (normalized)
          for(j = 0; TMAP_VSW_LIKELY(j < slen); ++j) { 
              h = __tmap_vsw_mm_load_si128(H1 + j); // h=H(i,j)
              h = __tmap_vsw16_mm_max_epi16(h, f); // h=H(i,j) = max{H(i,j), F(i,j)}
              h = __tmap_vsw16_mm_max_epi16(h, negative_infinity_mm); // bound with -inf
              __tmap_vsw_mm_store_si128(H1 + j, h); // save h to H(i,j) 
              h = __tmap_vsw16_mm_subs_epi16(h, pen_gapoe); // h=H(i,j)-gapo
              f = __tmap_vsw16_mm_subs_epi16(f, pen_gape); // f=F(i,j)-pen_gape
              f = __tmap_vsw16_mm_max_epi16(f, h); // f=F(i,j+1) = max{F(i,j)-pen_gape, H(i,j)-pen_gapoe}
              f = __tmap_vsw16_mm_max_epi16(f, negative_infinity_mm); // bound with -inf
              // check to see if h could have been updated by f?
              // NB: the comparison below will have some false positives, in
              // other words h == f a priori, but this is rare.
              uint16_t cmp_old, cmp_new;
              cmp_old = __tmap_vsw16_mm_movemask_epi16(__tmap_vsw16_mm_cmpeq_epi16(__tmap_vsw16_mm_subs_epi16(f, h), zero_mm));
              cmp_new = __tmap_vsw16_mm_movemask_epi16(__tmap_vsw16_mm_cmpeq_epi16(f, h));
              if(cmp_old != cmp_new) {
                  fprintf(stderr, "cmp=[%d,%d]\n", cmp_old, cmp_new);
                  tmap_error("bug encountered", Exit, OutOfRange);
              }
              cmp = cmp_old;
              if (TMAP_VSW_UNLIKELY(cmp == 0xffff)) goto end_loop; // ACK: goto statement
          }
          *?
      }
end_loop:
      */
#ifdef TMAP_VSW_DEBUG 
      fprintf(stderr, "H1 i=%d target[i]=%d", i, target[i]);
      for(k = l = 0; k < tmap_vsw16_values_per_128_bits; k++) { // for each start position in the stripe
          for(j = 0; j < slen; j++, l++) {
              fprintf(stderr, " %d:%d", l, ((tmap_vsw16_int_t*)(H1 + j))[k] - zero);
          }
      }
      fprintf(stderr, "\n");
#endif
      __tmap_vsw16_max(imax, max); // imax is the maximum number in max
      if(query->max_aln_score - query->max_edit_score < imax) { // overflow
          tmap_error("bug encountered", Exit, OutOfRange);
      }
      if(imax > gmax) { 
          gmax = imax; // global maximum score 
      }
      if(score_thr <= imax 
         && (imax > best || (1 == direction && imax == best))) { // potential best score
          tmap_vsw16_int_t *t;
          if(query_end_clip == 0) { // check the last
              j = (query->qlen-1) % slen; // stripe
              k = (query->qlen-1) / slen; // byte
              t = (tmap_vsw16_int_t*)(H1 + j);
#ifdef TMAP_VSW_DEBUG 
              fprintf(stderr, "j=%d k=%d slen=%d qlen=%d best=%d pot=%d\n", j, k, slen, query->qlen, best-zero, t[k]-zero);
#endif
              if((int32_t)t[k] > best || (1 == direction && (int32_t)t[k] == best)) { // found
                  (*query_end) = query->qlen-1;
                  (*target_end) = i;
                  best = t[k];
#ifdef TMAP_VSW_DEBUG 
                  fprintf(stderr, "FOUND A i=%d query->qlen-1=%d query_end=%d target_end=%d best=%d\n", i, query->qlen-1, *query_end, *target_end, best-zero);
#endif
              }
          }
          else { // check all
              t = (tmap_vsw16_int_t*)H1;
              (*target_end) = i;
              for(j = 0, *query_end = -1; TMAP_VSW_LIKELY(j < slen); ++j) { // for each stripe in the query
                  for(k = 0; k < tmap_vsw16_values_per_128_bits; k++, t++) { // for each cell in the stripe
                      if((int32_t)*t > best || (1 == direction && (int32_t)*t== best)) { // found
                          best = *t;
                          *query_end = j + ((k & (tmap_vsw16_values_per_128_bits-1)) * slen);
#ifdef TMAP_VSW_DEBUG 
                          fprintf(stderr, "FOUND B j=%d k=%d query_end=%d best=%d\n", j, k, *query_end, best-zero);
#endif
                      }
                  }
              }
              if(-1 == *query_end) {
                  tmap_error("bug encountered", Exit, OutOfRange);
              }
              //fprintf(stderr, "FOUND B i=%d imax=%d best=%d query_end=%d target_end=%d\n", i, imax-zero, best-zero, *query_end, *target_end);
          }
      }
      if(query->max_aln_score - query->max_edit_score < imax) { // overflow
          if(NULL != overflow) {
              *overflow = 1;
              return tmap_vsw16_min_value;
          }
          // When overflow is going to happen, subtract tmap_vsw16_mid_value from all scores. This heuristic
          // may miss the best alignment, but in practice this should happen very rarely.
          sum += tmap_vsw16_mid_value; 
          if(query->min_aln_score + tmap_vsw16_mid_value < gmax) gmax -= tmap_vsw16_mid_value;
          else gmax = query->min_aln_score;
          for(j = 0; TMAP_VSW_LIKELY(j < slen); ++j) {
              h = __tmap_vsw16_mm_subs_epi16(__tmap_vsw_mm_load_si128(H1 + j), reduce_mm);
              __tmap_vsw_mm_store_si128(H1 + j, h);
              e = __tmap_vsw16_mm_subs_epi16(__tmap_vsw_mm_load_si128(E + j), reduce_mm);
              __tmap_vsw_mm_store_si128(E + j, e);
          }
      }
      // check for underflow
      __tmap_vsw16_min(imin, min); // imin is the minimum number in min
      if(imin < tmap_vsw16_min_value - query->min_edit_score) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }
      S = H1; H1 = H0; H0 = S; // swap H0 and H1
  }
  if(tmap_vsw16_min_value == best) {
      (*query_end) = (*target_end) = -1;
      return best;
  }
  return best + sum - zero;
}
