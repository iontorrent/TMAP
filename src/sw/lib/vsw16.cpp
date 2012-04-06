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
#include "vsw.h"
#include "vsw16.h"

// TODO: remove query profile so it does not need to be filled in each time

vsw16_query_t *
vsw16_query_init(vsw16_query_t *prev, const uint8_t *query, int32_t qlen, 
                 int32_t query_start_clip, int32_t query_end_clip,
                 vsw_opt_t *opt)
{
  int32_t qlen_mem, slen, a;
  int32_t min_gap_score, min_mm_score;
  vsw16_int_t *t;

  // get the number of stripes
  slen = __vsw16_calc_slen(qlen); 

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
      qlen_mem = __vsw_calc_qlen_mem(slen);
      //prev = memalign(16, sizeof(vsw16_query_t) + 15 + qlen_mem * (VSW_ALPHABET_SIZE + 3), "prev"); // add three for H0, H1, and E
      prev = (vsw16_query_t*)malloc(sizeof(vsw16_query_t) + 15 + qlen_mem * (VSW_ALPHABET_SIZE + 3)); // add three for H0, H1, and E
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
  prev->max_aln_score = vsw16_max_value - prev->max_edit_score; // so it doesn't overflow
  // normalize
  prev->zero_aln_score = vsw16_min_value - prev->min_aln_score- prev->zero_aln_score; 
  prev->min_aln_score = vsw16_min_value - prev->min_aln_score; 
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
  prev->query_profile = (__m128i*)__vsw_16((size_t)prev + sizeof(vsw16_query_t)); // skip over the struct variables, align memory 
  prev->H0 = prev->query_profile + (prev->slen * VSW_ALPHABET_SIZE); // skip over the query profile
  prev->H1 = prev->H0 + prev->slen; // skip over H0
  prev->E = prev->H1 + prev->slen; // skip over H1

  // create the query profile
  t = (vsw16_int_t*)prev->query_profile;
  for(a = 0; a < VSW_ALPHABET_SIZE; a++) {
      int32_t i, k;
      for(i = 0; i < prev->slen; ++i) { // for each stripe
          // fill in this stripe
          for(k = i; k < prev->slen << vsw16_values_per_128_bits_log2; k += prev->slen) { //  for q_{i+1}, q_{2i+1}, ..., q_{2s+1}
              // NB: pad with zeros
              *t++ = ((k >= qlen) ? prev->min_edit_score : ((a == query[k]) ? opt->score_match : -opt->pen_mm)); 
          }
      }
  }
  //vsw16_query_print_query_profile(prev);
  return prev;
}

void
vsw16_query_destroy(vsw16_query_t *vsw)
{
  // all memory was allocated as one block
  free(vsw);
}


static inline int32_t
vsw16_sse2_dir_cmp(int32_t cur_score, int32_t cur_qe, int32_t cur_te,
               int32_t next_score, int32_t next_qe, int32_t next_te,
               int32_t dir)
{
  if(next_score < cur_score) return 0;
  if(cur_score < next_score) return 1;
  if(0 == dir) {
      if(next_qe < cur_qe) return 1;
      if(next_qe == cur_qe && next_te < cur_te) return 1;
  }
  else {
      if(cur_qe < next_qe) return 1;
      if(cur_qe == next_qe && cur_te < next_te) return 1;
  }
  return 0;
}

int32_t
vsw16_sse2_forward(vsw16_query_t *query, const uint8_t *target, int32_t tlen, 
                   int32_t query_start_clip, int32_t query_end_clip,
                   vsw_opt_t *opt, int16_t *query_end, int16_t *target_end,
                   int32_t direction, int32_t *overflow, int32_t *n_best, int32_t score_thr)
{
  int32_t slen, i, j, k, sum = 0;
#ifdef VSW_DEBUG
  int32_t l;
#endif
  uint16_t cmp;
  vsw16_int_t gmax, best;
  vsw16_int_t zero, imin = 0, imax = 0;
  __m128i zero_mm, zero_start_mm, negative_infinity_mm, positive_infinity_mm, reduce_mm;
  __m128i pen_gapoe, pen_gape, *H0, *H1, *E;

  // initialization
  // normalize these
  zero = query->zero_aln_score; // where the normalized zero alignment score occurs
  negative_infinity_mm = __vsw16_mm_set1_epi16(query->min_aln_score); // the minimum possible value
  positive_infinity_mm = __vsw16_mm_set1_epi16(query->max_aln_score); // the minimum possible value
  score_thr += zero; // for the scoring threshold
  // these are not normalized
  pen_gapoe = __vsw16_mm_set1_epi16(opt->pen_gapo + opt->pen_gape); // gap open penalty
  pen_gape = __vsw16_mm_set1_epi16(opt->pen_gape); // gap extend penalty
  zero_mm = __vsw16_mm_set1_epi16(zero); // where the normalized zero alignment score occurs
  zero_start_mm = __vsw16_mm_set1_epi16(zero); // where the normalized zero alignment score occurs
  gmax = vsw16_min_value;
  best = vsw16_min_value;
  if(NULL != n_best) (*n_best) = 0;
  reduce_mm = __vsw16_mm_set1_epi16(vsw16_mid_value);
  // vectors
  H0 = query->H0; 
  H1 = query->H1; 
  E = query->E; 
  slen = query->slen;
  if(NULL != overflow) *overflow = 0;
#ifdef VSW_DEBUG
  fprintf(stderr, "qlen=%d tlen=%d\n", query->qlen, tlen);
  fprintf(stderr, "query_start_clip=%d query_end_clip=%d\n",
          query_start_clip, query_end_clip);
#endif

  // select query end only
  // check stripe #: __vsw16_query_index_to_stripe_number(query->qlen, slen) 
  // check byte # __vsw16_query_index_to_byte_number(query->qlen, slen)

  if(0 == query_start_clip) {
      __m128i f;
      // Set start
      for(j = 0; j < slen; j++) {
          // initialize E to negative infinity 
          __vsw_mm_store_si128(E + j, negative_infinity_mm); // E(0,j)
      }
      // NB: setting the start == 0 will be done later
      // set the leading insertions 
      f = __vsw16_mm_set1_epi16(query->min_aln_score);
      f = __vsw16_mm_insert_epi16(f, zero, 0); 
      for(j = 0; j < vsw16_values_per_128_bits; ++j) {
          f = __vsw16_mm_insert(f , -opt->pen_gapo + -opt->pen_gape + (-opt->pen_gape * slen * j) + zero, j);
      }
      for(j = 0; j < slen; j++) {
          f = __vsw16_mm_subs_epi16(f, pen_gape); // f=F(i,j)-pen_gape
          f = __vsw16_mm_max_epi16(f, negative_infinity_mm); // bound with -inf
          __vsw_mm_store_si128(H0 + ((slen + j - 1) % slen), f); // H(0,j)
      }
  }
  else {
      // Initialize all zeros
      for(i = 0; i < slen; ++i) {
          __vsw_mm_store_si128(E + i, zero_mm);
          __vsw_mm_store_si128(H0 + i, zero_mm);
      }
  }

  // the core loop
  for(i = 0, sum = 0; i < tlen; i++) { // for each base in the target
      __m128i e, h, f, g, max, min, *S;

      max = __vsw16_mm_set1_epi16(query->min_aln_score); // max is negative infinity
      min = __vsw16_mm_set1_epi16(query->max_aln_score); // min is positive infinity
      S = query->query_profile + target[i] * slen; // s is the 1st score vector

      // load H(i-1,-1)
      h = __vsw_mm_load_si128(H0 + slen - 1); // the last stripe, which holds the j-1 
      if(UNLIKELY(0 < i)) { // only if we have previous results
          h = __vsw_mm_slli_si128(h, vsw16_shift_bytes); // shift left since x64 is little endian
      }
      h = __vsw16_mm_insert_epi16(h, zero, 0);

      // e does not need to be set

      // set F to -inf
      f = __vsw16_mm_set1_epi16(query->min_aln_score);
      // leading insertions
      g = __vsw16_mm_set1_epi16(query->min_aln_score);
      if(0 == query_start_clip) { 
          g = __vsw16_mm_set1_epi16(query->min_aln_score);
          g = __vsw16_mm_insert_epi16(g, -opt->pen_gapo + opt->pen_gape + zero, 0);
          for(j = 1; LIKELY(j < vsw16_values_per_128_bits); j++) {
              // NB: do not include the gap extend, that will be added below
              // BEFORE considering H
              g = __vsw16_mm_insert(g, -opt->pen_gapo + (-opt->pen_gape * (j * slen - 1)) + zero, j);
          }
      }

      for(j = 0; LIKELY(j < slen); ++j) { // for each stripe in the query
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
              h = __vsw16_mm_max_epi16(h, zero_start_mm);  
          }
          else {
              g = __vsw16_mm_subs_epi16(g, pen_gape); // leading insertion
#ifdef VSW_DEBUG 
              /*
                 fprintf(stderr, "h i=%d j=%d", i, j);
                 for(k = 0; k < vsw16_values_per_128_bits; k++) { // for each start position in the stripe
                 fprintf(stderr, " (%d,%d)", 
                 ((vsw16_int_t*)(&h))[k] - zero,
                 ((vsw16_int_t*)(&g))[k] - zero);
                 }
                 fprintf(stderr, "\n");
                 */
#endif
              h = __vsw16_mm_max_epi16(h, g);
          }
          // compute H(i,j); 
#ifdef VSW_DEBUG 
          /*
             __m128i s = __vsw_mm_load_si128(S + j);
             fprintf(stderr, "s i=%d j=%d", i, j);
             for(k = 0; k < vsw16_values_per_128_bits; k++) { // for each start position in the stripe
             fprintf(stderr, " %4d", ((vsw16_int_t*)(&s))[k]);
             }
             fprintf(stderr, "\n");
             fprintf(stderr, "h i=%d j=%d", i, j);
             for(k = 0; k < vsw16_values_per_128_bits; k++) { // for each start position in the stripe
             fprintf(stderr, " %4d", ((vsw16_int_t*)(&h))[k] - zero);
             }
             fprintf(stderr, "\n");
             */
#endif
          h = __vsw16_mm_adds_epi16(h, __vsw_mm_load_si128(S + j)); // h=H(i-1,j-1)+S(i,j)
          e = __vsw_mm_load_si128(E + j); // e=E(i,j)
          h = __vsw16_mm_max_epi16(h, e); // h=H(i,j) = max{E(i,j), H(i-1,j-1)+S(i,j)}
          h = __vsw16_mm_max_epi16(h, f); // h=H(i,j) = max{max{E(i,j), H(i-1,j-1)+S(i,j)}, F(i,j)}
          h = __vsw16_mm_max_epi16(h, negative_infinity_mm); // bound with -inf
          max = __vsw16_mm_max_epi16(max, h); // save the max values in this stripe versus the last
          min = __vsw16_mm_min_epi16(min, h); // save the max values in this stripe versus the last
          __vsw_mm_store_si128(H1 + j, h); // save h to H(i,j)
          // next, compute E(i+1,j)
          h = __vsw16_mm_subs_epi16(h, pen_gapoe); // h=H(i,j)-pen_gapoe
          e = __vsw16_mm_subs_epi16(e, pen_gape); // e=E(i,j)-pen_gape
          e = __vsw16_mm_max_epi16(e, h); // e=E(i+1,j) = max{E(i,j)-pen_gape, H(i,j)-pen_gapoe}
          e = __vsw16_mm_max_epi16(e, negative_infinity_mm); // bound with -inf
          __vsw_mm_store_si128(E + j, e); // save e to E(i+1,j)
          // now compute F(i,j+1)
          //h = __vsw16_mm_subs_epi16(h, pen_gapoe); // h=H(i,j)-pen_gapoe
          f = __vsw16_mm_subs_epi16(f, pen_gape); // f=F(i,j)-pen_gape
          f = __vsw16_mm_max_epi16(f, h); // f=F(i,j+1) = max{F(i,j)-pen_gape, H(i,j)-pen_gapoe}
          f = __vsw16_mm_max_epi16(f, negative_infinity_mm); // bound with -inf
          // get H(i-1,j) and prepare for the next j
          h = __vsw_mm_load_si128(H0 + j); // h=H(i-1,j)
      }

      // NB: we do not need to set E(i,j) as we disallow adjacent insertion and then deletion
      // iterate through each value stored in a stripe
      f = __vsw16_mm_set1_epi16(query->min_aln_score);
      // we require at most 'vsw16_values_per_128_bits' iterations to guarantee F propagation
      for(k = 0; LIKELY(k < vsw16_values_per_128_bits); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
          f = __vsw_mm_slli_si128(f, vsw16_shift_bytes); // since x86 is little endian
          f = __vsw16_mm_insert_epi16(f, query->min_aln_score, 0); // set F(i-1,-1)[0] as negative infinity (normalized)
          for(j = 0; LIKELY(j < slen); ++j) { 
              h = __vsw_mm_load_si128(H1 + j); // h=H(i,j)
              h = __vsw16_mm_max_epi16(h, f); // h=H(i,j) = max{H(i,j), F(i,j)}
              h = __vsw16_mm_max_epi16(h, negative_infinity_mm); // bound with -inf
              __vsw_mm_store_si128(H1 + j, h); // save h to H(i,j) 
              h = __vsw16_mm_subs_epi16(h, pen_gapoe); // h=H(i,j)-gapo
              f = __vsw16_mm_subs_epi16(f, pen_gape); // f=F(i,j)-pen_gape
              f = __vsw16_mm_max_epi16(f, negative_infinity_mm); // bound with -inf
              // check to see if h could have been updated by f?
              // NB: the comparison below will have some false positives, in
              // other words h == f a priori, but this is rare.
              cmp = __vsw16_mm_movemask_epi16(__vsw16_mm_cmpgt_epi16(f, h));
              if (LIKELY(cmp = 0x0000)) goto end_loop; // ACK: goto statement
              f = __vsw16_mm_max_epi16(f, h); // f=F(i,j+1) = max{F(i,j)-pen_gape, H(i,j)-pen_gapoe}
          }
      }
end_loop:
#ifdef VSW_DEBUG 
      fprintf(stderr, "H1 i=%d target[i]=%d", i, target[i]);
      for(k = l = 0; k < vsw16_values_per_128_bits; k++) { // for each start position in the stripe
          for(j = 0; j < slen; j++, l++) {
              fprintf(stderr, " %d:%d", l, ((vsw16_int_t*)(H1 + j))[k] - zero);
          }
      }
      fprintf(stderr, "\n");
#endif
      __vsw16_max(imax, max); // imax is the maximum number in max
      if(query->max_aln_score - query->max_edit_score < imax) { // overflow
          if(NULL != overflow) *overflow = 1;
          return vsw16_min_value;
      }
      if(imax > gmax) { 
          gmax = imax; // global maximum score 
      }
      if(score_thr <= imax && best <= imax) { // potential best score
          vsw16_int_t *t;
          //int32_t save_score = 0;
          //if(imax > best || (1 == direction && imax == best)) save_score = 1;
          if(query_end_clip == 0) { // check the last
              j = (query->qlen-1) % slen; // stripe
              k = (query->qlen-1) / slen; // byte
              t = (vsw16_int_t*)(H1 + j);
#ifdef VSW_DEBUG 
              fprintf(stderr, "j=%d k=%d slen=%d qlen=%d best=%d pot=%d\n", j, k, slen, query->qlen, best-zero, t[k]-zero);
#endif
              if(NULL != n_best) {
                  if(best == (int32_t)t[k]) { // duplicate best score
                      (*n_best)++;
                  }
                  else if(best < (int32_t)t[k]) {
                      (*n_best) = 1;
                  }
              }
              if(1 == vsw16_sse2_dir_cmp(best, (*query_end), (*target_end), (int32_t)t[k], query->qlen-1, i, direction)) {
                  (*query_end) = query->qlen-1;
                  (*target_end) = i;
                  best = t[k];
#ifdef VSW_DEBUG 
                  fprintf(stderr, "FOUND A i=%d query->qlen-1=%d query_end=%d target_end=%d best=%d\n", i, query->qlen-1, *query_end, *target_end, best-zero);
#endif
              }
          }
          else { // check all
              int32_t found_best = 0;
              t = (vsw16_int_t*)H1;
              if(best < imax || (best == imax && 1 == direction)) {
                  *query_end = -1;
              }
              for(j = 0; LIKELY(j < slen); ++j) { // for each stripe in the query
                  for(k = 0; k < vsw16_values_per_128_bits; k++, t++) { // for each cell in the stripe
                      if(NULL != n_best) {
                          if(best == (int32_t)*t) { // duplicate best score
                              (*n_best)++;
                          }
                          else if(best < (int32_t)*t) {
                              (*n_best) = 1;
                          }
                      }
                      int32_t cur_qe = j + ((k & (vsw16_values_per_128_bits-1)) * slen);
                      if(1 == vsw16_sse2_dir_cmp(best, (*query_end), (*target_end), (int32_t)*t, cur_qe, i, direction)) {
                          found_best = 1;
                          best = *t;
                          (*query_end) = cur_qe;
                          (*target_end) = i;
#ifdef VSW_DEBUG 
                          fprintf(stderr, "FOUND B j=%d k=%d query_end=%d best=%d\n", j, k, *query_end, best-zero);
#endif
                      }
                  }
              }
              if(1 == found_best && -1 == *query_end) {
                  if(-1 == *query_end) {
                      fprintf(stderr, "bug encountered\n");
                      exit(1);
                  }
              }
              //fprintf(stderr, "FOUND B i=%d imax=%d best=%d query_end=%d target_end=%d\n", i, imax-zero, best-zero, *query_end, *target_end);
          }
      }
      if(query->max_aln_score - query->max_edit_score < imax) { // overflow
          if(NULL != overflow) {
              *overflow = 1;
              return vsw16_min_value;
          }
          // When overflow is going to happen, subtract vsw16_mid_value from all scores. This heuristic
          // may miss the best alignment, but in practice this should happen very rarely.
          sum += vsw16_mid_value; 
          if(query->min_aln_score + vsw16_mid_value < gmax) gmax -= vsw16_mid_value;
          else gmax = query->min_aln_score;
          for(j = 0; LIKELY(j < slen); ++j) {
              h = __vsw16_mm_subs_epi16(__vsw_mm_load_si128(H1 + j), reduce_mm);
              __vsw_mm_store_si128(H1 + j, h);
              e = __vsw16_mm_subs_epi16(__vsw_mm_load_si128(E + j), reduce_mm);
              __vsw_mm_store_si128(E + j, e);
          }
      }
      // check for underflow
      __vsw16_min(imin, min); // imin is the minimum number in min
      if(imin < vsw16_min_value - query->min_edit_score) {
          fprintf(stderr, "bug encountered\n");
          exit(1);
      }
      S = H1; H1 = H0; H0 = S; // swap H0 and H1
  }
  if(vsw16_min_value == best) {
      (*query_end) = (*target_end) = -1;
      return best;
  }
  return best + sum - zero;
}
