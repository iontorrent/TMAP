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
#include <config.h>
#include "../util/tmap_alloc.h"
#include "../util/tmap_error.h"
#include "../util/tmap_definitions.h"
#include "tmap_sw.h"
#include "tmap_vsw16.h"
#include "tmap_vsw.h"

tmap_vsw_opt_t*
tmap_vsw_opt_init(int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape, int32_t score_thr)
{
  tmap_vsw_opt_t *opt;

  opt = tmap_malloc(sizeof(tmap_vsw_opt_t), "opt");
  opt->score_match = score_match; 
  opt->pen_mm = pen_mm;
  opt->pen_gapo = pen_gapo;
  opt->pen_gape = pen_gape;
  opt->score_thres = score_thr;

  return opt;
}

void
tmap_vsw_opt_destroy(tmap_vsw_opt_t *opt)
{
  free(opt);
}

tmap_vsw_query_t*
tmap_vsw_query_init(const uint8_t *query, int32_t qlen,
                    int32_t tlen,
                    int32_t query_start_clip, int32_t query_end_clip,
                    tmap_vsw_opt_t *opt)
{
  tmap_vsw_query_t *vsw_query;
  vsw_query = tmap_calloc(1, sizeof(tmap_vsw_query_t), "query");
  vsw_query->query16 = tmap_vsw16_query_init(vsw_query->query16, query, qlen, query_start_clip, query_end_clip, opt);
  return vsw_query;
}

void
tmap_vsw_query_destroy(tmap_vsw_query_t *query)
{
  if(NULL == query) return;
  if(NULL != query->query16) tmap_vsw16_query_destroy(query->query16);
  free(query);
}

int32_t
tmap_vsw_sse2(tmap_vsw_query_t *vsw_query,
              const uint8_t *query, int32_t qlen,
              uint8_t *target, int32_t tlen, 
              int32_t query_start_clip, int32_t query_end_clip,
              tmap_vsw_opt_t *opt, 
              int16_t *score_fwd, int16_t *score_rev,
              int16_t *query_start, int16_t *query_end,
              int16_t *target_start, int16_t *target_end,
              int32_t *overflow, int32_t *n_best, int32_t score_thr, int32_t is_rev)
{
  uint8_t *tmp_target;
  int32_t tmp_tlen, found_forward = 1;
#ifdef TMAP_VSW_DEBUG
  int32_t i;
#endif
  // TODO: check potential overflow
  // TODO: check that gap penalties will not result in an overflow
  // TODO: check that the max/min alignment score do not result in an overflow


#ifdef TMAP_VSW_DEBUG
  fprintf(stderr, "in %s is_rev=%d\n", __func__, is_rev);
  fprintf(stderr, "query_start_clip=%d\n", query_start_clip);
  fprintf(stderr, "query_end_clip=%d\n", query_end_clip);
  for(i=0;i<qlen;i++) {
      fputc("ACGTN"[query[i]], stderr);
  }
  fputc('\n', stderr);
  for(i=0;i<tlen;i++) {
      fputc("ACGTN"[target[i]], stderr);
  }
  fputc('\n', stderr);
#endif

  if(0 == is_rev) {
      // init forward
      vsw_query->query16 = tmap_vsw16_query_init(vsw_query->query16, query, qlen, query_start_clip, query_end_clip, opt);

      // run forward
      (*score_fwd) = tmap_vsw16_sse2_forward(vsw_query->query16, target, tlen,
                                                  query_start_clip, query_end_clip,
                                                  opt, query_end, target_end,
                                                  0, overflow, n_best, score_thr);

      // check forward results
      if(NULL != overflow && 1 == *overflow) {
          found_forward = 0;
      }
      else if((*score_fwd) < score_thr) {
          found_forward = 0;
      }
      else if(((*query_end) == (*query_start) || (*target_end) == (*target_start))
              && (*score_fwd) <= 0) {
          found_forward = 0;
      }
      else if(-1 == (*query_end)) {
          tmap_bug();
      }

      // return if we found no legal/good forward results
      if(0 == found_forward) {
          (*query_end) = (*query_start) = 0;
          (*target_end) = (*target_start) = 0;
          (*score_fwd) = (*score_rev) = INT16_MIN;
          if(NULL != n_best) (*n_best) = 0;
          return INT32_MIN;
      }

      (*query_start) = (*target_start) = 0;
      (*score_rev) = INT16_MIN;

#ifdef TMAP_VSW_DEBUG
      fprintf(stderr, "(*score_fwd)=%d (*score_rev)=%d\n",
              (*score_fwd), (*score_rev));
      fprintf(stderr, "{?-%d] {?-%d}\n",
              (*query_end),
              (*target_end));
#endif

      return (*score_fwd);
  }
  else {

      // reverse compliment
      tmp_target = target;
      tmp_tlen = tlen;
      tmap_reverse_compliment_int(tmp_target, tmp_tlen);

#ifdef TMAP_VSW_DEBUG
      fprintf(stderr, "{?-%d] {?-%d}\n",
              (*query_end),
              (*target_end));
      fprintf(stderr, "qlen=%d tlen=%d\n", qlen, tlen);
      for(i=0;i<qlen;i++) {
          fputc("ACGTN"[query[i]], stderr);
      }
      fputc('\n', stderr);
      for(i=0;i<tlen;i++) {
          fputc("ACGTN"[target[i]], stderr);
      }
      fputc('\n', stderr);
#endif

      // init reverse
      // NB: reverse the clipping
      vsw_query->query16 = tmap_vsw16_query_init(vsw_query->query16, query, qlen, query_end_clip, query_start_clip, opt);

      // run the reverse
      // NB: we do not start clip since we know here the alignment should begin on
      // the reverse
      (*score_rev) = tmap_vsw16_sse2_forward(vsw_query->query16, target, tlen,
                                                  0/*query_end_clip*/, query_start_clip,
                                                  opt, query_start, target_start,
                                                  1, overflow, n_best, score_thr);

#ifdef TMAP_VSW_DEBUG
      fprintf(stderr, "is_rev=%d (*score_fwd)=%d (*score_rev)=%d\n",
              is_rev, (*score_fwd), (*score_rev));
#endif

      // check reverse results
      if(1 == *overflow) {
          (*query_end) = (*query_start) = 0;
          (*target_end) = (*target_start) = 0;
          (*score_fwd) = (*score_rev) = INT16_MIN;
          if(NULL != n_best) (*n_best) = 0;
          return INT32_MIN; 
      }
      else if((*score_fwd) != (*score_rev)) {
          fprintf(stderr, "{%d-%d} {%d-%d}\n",
                  (*query_start), (*query_end),
                  (*target_start), (*target_end));
          fprintf(stderr, "(*score_fwd)=%d (*score_rev)=%d\n",
                  (*score_fwd), (*score_rev));
          tmap_bug();
      }

#ifdef TMAP_VSW_DEBUG
      fprintf(stderr, "{%d-%d} {%d-%d}\n",
              (*query_start), (*query_end),
              (*target_start), (*target_end));
#endif

      // adjust query_start and query_end
      (*query_start) = qlen - (*query_start) - 1;
      (*target_start) = tlen - (*target_start) - 1;

#ifdef TMAP_VSW_DEBUG
      fprintf(stderr, "{%d-%d} {%d-%d}\n",
              (*query_start), (*query_end),
              (*target_start), (*target_end));
#endif
      // reverse compliment (is this necessary?)
      tmap_reverse_compliment_int(tmp_target, tmp_tlen);

      // return
      return (*score_fwd);
  }
}
