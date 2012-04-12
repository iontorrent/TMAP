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
#include "tmap_vsw_definitions.h"
#include "tmap_vsw.h"

tmap_vsw_t*
tmap_vsw_init(const uint8_t *query, int32_t qlen,
                    int32_t query_start_clip, int32_t query_end_clip,
                    int32_t type,
                    tmap_vsw_opt_t *opt)
{
  tmap_vsw_t *vsw = NULL;
  vsw = tmap_calloc(1, sizeof(tmap_vsw_t), "vsw");
  vsw->type = type;
  vsw->query_start_clip = query_start_clip;
  vsw->query_end_clip = query_end_clip;
  vsw->opt = opt;
  vsw->algorithm = tmap_vsw_wrapper_init(type);
  vsw->algorithm_default = tmap_vsw_wrapper_init(1);
  return vsw;
}

void
tmap_vsw_destroy(tmap_vsw_t *vsw)
{
  if(NULL == vsw) return;
  tmap_vsw_wrapper_destroy(vsw->algorithm);
  tmap_vsw_wrapper_destroy(vsw->algorithm_default);
  free(vsw);
}

#ifdef TMAP_VSW_DEBUG_CMP
static void
tmap_vsw_process_compare(tmap_vsw_t *vsw,
                         const uint8_t *query, int32_t qlen,
                         uint8_t *target, int32_t tlen, 
                         tmap_vsw_result_t *result,
                         int32_t *overflow, int32_t score_thr, int32_t is_rev)
{
  int32_t query_end_type0, target_end_type0, n_best_type0, score_type0;
  int32_t query_end_type1, target_end_type1, n_best_type1, score_type1;

  query_end_type0 = target_end_type0 = n_best_type0 = score_type0 = 0;
  query_end_type1 = target_end_type1 = n_best_type1 = score_type1 = 0;
  
  // baseline
  tmap_vsw_wrapper_process(vsw->algorithm_default,
                           target, tlen, 
                           query, qlen, 
                           vsw->opt->score_match,
                           -vsw->opt->pen_mm,
                           -vsw->opt->pen_gapo,
                           -vsw->opt->pen_gape,
                           is_rev, 
                           vsw->query_start_clip, vsw->query_end_clip, 
                           &score_type0, &target_end_type0, &query_end_type0, &n_best_type0);
  // current
  tmap_vsw_wrapper_process(vsw->algorithm,
                           target, tlen, 
                           query, qlen, 
                           vsw->opt->score_match,
                           -vsw->opt->pen_mm,
                           -vsw->opt->pen_gapo,
                           -vsw->opt->pen_gape,
                           is_rev, 
                           vsw->query_start_clip, vsw->query_end_clip, 
                           &score_type1, &target_end_type1, &query_end_type1, &n_best_type1);

  if(tlen <= target_end_type0) tmap_bug();
  if(qlen <= query_end_type0) tmap_bug();
  /*
  if(tlen <= target_end_type1) tmap_bug();
  if(qlen <= query_end_type1) tmap_bug();
  */

  if(score_type0 != score_type1
     || target_end_type0 != target_end_type1
     || query_end_type0 != query_end_type1
     || n_best_type0 != n_best_type1) {
      int32_t i;
      fprintf(stderr, "in %s is_rev=%d\n", __func__, is_rev);
      fprintf(stderr, "query_start_clip=%d\n", vsw->query_start_clip);
      fprintf(stderr, "query_end_clip=%d\n", vsw->query_end_clip);
      fprintf(stderr, "qlen=%d tlen=%d\n", qlen, tlen);
      for(i=0;i<qlen;i++) {
          fputc("ACGTN"[query[i]], stderr);
      }
      fputc('\n', stderr);
      for(i=0;i<tlen;i++) {
          fputc("ACGTN"[target[i]], stderr);
      }
      fputc('\n', stderr);
      fprintf(stderr, "tlen=%d qlen=%d score=[%d,%d] target_end=[%d,%d] query_end=[%d,%d] n_best=[%d,%d]\n",
              tlen, qlen,
              score_type0, score_type1,
              target_end_type0, target_end_type1,
              query_end_type0, query_end_type1,
              n_best_type0, n_best_type1);
      do {
          tmap_sw_param_t ap;
          tmap_sw_path_t *path = NULL;
          int32_t i, matrix[25], n_cigar, score, path_len;
          uint32_t *cigar = NULL;
          ap.matrix=matrix;
          for(i=0;i<25;i++) { 
              ap.matrix[i] = -(vsw->opt)->pen_mm; 
          } 
          for(i=0;i<4;i++) { 
              ap.matrix[i*5+i] = vsw->opt->score_match; 
          } 
          ap.gap_open = vsw->opt->pen_gapo; ap.gap_ext = vsw->opt->pen_gape; 
          ap.gap_end = vsw->opt->pen_gape; 
          ap.row = 5; 
          path_len = 0;
          path = tmap_calloc(1024, sizeof(tmap_sw_path_t), "path");
          score = tmap_sw_clipping_core((uint8_t*)target, tlen, (uint8_t*)query, qlen, &ap,
                                        vsw->query_start_clip, vsw->query_end_clip, 
                                        path, &path_len, is_rev);
          // print out the path
          cigar = tmap_sw_path2cigar(path, path_len, &n_cigar);
          fprintf(stderr, "tmap_sw_clipping_core score=%d\n", score);
          for(i=0;i<n_cigar;i++) {
              fprintf(stderr, "%d%c", cigar[i]>>4, "MIDNSHP"[cigar[i]&0xf]);
          }
          fputc('\n', stderr);
          free(path);
          free(cigar);
      } while(0);
      // try the opposite direction
      // baseline
      tmap_vsw_wrapper_process(vsw->algorithm_default,
                               target, tlen, 
                               query, qlen, 
                               vsw->opt->score_match,
                               -vsw->opt->pen_mm,
                               -vsw->opt->pen_gapo,
                               -vsw->opt->pen_gape,
                               1-is_rev, 
                               vsw->query_start_clip, vsw->query_end_clip, 
                               &score_type0, &target_end_type0, &query_end_type0, &n_best_type0);
      fprintf(stderr, "baseline opposite reverse tlen=%d qlen=%d score=[%d,%d] target_end=[%d,%d] query_end=[%d,%d] n_best=[%d,%d]\n",
              tlen, qlen,
              score_type0, score_type1,
              target_end_type0, target_end_type1,
              query_end_type0, query_end_type1,
              n_best_type0, n_best_type1);
      // top coder
      for(i=0;i<tlen;i++) fputc("ACGTN"[target[i]], stderr);
      fputc('\t', stderr);
      for(i=0;i<qlen;i++) fputc("ACGTN"[query[i]], stderr);
      fputc('\t', stderr);
      fprintf(stderr, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
              vsw->query_start_clip, vsw->query_end_clip, 
              vsw->opt->score_match,
              -vsw->opt->pen_mm,
              -vsw->opt->pen_gapo,
              -vsw->opt->pen_gape,
              is_rev,
              -1, -1, -1, -1);
      fprintf(stderr, "Error: algorithms produced different results!\n");
      exit(1);
  }
}
#endif

static int32_t
tmap_vsw_process(tmap_vsw_t *vsw,
              const uint8_t *query, int32_t qlen,
              uint8_t *target, int32_t tlen, 
              tmap_vsw_result_t *result,
              int32_t *overflow, int32_t score_thr, 
              int32_t is_rev, int32_t direction)
{
  int32_t found_forward = 1, query_end, target_end, n_best, score = INT32_MIN;
#ifdef TMAP_VSW_DEBUG
  int32_t i;
#endif
  // TODO: check potential overflow
  // TODO: check that gap penalties will not result in an overflow
  // TODO: check that the max/min alignment score do not result in an overflow

#ifdef TMAP_VSW_DEBUG
  fprintf(stderr, "in %s is_rev=%d\n", __func__, is_rev);
  fprintf(stderr, "query_start_clip=%d\n", vsw->query_start_clip);
  fprintf(stderr, "query_end_clip=%d\n", vsw->query_end_clip);
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

  // update based on current problem
  query_end = target_end = n_best = 0;
  if(NULL != overflow) (*overflow) = 0;

  if(tlen <= tmap_vsw_wrapper_get_max_tlen(vsw->algorithm)
     && qlen <= tmap_vsw_wrapper_get_max_qlen(vsw->algorithm)) {
      tmap_vsw_wrapper_process(vsw->algorithm,
                               target, tlen, 
                               query, qlen, 
                               vsw->opt->score_match,
                               -vsw->opt->pen_mm,
                               -vsw->opt->pen_gapo,
                               -vsw->opt->pen_gape,
                               direction, 
                               vsw->query_start_clip, vsw->query_end_clip, 
                               &score, &target_end, &query_end, &n_best);
#ifdef TMAP_VSW_DEBUG_CMP
      tmap_vsw_process_compare(vsw,
                               query, qlen,
                               target, tlen, 
                               result,
                               overflow, score_thr, direction);
#endif
  }
  else { // try the default
      tmap_vsw_wrapper_process(vsw->algorithm_default,
                               target, tlen, 
                               query, qlen, 
                               vsw->opt->score_match,
                               -vsw->opt->pen_mm,
                               -vsw->opt->pen_gapo,
                               -vsw->opt->pen_gape,
                               direction, 
                               vsw->query_start_clip, vsw->query_end_clip, 
                               &score, &target_end, &query_end, &n_best);
  }
  if(score < score_thr || 0 == n_best) {
      query_end = target_end = -1;
      n_best = 0;
      score = INT32_MIN;
  }
  
  if(0 == is_rev) {

      result->query_end = query_end;
      result->target_end = target_end;
      result->n_best = n_best;
      result->score_fwd = score;

      // check forward results
      if(NULL != overflow && 1 == (*overflow)) {
          found_forward = 0;
      }
      else if(result->score_fwd <- score_thr) {
          found_forward = 0;
      }
      else if((result->query_end == result->query_start || result->target_end == result->target_start)
              && result->score_fwd <= 0) {
          found_forward = 0;
      }
      else if(n_best <= 0) {
          found_forward = 0;
      }
      else if(-1 == result->query_end) {
          tmap_bug();
      }

      // return if we found no legal/good forward results
      if(0 == found_forward) {
          result->query_end = result->query_start = 0;
          result->target_end = result->target_start = 0;
          result->score_fwd = result->score_rev = INT16_MIN;
          result->n_best = 0;
          return INT32_MIN;
      }

      result->query_start = result->target_start = 0;
      result->score_rev = INT16_MIN;

#ifdef TMAP_VSW_DEBUG
      fprintf(stderr, "result->score_fwd=%d result->score_rev=%d\n",
              result->score_fwd, result->score_rev);
      fprintf(stderr, "{?-%d] {?-%d}\n",
              result->query_end,
              result->target_end);
#endif

      return result->score_fwd;
  }
  else {

      result->query_start = qlen - query_end - 1;
      result->target_start = tlen - target_end - 1;
      result->n_best = n_best;
      result->score_rev = score;

#ifdef TMAP_VSW_DEBUG
      fprintf(stderr, "is_rev=%d result->score_fwd=%d result->score_rev=%d\n",
              is_rev, result->score_fwd, result->score_rev);
#endif

      // check reverse results
      if(NULL != overflow && 1 == (*overflow)) {
          result->query_end = result->query_start = 0;
          result->target_end = result->target_start = 0;
          result->score_fwd = result->score_rev = INT16_MIN;
          result->n_best = 0;
          return INT32_MIN; 
      }
      else if(result->score_fwd != result->score_rev) { // something went wrong... FIXME
          // use the default
          if(0 != vsw->type) {
              tmap_vsw_wrapper_process(vsw->algorithm_default,
                                       target, tlen, 
                                       query, qlen, 
                                       vsw->opt->score_match,
                                       -vsw->opt->pen_mm,
                                       -vsw->opt->pen_gapo,
                                       -vsw->opt->pen_gape,
                                       direction, 
                                       vsw->query_start_clip, vsw->query_end_clip, 
                                       &score, &target_end, &query_end, &n_best);
              result->query_start = qlen - query_end - 1;
              result->target_start = tlen - target_end - 1;
              result->n_best = n_best;
              result->score_rev = score;
              if(result->score_fwd != result->score_rev) { // something went wrong... again...
                  // ignore
                  result->query_end = result->query_start = 0;
                  result->target_end = result->target_start = 0;
                  result->score_fwd = result->score_rev = INT16_MIN;
                  result->n_best = 0;
                  return INT32_MIN;
              }
              return score;
          }
          // Bug!
          fprintf(stderr, "{%d-%d} {%d-%d}\n",
                  result->query_start, result->query_end,
                  result->target_start, result->target_end);
          fprintf(stderr, "result->score_fwd=%d result->score_rev=%d\n",
                  result->score_fwd, result->score_rev);
          fprintf(stderr, "score=%d\n", score);
          /*
          tmap_sw_param_t ap;
          int32_t i, matrix[25];
          ap.matrix=matrix;
          for(i=0;i<25;i++) { 
              ap.matrix[i] = -(vsw->opt)->pen_mm; 
          } 
          for(i=0;i<4;i++) { 
              ap.matrix[i*5+i] = vsw->opt->score_match; 
          } 
          ap.gap_open = vsw->opt->pen_gapo; ap.gap_ext = vsw->opt->pen_gape; 
          ap.gap_end = vsw->opt->pen_gape; 
          ap.row = 5; 
          score = tmap_sw_clipping_core((uint8_t*)target, tlen, (uint8_t*)query, qlen, &ap,
                                        vsw->query_start_clip, vsw->query_end_clip, 
                                        NULL, NULL, direction);
          fprintf(stderr, "score=%d\n", score);
          */
          tmap_bug();
      }

      return result->score_fwd;
  }
}

int32_t
tmap_vsw_process_fwd(tmap_vsw_t *vsw,
              const uint8_t *query, int32_t qlen,
              uint8_t *target, int32_t tlen, 
              tmap_vsw_result_t *result,
              int32_t *overflow, int32_t score_thr, int32_t direction)
{
  return tmap_vsw_process(vsw, query, qlen, target, tlen, result, overflow, score_thr, 0, direction);
}

int32_t
tmap_vsw_process_rev(tmap_vsw_t *vsw,
              const uint8_t *query, int32_t qlen,
              uint8_t *target, int32_t tlen, 
              tmap_vsw_result_t *result,
              int32_t *overflow, int32_t score_thr, int32_t direction)
{
  return tmap_vsw_process(vsw, query, qlen, target, tlen, result, overflow, score_thr, 1, direction);
}
