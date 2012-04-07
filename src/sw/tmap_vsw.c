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
#ifdef TMAP_VSW_DEBUG_CMP
  // TODO
#endif
  return vsw;
}

void
tmap_vsw_destroy(tmap_vsw_t *vsw)
{
  if(NULL == vsw) return;
  tmap_vsw_wrapper_destroy(vsw->algorithm);
  tmap_vsw_wrapper_destroy(vsw->algorithm_default);
#ifdef TMAP_VSW_DEBUG_CMP
  // TODO
#endif
  free(vsw);
}

/*
void
tmap_vsw_get_max(tmap_vsw_t *vsw, int32_t *max_qlen, int32_t *max_tlen)
{
}
*/

int32_t
tmap_vsw_update(tmap_vsw_t *vsw, const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen)
{
  // TODO
  /*
  int32_t max_qlen=0, max_tlen=0;
#ifdef TMAP_VSW_DEBUG_CMP
  vsw->s0 = tmap_vsw_data_update_s0(vsw->s0, query, qlen, target, tlen);
  vsw->s1 = tmap_vsw_data_update_s1(vsw->s1, query, qlen, target, tlen);
  vsw->s3 = tmap_vsw_data_update_s3(vsw->s3, query, qlen, target, tlen);
#endif
  tmap_vsw_get_max(vsw, &max_qlen, &max_tlen); // get the maximum qlen/tlen supported
  if(1 == vsw->use_default && (max_qlen < qlen || max_tlen < tlen)) { // not suported
      // update the default only
      vsw->default_s = tmap_vsw_data_update_s0(vsw->default_s, query, qlen, target, tlen);
      return 0;
  }
  else {
      // update the current algorithm
      switch(vsw->type) {
        case TMAP_VSW_TYPE_S0:
          vsw->data.s0 = tmap_vsw_data_update_s0(vsw->data.s0, query, qlen, target, tlen);
          break;
        case TMAP_VSW_TYPE_S1:
          vsw->data.s1 = tmap_vsw_data_update_s1(vsw->data.s1, query, qlen, target, tlen);
          break;
        case TMAP_VSW_TYPE_S3:
          vsw->data.s3 = tmap_vsw_data_update_s3(vsw->data.s3, query, qlen, target, tlen);
          break;
        default:
          tmap_bug();
          break;
      }
      return 1;
  }
  */
  return 1;
}
  
#ifdef TMAP_VSW_DEBUG_CMP
// TODO
/*
int32_t
tmap_vsw_process_compare2(tmap_vsw_t *vsw,
                         const uint8_t *query, int32_t qlen,
                         uint8_t *target, int32_t tlen, 
                         tmap_vsw_result_t *result,
                         int32_t *overflow, int32_t score_thr, int32_t is_rev,
                         int32_t *query_end, int32_t *target_end, int32_t *n_best, int32_t *score,
                         int32_t type)
{
  switch(type) {
    case TMAP_VSW_TYPE_S0:
      (*score) = tmap_vsw_process_s0(vsw->s0, 
                                     query, qlen, target, tlen, 
                                     vsw->query_start_clip, vsw->query_end_clip, vsw->opt,
                                     is_rev, score_thr, 
                                     query_end, target_end, n_best, overflow);
      break;
    case TMAP_VSW_TYPE_S1:
      (*score) = tmap_vsw_process_s1(vsw->s1, 
                                     query, qlen, target, tlen, 
                                     vsw->query_start_clip, vsw->query_end_clip, vsw->opt,
                                     is_rev, score_thr, 
                                     query_end, target_end, n_best, overflow);
      break;
    case TMAP_VSW_TYPE_S3:
      (*score) = tmap_vsw_process_s3(vsw->s3, 
                                     query, qlen, target, tlen, 
                                     vsw->query_start_clip, vsw->query_end_clip, vsw->opt,
                                     is_rev, score_thr, 
                                     query_end, target_end, n_best, overflow);
      break;
    default:
      return 0;
  }
  return 1;
}

static void
tmap_vsw_process_compare(tmap_vsw_t *vsw,
                         const uint8_t *query, int32_t qlen,
                         uint8_t *target, int32_t tlen, 
                         tmap_vsw_result_t *result,
                         int32_t *overflow, int32_t score_thr, int32_t is_rev,
                         int32_t from, int32_t to)
{
  int32_t query_end_type0, target_end_type0, n_best_type0, score_type0;
  int32_t query_end_type1, target_end_type1, n_best_type1, score_type1;
  int32_t baseline = from;

  // base line
  tmap_vsw_process_compare2(vsw, query, qlen, target, tlen, result,
                            overflow, score_thr, is_rev,
                            &query_end_type0, &target_end_type0, &n_best_type0, &score_type0,
                            from);

  // go through the rest
  from++;
  while(from <= to) {
      if(1 == tmap_vsw_process_compare2(vsw, query, qlen, target, tlen, result,
                                        overflow, score_thr, is_rev,
                                        &query_end_type1, &target_end_type1, &n_best_type1, &score_type1,
                                        from)) {
          if(score_type0 != score_type1 
             || query_end_type0 != query_end_type1
             || target_end_type0 != target_end_type1
             || n_best_type0 != n_best_type1) {
              int32_t i;
              fprintf(stderr, "baseline=%d current-type=%d\n",
                      baseline, from);

              fprintf(stderr, "QSC=%d QEC=%d is_rev=%d\n",
                      vsw->query_start_clip,
                      vsw->query_end_clip,
                      is_rev);
              for(i=0;i<qlen;i++) {
                  fputc("ACGTN"[query[i]], stderr);
              }
              fputc('\n', stderr);
              for(i=0;i<tlen;i++) {
                  fputc("ACGTN"[target[i]], stderr);
              }
              fputc('\n', stderr);
              fprintf(stderr, "score=[%d,%d] QE=[%d,%d] TE=[%d,%d] NB=[%d,%d]\n",
                      score_type0, score_type1,
                      query_end_type0, query_end_type1,
                      target_end_type0, target_end_type1,
                      n_best_type0, n_best_type1);
              tmap_bug();
          }
      }
      from++;
  }
}
*/
#endif

int32_t
tmap_vsw_process(tmap_vsw_t *vsw,
              const uint8_t *query, int32_t qlen,
              uint8_t *target, int32_t tlen, 
              tmap_vsw_result_t *result,
              int32_t *overflow, int32_t score_thr, int32_t is_rev)
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
  /*
  if(0 == tmap_vsw_update(vsw, query, qlen, target, tlen)) { // use the default
      tmap_bug();
      // TODO default algorithm
  }
  */
  if(tmap_vsw_wrapper_process(vsw->algorithm,
                           target, tlen, 
                           query, qlen, 
                           vsw->opt->score_match,
                           -vsw->opt->pen_mm,
                           -vsw->opt->pen_gapo,
                           -vsw->opt->pen_gape,
                           is_rev, 
                           vsw->query_start_clip, vsw->query_end_clip, 
                           &score, &target_end, &query_end, &n_best) <= INT16_MIN) {
      // try the default
      if(0 != vsw->type) {
          tmap_vsw_wrapper_process(vsw->algorithm_default,
                                   target, tlen, 
                                   query, qlen, 
                                   vsw->opt->score_match,
                                   -vsw->opt->pen_mm,
                                   -vsw->opt->pen_gapo,
                                   -vsw->opt->pen_gape,
                                   is_rev, 
                                   vsw->query_start_clip, vsw->query_end_clip, 
                                   &score, &target_end, &query_end, &n_best);
      }
  }
  if(score < score_thr) {
      query_end = target_end = -1;
      n_best = 0;
      score = INT32_MIN;
  }

#ifdef TMAP_VSW_DEBUG_CMP
  // TODO
  /*
  tmap_vsw_process_compare(vsw,
                         query, qlen,
                         target, tlen, 
                         result,
                         overflow, score_thr, is_rev,
                         TMAP_VSW_TYPE_S0, TMAP_VSW_TYPE_S3);
                         */
#endif
  
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
                                       is_rev, 
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
                                        NULL, NULL, is_rev);
          fprintf(stderr, "score=%d\n", score);
          */
          tmap_bug();
      }

      return result->score_fwd;
  }
}
