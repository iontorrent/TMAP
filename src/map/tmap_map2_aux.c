/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
#include "../util/tmap_alloc.h"
#include "../util/tmap_sort.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwtl.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_sa.h"
#include "../sw/tmap_sw.h"
#include "../sw/tmap_fsw.h"
#include "tmap_map_util.h"
#include "tmap_map2.h"
#include "tmap_map2_core.h"
#include "tmap_map2_chain.h"
#include "tmap_map2_aux.h"

#define __left_lt(a, b) ((a).end > (b).end)
TMAP_SORT_INIT(hit, tmap_map2_hit_t, __left_lt)

#define __hitG_lt(a, b) ((a).G > (b).G)
TMAP_SORT_INIT(hitG, tmap_map2_hit_t, __hitG_lt)

#define TMAP_MAP2_AUX_IS 0

#define __check_softclip(_softclip_type, _strand, _cigar, _n_cigar) do { \
    int32_t _type = _softclip_type; \
    if(1 == _strand) { \
        if(_type == TMAP_MAP_UTIL_SOFT_CLIP_LEFT) { \
            _type = TMAP_MAP_UTIL_SOFT_CLIP_RIGHT; \
        } \
        else if(_type == TMAP_MAP_UTIL_SOFT_CLIP_RIGHT) { \
            _type = TMAP_MAP_UTIL_SOFT_CLIP_LEFT; \
        } \
    } \
    switch(_type) { \
      case TMAP_MAP_UTIL_SOFT_CLIP_ALL: \
        break; \
      case TMAP_MAP_UTIL_SOFT_CLIP_LEFT: \
        if(BAM_CSOFT_CLIP == TMAP_SW_CIGAR_OP(_cigar[_n_cigar-1])) { \
            tmap_file_fprintf(tmap_file_stderr, "the strand is %d\n", strand); \
            tmap_error("found right soft clip", Warn, OutOfRange); \
        } \
        break; \
      case TMAP_MAP_UTIL_SOFT_CLIP_RIGHT: \
        if(BAM_CSOFT_CLIP == TMAP_SW_CIGAR_OP(_cigar[0])) { \
            tmap_file_fprintf(tmap_file_stderr, "the strand is %d\n", strand); \
            tmap_error("found left soft clip", Warn, OutOfRange); \
        } \
        break; \
      case TMAP_MAP_UTIL_SOFT_CLIP_NONE: \
        if(BAM_CSOFT_CLIP == TMAP_SW_CIGAR_OP(_cigar[_n_cigar-1])) { \
            tmap_file_fprintf(tmap_file_stderr, "the strand is %d\n", strand); \
            tmap_error("found right soft clip", Warn, OutOfRange); \
        } \
        if(BAM_CSOFT_CLIP == TMAP_SW_CIGAR_OP(_cigar[0])) { \
            tmap_file_fprintf(tmap_file_stderr, "the strand is %d\n", strand); \
            tmap_error("found left soft clip", Warn, OutOfRange); \
        } \
        break; \
      default: \
        tmap_error("soft clipping type was not recognized", Warn, OutOfRange); \
        break; \
    } \
} while(0)


int32_t
tmap_map2_aux_resolve_duphits(const tmap_bwt_t *bwt, const tmap_sa_t *sa, tmap_map2_aln_t *b, int32_t IS, int32_t min_as)
{
  int32_t i, j, n;
  if(b->n == 0) return 0;
  if(NULL != bwt) { // convert to chromosomal coordinates if suitable
      tmap_map2_aln_t *tmp_b;
      // copy over
      tmp_b = tmap_map2_aln_init();
      tmp_b->hits = b->hits;
      tmp_b->n = b->n;
      tmp_b->max = b->max;
      // nullify
      b->hits = NULL;
      b->n = b->max = 0;
      // go through tmp hits
      for(i = n = 0; i < tmp_b->n; ++i) {
          tmap_map2_hit_t *p = tmp_b->hits + i;
          if(0 != p->n_cigar) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }
          if(NULL != p->cigar) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }
          if(p->l - p->k + 1 <= IS) n += p->l - p->k + 1;
          else if(p->G > min_as) ++n;
      }
      // realloc
      tmap_map2_aln_realloc(b, n);
      b->n = n;
      // copy over
      for(i = j = 0; i < tmp_b->n; ++i) {
          tmap_map2_hit_t *p = tmp_b->hits + i;
          if(p->l - p->k + 1 <= IS) {
              uint32_t k;
              for(k = p->k; k <= p->l; ++k) {
                  b->hits[j] = *p;
                  b->hits[j].k = tmap_sa_pac_pos(sa, bwt, k);
                  b->hits[j].l = 0;
                  // do not copy over cigar
                  b->hits[j].cigar = NULL;
                  b->hits[j].n_cigar = 0;
                  ++j;
              }
          } else if(p->G > min_as) {
              b->hits[j] = *p;
              p->cigar = NULL; // nullify
              b->hits[j].k = tmap_sa_pac_pos(sa, bwt, p->k);
              b->hits[j].l = 0;
              b->hits[j].flag |= 1;
              // do not copy over cigar
              b->hits[j].cigar = NULL;
              b->hits[j].n_cigar = 0;
              ++j;
          }
      }
      tmap_map2_aln_destroy(tmp_b);
  }
  tmap_sort_introsort(hitG, b->n, b->hits);
  for(i = 1; i < b->n; ++i) {
      tmap_map2_hit_t *p = b->hits + i;
      if(p->G <= min_as) break;
      for(j = 0; j < i; ++j) {
          tmap_map2_hit_t *q = b->hits + j;
          int32_t compatible = 1;
          if(q->G == 0) continue;
          if(p->l == 0 && q->l == 0) {
              int32_t qol = (p->end < q->end? p->end : q->end) - (p->beg > q->beg? p->beg : q->beg);
              if(qol < 0) qol = 0;
              if((double)qol / (p->end - p->beg) > TMAP_MAP2_MASK_LEVEL 
                 || (double)qol / (q->end - q->beg) > TMAP_MAP2_MASK_LEVEL) {
                  int64_t tol = (int64_t)(p->k + p->len < q->k + q->len? p->k + p->len : q->k + q->len)
                    - (int64_t)(p->k > q->k? p->k : q->k);
                  if((double)tol / p->len > TMAP_MAP2_MASK_LEVEL || (double)tol / q->len > TMAP_MAP2_MASK_LEVEL)
                    compatible = 0;
              }
          }
          if(!compatible) {
              p->G = TMAP_MAP2_MINUS_INF; 
              break;
          }
      }
  }
  n = i;
  for(i = j = 0; i < n; ++i) {
      if(b->hits[i].G <= min_as) continue;
      if(i != j) b->hits[j++] = b->hits[i];
      else ++j;
  }
  b->n = j;
  return b->n;
}

/*
   static int32_t 
   tmap_map2_aux_resolve_query_overlaps(tmap_map2_aln_t *b, double mask_level, int32_t min_as)
   {
   int32_t i, j, n;
   if(b->n == 0) return 0;
   tmap_sort_introsort(hitG, b->n, b->hits);
   { // choose a random one
   int G0 = b->hits[0].G;
   for (i = 1; i < b->n; ++i)
   if (b->hits[i].G != G0) break;
   j = (int)(i * drand48());
   if (j) {
   tmap_map2_hit_t tmp;
   tmp = b->hits[0]; b->hits[0] = b->hits[j]; b->hits[j] = tmp;
   }
   }
   for(i = 1; i < b->n; ++i) {
   tmap_map2_hit_t *p = b->hits + i;
   int32_t all_compatible = 1;
   if(p->G == min_as) break;
   for(j = 0; j < i; ++j) {
   tmap_map2_hit_t *q = b->hits + j;
   int64_t tol = 0;
   int32_t qol, compatible = 0;
   double fol;
   if(q->G == min_as) continue;
   qol = (p->end < q->end? p->end : q->end) - (p->beg > q->beg? p->beg : q->beg);
   if(qol < 0) qol = 0;
   if(p->l == 0 && q->l == 0) {
   tol = (int64_t)(p->k + p->len < q->k + q->len? p->k + p->len : q->k + q->len)
   - (p->k > q->k? p->k : q->k);
   if(tol < 0) tol = 0;
   }
   fol = (double)qol / (p->end - p->beg < q->end - q->beg? p->end - p->beg : q->end - q->beg);
   if(fol < mask_level || (tol > 0 && qol < p->end - p->beg && qol < q->end - q->beg)) compatible = 1;
   if(!compatible) {
   if(q->G2 < p->G) q->G2 = p->G;
   all_compatible = 0;
   }
   }
   if(!all_compatible) p->G = min_as;
   }
   n = i;
   for(i = j = 0; i < n; ++i) {
   if(b->hits[i].G == min_as) continue;
   if(i != j) b->hits[j++] = b->hits[i];
   else ++j;
   }
   b->n = j;
   return j;
   }
   */

static void
tmap_map2_hit_nullify(tmap_map2_hit_t * hit)
{
  hit->k = hit->l = hit->flag = hit->n_seeds = 0;
  hit->len = hit->G = hit->G2 = hit->beg = hit->end = 0;
  hit->n_cigar = 0;
  hit->cigar = NULL;
}

tmap_map2_aln_t*
tmap_map2_aln_init()
{
  return tmap_calloc(1, sizeof(tmap_map2_aln_t), "a");
}

void 
tmap_map2_aln_realloc(tmap_map2_aln_t *a, int32_t n)
{
  int32_t i;
  if(NULL == a) return;
  // free old cigars
  if(n < a->n) {
      for(i=n;i<a->n;i++) {
          free(a->hits[i].cigar);
          tmap_map2_hit_nullify(&a->hits[i]);
      }
      a->n = n;
  }
  else if(a->max < n) { // allocate more memory
      i = a->max; // save for init
      a->max = (0 == a->max && n < 4) ? 4 : tmap_roundup32(n);
      // resize
      a->hits = tmap_realloc(a->hits, sizeof(tmap_map2_hit_t) * a->max, "a->hits");
      // init
      while(i < a->max) {
          tmap_map2_hit_nullify(&a->hits[i]);
          i++;
      }
  }
}

void 
tmap_map2_aln_destroy(tmap_map2_aln_t *a)
{
  int32_t i;
  if(NULL == a) return;
  for(i=0;i<a->n;i++) {
      free(a->hits[i].cigar);
  }
  free(a->hits);
  free(a);
}

// Note: this is the reverse, not the reverse compliment
#define tmap_map2_rseq_i(_refseq, _i) (tmap_refseq_seq_i(_refseq, _refseq->len-_i-1))

#define tmap_map2_aux_reverse_query(_query, _ql) \
  for(i=0;i<(_ql>>1);i++) { \
      uint8_t tmp = _query[i]; \
      _query[i] = _query[_ql-1-i]; \
      _query[_ql-1-i] = tmp; \
  }

static void 
tmap_map2_aux_extend_left(tmap_map_opt_t *opt, tmap_map2_aln_t *b, 
                          uint8_t *query,
                          int32_t query_length,
                          tmap_refseq_t *refseq,
                          uint8_t is_rev, uint8_t strand, uint8_t *_mem,
                          int32_t softclip_type)
{
  int32_t i, matrix[25];
  uint32_t k;
  uint8_t *target = NULL;
  int32_t target_length, to_fit;
  tmap_sw_param_t par;

  // since we reverse the query, we must reverse the soft-clipping
  // we do not need to reverse it twice
  //fprintf(stderr, "%s 1 softclip_type=%d is_rev=%d strand=%d\n", __func__, softclip_type, is_rev, strand);
  softclip_type = __tmap_map_util_reverse_soft_clipping(softclip_type);
  switch(softclip_type) {
    case TMAP_MAP_UTIL_SOFT_CLIP_ALL:
    case TMAP_MAP_UTIL_SOFT_CLIP_RIGHT:
      to_fit = 0;
      break;
    case TMAP_MAP_UTIL_SOFT_CLIP_LEFT:
    case TMAP_MAP_UTIL_SOFT_CLIP_NONE:
    default:
      to_fit = 1;
      break;
  }
  //fprintf(stderr, "%s 2 softclip_type=%d is_rev=%d strand=%d to_fit=%d\n", __func__, softclip_type, is_rev, strand, to_fit);

  par.matrix = matrix;
  __gen_ap(par, opt);
  // sort according to the descending order of query end
  tmap_sort_introsort(hit, b->n, b->hits);
  target_length = ((query_length + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + query_length;
  target = tmap_calloc(target_length, sizeof(uint8_t), "target");
  tmap_map2_aux_reverse_query(query, query_length); // reverse the query
  // core loop
  for(i = 0; i < b->n; ++i) {
      tmap_map2_hit_t *p = b->hits + i;
      int32_t lt = ((p->beg + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + query_length;
      int32_t score, j;
      tmap_sw_path_t path;
      /*
      fprintf(stderr, "%s before p->G=%d p->len=%d p->beg=%d p->end=%d p->k=%u p->l=%u score=%d strand=%d is_rev=%d\n",
              __func__, p->G, p->len, p->beg, p->end, p->k, p->l, score, strand, is_rev);
              */
      if(target_length < lt) tmap_error("target_length < lt", Exit, OutOfRange);
      p->n_seeds = 1;
      if(p->l || p->k == 0) {
          // we want to remove this mapping since it will cause improper
          // soft-clipping
          if(softclip_type != TMAP_MAP_UTIL_SOFT_CLIP_ALL) {
              p->G = TMAP_MAP2_MINUS_INF; 
          }
          continue;
      }
      if(0 == p->beg) continue; // no more base to extend
      for(j = score = 0; j < i; ++j) {
          tmap_map2_hit_t *q = b->hits + j;
          if(q->beg <= p->beg && q->k <= p->k && q->k + q->len >= p->k + p->len) {
              if(q->n_seeds < (1<<14) - 2) ++q->n_seeds;
              ++score;
          }
      }
      if(0 < score) { // contained in a previous alignment
          p->G = TMAP_MAP2_MINUS_INF; 
          continue;
      }
      if(p->k < lt) lt = p->k; // so we do not go off the beginning of refseq
      if(is_rev) {
          k = refseq->len - p->k + 1;
          lt = tmap_refseq_subseq(refseq, k, lt, target);
      }
      else {
          k = p->k;
          if(0 < lt) k = (k < lt) ? 1 : (k - lt + 1);
          lt = tmap_refseq_subseq(refseq, k, lt, target);
          // reverse it
          for(j=0;j<lt>>1;j++) {
              k = target[j]; 
              target[j] = target[lt-j-1];
              target[lt-j-1] = k;
          }
      }

      /*
      fprintf(stderr, "\nIn %s is_rev=%d strand=%d to_fit=%d p->beg=%d query_length=%d\n", 
              __func__, is_rev, strand, to_fit, p->beg, query_length);
      for(j=0;j<p->beg;j++) {
          fputc("ACGTN"[query[query_length-p->beg+j]], stderr);
      }
      fputc('\n', stderr);
      for(j=0;j<lt;j++) {
          fputc("ACGTN"[target[j]], stderr);
      }
      fputc('\n', stderr);
      */

      if(0 == to_fit) {
          score = tmap_sw_extend_core(target, lt, query + query_length - p->beg, p->beg, &par, &path, 0, p->G, _mem);
      }
      else {
          score = tmap_sw_extend_fitting_core(target, lt, query + query_length - p->beg, p->beg, &par, &path, 0, p->G, _mem);
      }

      if(1 == to_fit || (opt->score_thr < score && p->G <= score)) {
          p->G = score;
          p->len += path.i;
          p->beg -= path.j;
          p->k -= path.i;
      }
      /*
      fprintf(stderr, "%s after p->G=%d p->len=%d p->beg=%d p->end=%d p->k=%u\n",
              __func__, p->G, p->len, p->beg, p->end, p->k);
              */
  }
  tmap_map2_aux_reverse_query(query, query_length); // reverse back the query
  free(target);
}

static void 
tmap_map2_aux_extend_right(tmap_map_opt_t *opt, tmap_map2_aln_t *b, 
                           uint8_t *query,
                           int32_t query_length,
                           tmap_refseq_t *refseq,
                           uint8_t is_rev, uint8_t strand, uint8_t *_mem,
                           int32_t softclip_type)
{
  int32_t i, matrix[25];
  uint32_t k;
  uint8_t *target = NULL;
  int32_t target_length, to_fit;
  tmap_sw_param_t par;
  
  switch(softclip_type) {
    case TMAP_MAP_UTIL_SOFT_CLIP_ALL:
    case TMAP_MAP_UTIL_SOFT_CLIP_RIGHT:
      to_fit = 0;
      break;
    case TMAP_MAP_UTIL_SOFT_CLIP_LEFT:
    case TMAP_MAP_UTIL_SOFT_CLIP_NONE:
    default:
      to_fit = 1;
      break;
  }
  //fprintf(stderr, "%s softclip_type=%d is_rev=%d strand=%d to_fit=%d\n", __func__, softclip_type, is_rev, strand, to_fit);

  par.matrix = matrix;
  __gen_ap(par, opt);
  target_length = ((query_length + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + query_length;
  target = tmap_calloc(target_length, sizeof(uint8_t), "target");
  for(i = 0; i < b->n; ++i) {
      tmap_map2_hit_t *p = b->hits + i;
      int32_t lt = ((query_length - p->beg + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + query_length;
      int32_t j, score;
      tmap_sw_path_t path;
      /*
      fprintf(stderr, "%s before p->G=%d p->len=%d p->beg=%d p->end=%d p->k=%u score=%d strand=%d is_rev=%d\n",
              __func__, p->G, p->len, p->beg, p->end, p->k, score, strand, is_rev);
              */
      if(p->l) continue;
      if(query_length == p->end) continue; // no more base to extend
      if(is_rev) {
          k = refseq->len - p->k; // one-based
          if(0 < lt) { 
              if(k < lt) {
                  lt = k;
                  k = 1;
              }
              else {
                  k = k - lt + 1;
              }
          }
          lt = tmap_refseq_subseq(refseq, k, lt, target);
          // reverse it
          for(j=0;j<lt>>1;j++) {
              k = target[j]; 
              target[j] = target[lt-j-1];
              target[lt-j-1] = k;
          }
      }
      else {
          k = p->k + 1; // one-based
          lt = tmap_refseq_subseq(refseq, k, lt, target);
      }

      /*
      fprintf(stderr, "\nIn %s is_rev=%d strand=%d to_fit=%d p->beg=%d query_length=%d\n", 
              __func__, is_rev, strand, to_fit, p->beg, query_length);
      for(j=0;j<query_length - p->beg;j++) {
          fputc("ACGTN"[query[p->beg+j]], stderr);
      }
      fputc('\n', stderr);
      for(j=0;j<lt;j++) {
          fputc("ACGTN"[target[j]], stderr);
      }
      fputc('\n', stderr);
      */
      if(0 == to_fit) {
          score = tmap_sw_extend_core(target, lt, query + p->beg, query_length - p->beg, &par, &path, NULL, 1, _mem);
      }
      else {
          score = tmap_sw_extend_fitting_core(target, lt, query + p->beg, query_length - p->beg, &par, &path, NULL, 1, _mem);
      }

      if(1 == to_fit || (opt->score_thr < score && p->G <= score)) {
          p->G = score;
          p->len = path.i;
          p->end = path.j + p->beg;
      }
      /*
      fprintf(stderr, "%s after p->G=%d p->len=%d p->beg=%d p->end=%d p->k=%u\n",
              __func__, p->G, p->len, p->beg, p->end, p->k);
              */
  }
  free(target);
}

/* generate CIGAR array(s) in b->cigar[] */
static void 
tmap_map2_aux_gen_cigar(tmap_map_opt_t *opt, uint8_t *queries[2], 
                        int32_t query_length, tmap_refseq_t *refseq, tmap_map2_aln_t *b)
{
  uint8_t *target = NULL;
  int32_t i, n, matrix[25], target_len;
  tmap_sw_param_t par;
  tmap_sw_path_t *path;
  tmap_map_sam_t tmp_sam;

  par.matrix = matrix;
  __gen_ap(par, opt);
  target_len = ((query_length + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + query_length; // maximum possible target length
  target = tmap_calloc(target_len, sizeof(uint8_t), "target");
  path = tmap_calloc(target_len + query_length, sizeof(tmap_sw_path_t), "path");
  // memory clean up for b
  tmap_map2_aln_realloc(b, b->n);
  // generate CIGAR
  for(i = n = 0; i < b->n; ++i) {
      tmap_map2_hit_t *p = b->hits + i;
      tmap_map2_hit_t *q = b->hits + n;
      uint8_t *query;
      uint32_t seqid, coor;
      int32_t path_len, beg, end, added;
      int8_t strand;

      if(p->l) continue;

      strand = (p->flag & 0x10) ? 1 : 0;

      // adjust for contig boundaries
      if(tmap_refseq_pac2real(refseq, p->k, p->len, &seqid, &coor) <= 0) {
          if(1 == strand) { // reverse
              if(tmap_refseq_pac2real(refseq, p->k + p->len - 1, 1, &seqid, &coor) <= 0) {
                  // do nothing, this should fail later
              }
              else {
                  // move to the contig and position
                  p->k = refseq->annos[seqid].offset+1;
                  p->len = (target_len < refseq->annos[seqid].len) ? target_len : refseq->annos[seqid].len;
              }
          }
          else {
              if(tmap_refseq_pac2real(refseq, p->k, 1, &seqid, &coor) <= 0) {
                  // do nothing, this should fail later
              }
              else {
                  // move to the contig and position
                  p->k = refseq->annos[seqid].offset+1;
                  p->len = (target_len < refseq->annos[seqid].len) ? target_len : refseq->annos[seqid].len;
              }
          }
      }

      beg = (1 == strand) ? (query_length - p->end) : p->beg;
      end = (1 == strand) ? (query_length - p->beg) : p->end;

      // get more reference
      if(target_len < p->len) {
          target_len = p->len;
          target = tmap_realloc(target, target_len * sizeof(uint8_t), "target");
          path = tmap_realloc(path, (target_len + query_length) * sizeof(tmap_sw_path_t), "path");
      }

      query = queries[strand] + beg;
      if(p->len != tmap_refseq_subseq(refseq, p->k, p->len, target)) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }

      added = tmap_map_util_sw(&tmp_sam,
                               target, p->len,
                               query, end-beg,
                               0, 0,
                               &par, path, &path_len,
                               opt->score_thr, opt->softclip_type, strand);

      if(1 == added) {
          if(0 == path_len) {
              tmap_error("0 == path_len", Exit, OutOfRange);
          }
          if(0 == tmp_sam.n_cigar) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }


          // copy over
          if(n < i) {
              (*q) = (*p);
          }

          // adjust the alignment length
          q->len = path[0].i - path[path_len-1].i + 1;
          q->k += tmp_sam.pos; // adjust the alignment start
          q->cigar = tmp_sam.cigar;
          q->n_cigar = tmp_sam.n_cigar;
          q->G = tmp_sam.score;
          
          __check_softclip(opt->softclip_type, strand, q->cigar, q->n_cigar);

          /*
          fprintf(stderr, "strand=%d beg=%d end=%d p->beg=%d p->end=%d query_length=%d\n", 
                  strand, beg, end, p->beg, p->end, query_length);
                  */

          // add latent soft clipping at the front
          if(0 < beg){
              if(BAM_CSOFT_CLIP == TMAP_SW_CIGAR_OP(q->cigar[0])) {
                  TMAP_SW_CIGAR_ADD_LENGTH(q->cigar[0], beg);
              }
              else {
                  q->cigar = tmap_realloc(q->cigar, sizeof(uint32_t) * (q->n_cigar + 1), "b->cigar");
                  memmove(q->cigar + 1, q->cigar, q->n_cigar * sizeof(uint32_t));
                  TMAP_SW_CIGAR_STORE(q->cigar[0], BAM_CSOFT_CLIP, beg);
                  ++q->n_cigar;
              }
          }
          // soft clipping at the end
          if(end < query_length) { 
              if(BAM_CSOFT_CLIP == TMAP_SW_CIGAR_OP(q->cigar[q->n_cigar-1])) {
                  TMAP_SW_CIGAR_ADD_LENGTH(q->cigar[q->n_cigar-1], (query_length - end));
              }
              else {
                  q->cigar = tmap_realloc(q->cigar, sizeof(uint32_t) * (q->n_cigar + 1), "b->cigar");
                  TMAP_SW_CIGAR_STORE(q->cigar[q->n_cigar], BAM_CSOFT_CLIP, (query_length-end));
                  ++q->n_cigar;
              }
          }

          // check soft clipping
          __check_softclip(opt->softclip_type, strand, q->cigar, q->n_cigar);

          n++;
      }
  }
  free(target); free(path);

  // reallocate
  tmap_map2_aln_realloc(b, n);
}

static void 
tmap_map2_aux_merge_hits(tmap_map2_aln_t *b[2], int32_t l, int32_t is_reverse, int32_t softclip_type)
{
  int32_t i;
  tmap_map2_aln_realloc(b[0], b[0]->n + b[1]->n);
  for(i = 0; i < b[1]->n; ++i) {
      tmap_map2_hit_t *p = b[0]->hits + b[0]->n + i;
      *p = b[1]->hits[i];
      if(is_reverse) {
          int32_t x = p->beg;
          p->beg = l - p->end;
          p->end = l - x;
          p->flag |= 0x10;
      }
  }
  b[0]->n += b[1]->n;
  tmap_map2_aln_destroy(b[1]);
  b[1] = NULL;

  if(TMAP_MAP_UTIL_SOFT_CLIP_NONE == softclip_type) { // flag non-global hits
      for(i=0;i<b[0]->n;i++) {
          tmap_map2_hit_t *p = b[0]->hits + i;
          if(0 != p->beg || l != p->end) {
              p->G = TMAP_MAP2_MINUS_INF;
          }
      }
  }
}

static tmap_map2_aln_t *
tmap_map2_aux_aln(tmap_map_opt_t *opt, tmap_refseq_t *refseq, 
                  tmap_bwt_t *target_bwt, tmap_sa_t *target_sa,
                  tmap_string_t *seq[2], int32_t is_rev, tmap_map2_global_mempool_t *pool)
{
  tmap_map2_aln_t *b[2], **bb[2];
  int32_t k, softclip_type;

  for(k = 0; k < 2; ++k) {
      tmap_bwtl_t *query = tmap_bwtl_seq2bwtl(seq[k]->l, (uint8_t*)seq[k]->s);
      bb[k] = tmap_map2_core_aln(opt, query, target_bwt, target_sa, pool);
      tmap_bwtl_destroy(query);
  }
  b[0] = bb[0][1]; b[1] = bb[1][1]; // bb[*][1] are "narrow SA hits"
  tmap_map2_chain_filter(opt, seq[0]->l, b);
  // TODO: could we skip all of this extension and just generate the alignments
  // directly?
  for(k = 0; k < 2; ++k) {
      softclip_type = opt->softclip_type;
      if(k ^ is_rev) { // one or the other, but not both
          softclip_type = __tmap_map_util_reverse_soft_clipping(softclip_type);
      }
      // This causes problems with soft-clipping
      /*
      tmap_map2_aux_extend_left(opt, bb[k][1], (uint8_t*)seq[k]->s, seq[k]->l, refseq, is_rev, k, pool->aln_mem, softclip_type);
      tmap_map2_aux_merge_hits(bb[k], seq[k]->l, 0, TMAP_MAP_UTIL_SOFT_CLIP_ALL); // bb[k][1] and bb[k][0] are merged into bb[k][0]
      tmap_map2_aux_resolve_duphits(NULL, NULL, bb[k][0], TMAP_MAP2_AUX_IS, (TMAP_MAP_UTIL_SOFT_CLIP_NONE == softclip_type) ? TMAP_MAP2_MINUS_INF : 0);
      tmap_map2_aux_extend_right(opt, bb[k][0], (uint8_t*)seq[k]->s, seq[k]->l, refseq, is_rev, k, pool->aln_mem, softclip_type);
      */
      // Note, we are examining all hits, not just extending left on narrow, and
      // right on all...
      tmap_map2_aux_extend_left(opt, bb[k][0], (uint8_t*)seq[k]->s, seq[k]->l, refseq, is_rev, k, pool->aln_mem, softclip_type);
      tmap_map2_aux_resolve_duphits(NULL, NULL, bb[k][0], TMAP_MAP2_AUX_IS, (TMAP_MAP_UTIL_SOFT_CLIP_NONE == softclip_type) ? TMAP_MAP2_MINUS_INF : 0);
      tmap_map2_aux_extend_right(opt, bb[k][0], (uint8_t*)seq[k]->s, seq[k]->l, refseq, is_rev, k, pool->aln_mem, softclip_type);
      tmap_map2_aln_destroy(bb[k][1]); // ignore not so repetitive hits, since we want all hits
      b[k] = bb[k][0];
      free(bb[k]);		
  }
  tmap_map2_aux_merge_hits(b, seq[0]->l, 1, softclip_type); // b[1] and b[0] are merged into b[0]
  // Note: this will give duplicate mappings
  //tmap_map2_aux_resolve_query_overlaps(b[0], opt->mask_level, (TMAP_MAP_UTIL_SOFT_CLIP_NONE == opt->softclip_type) ? TMAP_MAP2_MINUS_INF : 0);

  return b[0];
}

/* set ->flag to records the origin of the hit (to forward bwt or reverse bwt) */
static void 
tmap_map2_aux_flag_fr(tmap_map2_aln_t *b[2])
{
  int32_t i, j;
  if(NULL != b[0]) {
      for(i = 0; i < b[0]->n; ++i) {
          tmap_map2_hit_t *p = b[0]->hits + i;
          p->flag |= 0x10000;
      }
  }
  if(NULL != b[1]) {
      for(i = 0; i < b[1]->n; ++i) {
          tmap_map2_hit_t *p = b[1]->hits + i;
          p->flag |= 0x20000;
      }
  }
  if(NULL != b[0] && NULL != b[1] && 0 < b[0]->n && 0 < b[1]->n) {
      for(i = 0; i < b[0]->n; ++i) {
          tmap_map2_hit_t *p = b[0]->hits + i;
          for(j = 0; j < b[1]->n; ++j) {
              tmap_map2_hit_t *q = b[1]->hits + j;
              if(q->beg == p->beg && q->end == p->end && q->k == p->k && q->len == p->len && q->G == p->G) {
                  q->flag |= 0x30000; p->flag |= 0x30000;
                  break;
              }
          }
      }
  }
}

/*
static int32_t 
tmap_map2_aux_fix_cigar(tmap_refseq_t *refseq, tmap_map2_hit_t *p, int32_t n_cigar, uint32_t *cigar)
{
  int32_t lq;
  int32_t x, y, i;
  uint32_t coor, seqid, refl;

  if(0 == tmap_refseq_pac2real(refseq, p->k, p->len, &seqid, &coor)) {
      return -1;
  }
  if(n_cigar == 0) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }

  refl = refseq->annos[seqid].len;
  x = coor-1; y = 0;
  // test if the alignment goes beyond the boundary
  for(i = 0; i < n_cigar; ++i) {
      int32_t op = TMAP_SW_CIGAR_OP(cigar[i]);
      int32_t ln = TMAP_SW_CIGAR_LENGTH(cigar[i]);
      if(op == 1 || op == 4 || op == 5) y += ln;
      else if(op == 2) x += ln;
      else x += ln, y += ln;
  }
  lq = y; // length of the query sequence
  if(x > refl) { // then fix it
      int32_t j, nc, mq[2], nlen[2];
      uint32_t *cn, kk = 0;
      nc = mq[0] = mq[1] = nlen[0] = nlen[1] = 0;
      cn = tmap_calloc(n_cigar + 3, sizeof(uint32_t), "cn");
      x = coor-1; y = 0;
      for(i = j = 0; i < n_cigar; ++i) {
          int32_t op = TMAP_SW_CIGAR_OP(cigar[i]);
          int32_t ln = TMAP_SW_CIGAR_LENGTH(cigar[i]);
          if(BAM_CINS == op || BAM_CSOFT_CLIP == op || BAM_CHARD_CLIP == op) { // ins or clipping
              y += ln;
              cn[j++] = cigar[i];
          } else if(BAM_CDEL == op) { // del
              if(x + ln >= refl && nc == 0) {
                  TMAP_SW_CIGAR_STORE(cn[j++], BAM_CSOFT_CLIP, (uint32_t)(lq-y));
                  nc = j;
                  TMAP_SW_CIGAR_STORE(cn[j++], BAM_CSOFT_CLIP, (uint32_t)y);
                  kk = p->k + (x + ln - refl);
                  nlen[0] = x - (coor-1);
                  nlen[1] = p->len - nlen[0] - ln;
              } else cn[j++] = cigar[i];
              x += ln;
          } else if(op == BAM_CMATCH) { // match
              if(x + ln >= refl && nc == 0) {
                  // FIXME: not consider a special case where a split right between M and I
                  TMAP_SW_CIGAR_STORE(cn[j++], BAM_CMATCH, (uint32_t)(refl-x)); // write M
                  TMAP_SW_CIGAR_STORE(cn[j++], BAM_CSOFT_CLIP, (uint32_t)(lq-y-(refl-x))); // write S
                  nc = j;
                  mq[0] += refl - x;
                  TMAP_SW_CIGAR_STORE(cn[j++], BAM_CSOFT_CLIP, (uint32_t)(y+(refl-x))); // write S
                  if(x + ln - refl) TMAP_SW_CIGAR_STORE(cn[j++], BAM_CMATCH, (uint32_t)(x+ln-refl)); // write M
                  mq[1] += x + ln - refl;
                  kk = refseq->annos[seqid].offset + refl;
                  nlen[0] = refl - (coor-1);
                  nlen[1] = p->len - nlen[0];
              } else {
                  cn[j++] = cigar[i];
                  mq[nc?1:0] += ln;
              }
              x += ln; y += ln;
          }
      }
      if(mq[0] > mq[1]) { // then take the first alignment
          n_cigar = nc;
          memcpy(cigar, cn, 4 * nc);
          p->len = nlen[0];
      } else {
          p->k = kk; p->len = nlen[1];
          n_cigar = j - nc;
          memcpy(cigar, cn + nc, 4 * (j - nc));
      }
      free(cn);
  }
  return n_cigar;
}
*/

static tmap_map_sams_t *
tmap_map1_aux_store_hits(tmap_refseq_t *refseq, tmap_map_opt_t *opt, 
                         tmap_map2_aln_t *aln)
{
  int32_t i, j;
  tmap_map_sams_t *sams = NULL;

  if(NULL == aln) return NULL;

  sams = tmap_map_sams_init();
  tmap_map_sams_realloc(sams, aln->n);

  for(i=j=0;i<aln->n;i++) {
      tmap_map2_hit_t *p = aln->hits + i;
      uint32_t seqid = 0, coor = 0;
      int32_t qual;
      tmap_map_sam_t *sam = &sams->sams[j];

      if(p->l == 0) {
          /*
          aln->hits[i].n_cigar = tmap_map2_aux_fix_cigar(refseq, p, aln->hits[i].n_cigar, aln->hits[i].cigar);
          */
          if(aln->hits[i].n_cigar <= 0) {
              continue; // no cigar
          }
          if(tmap_refseq_pac2real(refseq, p->k, p->len, &seqid, &coor) <= 0) {
              continue; // spanning two or more chromosomes
          }
      }

      sam->strand = (p->flag & 0x10) ? 1 : 0; // strand
      sam->seqid = seqid;
      sam->pos = coor-1; // make it zero-based
      if(p->l == 0) {
          // estimate mapping quality
          double c = 1.0;	
          int32_t subo = (p->G2 > opt->score_thr) ? p->G2 : opt->score_thr;
          if(p->flag>>16 == 1 || p->flag>>16 == 2) c *= .5;
          if(p->n_seeds < 2) c *= .2;
          qual = (int)(c * (p->G - subo) * (250.0 / p->G + 0.03 / opt->score_match) + .499);
          if(qual > 250) qual = 250;
          if((p->flag & 0x1)) {
              qual = 0;
              p->G2 = p->G; // Note: the flag indicates a repetitive match, so we need to update the sub-optimal score
          }
          sam->mapq = qual;

          // copy cigar memory
          sam->n_cigar = aln->hits[i].n_cigar;
          sam->cigar = aln->hits[i].cigar;
          aln->hits[i].n_cigar = 0;
          aln->hits[i].cigar = NULL;
      } 
      else {
          sam->mapq = 0;
          sam->n_cigar = 0;
          sam->cigar = NULL;
      }
      sam->algo_id = TMAP_MAP_ALGO_MAP2;
      sam->algo_stage = 0;
      sam->score = p->G;
      sam->score_subo = p->G2;
      // auxiliary data
      tmap_map_sam_malloc_aux(sam, TMAP_MAP_ALGO_MAP2);
      sam->aux.map2_aux->XE = p->n_seeds;
      sam->aux.map2_aux->XF = p->flag >> 16;
      if(p->l) {
          sam->aux.map2_aux->XI = p->l - p->k + 1;
      }
      else {
          sam->aux.map2_aux->XI = 0;
      }
      j++;
  }
  if(j != aln->n) {
      tmap_map_sams_realloc(sams, j);
  }

  return sams;
}

tmap_map_sams_t *
tmap_map2_aux_core(tmap_map_opt_t *_opt,
                   tmap_seq_t *seqs[4],
                   tmap_refseq_t *refseq,
                   tmap_bwt_t *bwt[2],
                   tmap_sa_t *sa[2],
                   tmap_map2_global_mempool_t *pool)
{
  tmap_map_opt_t opt;
  tmap_seq_t *orig_seq = NULL;
  tmap_string_t *seq[2]={NULL, NULL};
  tmap_string_t *rseq[2]={NULL, NULL};
  tmap_map_sams_t *sams = NULL;
  tmap_map2_aln_t *b[2]={NULL,NULL};
  tmap_string_t *bases = NULL;
  uint8_t *_seq[2];
  int32_t i, k, l, num_n;

  opt = (*_opt);

  // sequence length
  bases = tmap_seq_get_bases(seqs[0]);
  l = bases->l;

  // set opt->score_thr
  opt.score_thr = _opt->score_thr; // reset opt->score_thr
  if(opt.score_thr < log(l) * opt.length_coef) opt.score_thr = (int)(log(l) * opt.length_coef + .499);
  if(pool->max_l < l) { // then enlarge working space for tmap_sw_extend_core()
      int32_t tmp = ((l + 1) / 2 * opt.score_match + opt.pen_gape) / opt.pen_gape + l;
      pool->max_l = l;
      pool->aln_mem = tmap_realloc(pool->aln_mem, sizeof(uint8_t) * (tmp + 2) * 24, "pool->aln_mem");
  }

  // set opt->bw
  opt.bw = _opt->bw;
  k = (l * opt.score_match - 2 * opt.pen_gapo) / (2 * opt.pen_gape + opt.score_match);
  i = (l * opt.score_match - opt.score_match - opt.score_thr) / opt.pen_gape;
  if(k > i) k = i;
  if(k < 1) k = 1; // I do not know if k==0 causes troubles
  opt.bw = _opt->bw < k ? _opt->bw: k;

  // get the number of Ns
  for(i=num_n=0;i<l;i++) {
      uint8_t c = (uint8_t)tmap_nt_char_to_int[(int)bases->s[i]];
      if(c >= 4) num_n++; // FIXME: ambiguous bases are not properly handled
  }

  // will we always be lower than the score threshold
  if((l*opt.score_match) + (num_n*opt.pen_mm) < opt.score_thr) {
      return tmap_map_sams_init(0);
  }
  
  // save sequences
  seq[0] = tmap_seq_get_bases(seqs[0]); 
  seq[1] = tmap_seq_get_bases(seqs[1]);
  rseq[0] = tmap_seq_get_bases(seqs[2]); 
  rseq[1] = tmap_seq_get_bases(seqs[3]);

  // handle ambiguous bases
  if(0 < num_n) {
      // save original to de-randomize later
      orig_seq = tmap_seq_clone(seqs[0]);

      // randomize
      for(i=0;i<l;i++) {
          uint8_t c = (uint8_t)bases->s[i];
          if(c >= 4) {
              c = (int)(drand48() * 4); // FIXME: ambiguous bases are not properly handled
              seq[0]->s[i] = c; // original
              seq[1]->s[l-1-i] = 3 - c; // reverse compliment
              rseq[0]->s[l-1-i] = c; // reverse 
              rseq[1]->s[i] = 3 - c; // compliment
          }
      }
  }

  // alignment
  b[0] = tmap_map2_aux_aln(&opt, refseq, bwt[0], sa[0], seq, 0, pool);
  for(k = 0; k < b[0]->n; ++k) {
      if(b[0]->hits[k].n_seeds < opt.seeds_rev) break;
  } 
  if(k < b[0]->n) {
      b[1] = tmap_map2_aux_aln(&opt, refseq, bwt[1], sa[1], rseq, 1, pool);
      for(i = 0; i < b[1]->n; ++i) {
          tmap_map2_hit_t *p = b[1]->hits + i;
          int x = p->beg;
          p->beg = l - p->end;
          p->end = l - x;
          if(p->l == 0) {
              p->k = refseq->len - (p->k + p->len);
          }
      }
      tmap_map2_aux_merge_hits(b, l, 0, opt.softclip_type);
      tmap_map2_aux_resolve_duphits(NULL, NULL, b[0], TMAP_MAP2_AUX_IS, (0 == opt.softclip_type) ? 0 : TMAP_MAP2_MINUS_INF);
      // Note: this will give duplicate mappings
      //tmap_map2_aux_resolve_query_overlaps(b[0], opt.mask_level, (TMAP_MAP_UTIL_SOFT_CLIP_NONE == opt->softclip_type) ? TMAP_MAP2_MINUS_INF : 0);
  } else b[1] = 0;
      
  // set the flag to forward/reverse
  tmap_map2_aux_flag_fr(b);
  
  // make one-based for pac2real
  for(i = 0; i < b[0]->n; ++i) {
      b[0]->hits[i].k++;
  }

  // generate CIGAR and print SAM
  _seq[0] = (uint8_t*)seq[0]->s;
  _seq[1] = (uint8_t*)seq[1]->s;
  tmap_map2_aux_gen_cigar(&opt, _seq, l, refseq, b[0]);
  sams = tmap_map1_aux_store_hits(refseq, &opt, b[0]);

  // remove duplicate alignments
  tmap_map_util_remove_duplicates(sams, opt.dup_window);

  // free
  tmap_map2_aln_destroy(b[0]);

  // revert ambiguous bases 
  if(0 < num_n) {
      // de-randomize
      bases = tmap_seq_get_bases(orig_seq);
      for(i=0;i<l;i++) {
          uint8_t c = (uint8_t)bases->s[i];
          if(c >= 4) { 
              seq[0]->s[i] = c; // original
              seq[1]->s[l-1-i] = 3 - c; // reverse compliment
              rseq[0]->s[l-1-i] = c; // reverse 
              rseq[1]->s[i] = 3 - c; // compliment
          }
      }
      tmap_seq_destroy(orig_seq);
  }

  return sams;
}
