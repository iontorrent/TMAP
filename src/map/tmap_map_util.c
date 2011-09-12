/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "../util/tmap_alloc.h"
#include "../util/tmap_error.h"
#include "../util/tmap_sam_print.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sort.h"
#include "../util/tmap_definitions.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_sa.h"
#include "../sw/tmap_sw.h"
#include "../sw/tmap_fsw.h"
#include "../sw/tmap_vsw.h"
#include "tmap_map_opt.h"
#include "tmap_map_util.h"

#define tmap_map_util_reverse_query(_query, _ql, _i) \
  for(_i=0;_i<(_ql>>1);_i++) { \
      uint8_t _tmp = _query[_i]; \
      _query[_i] = _query[_ql-_i-1]; \
      _query[_ql-_i-1] = _tmp; \
  }

// sort by strand, min-seqid, min-position
#define __tmap_map_sam_sort_coord_lt(a, b) (  ((a).strand < (b).strand) \
                                            || ( (a).strand == (b).strand && (a).seqid < (b).seqid) \
                                            || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos < (b).pos ) \
                                            ? 1 : 0 )

#define __tmap_map_sam_sort_coord_end_lt(a, b) (  ((a).strand < (b).strand) \
                                            || ( (a).strand == (b).strand && (a).seqid < (b).seqid) \
                                            || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos + (a).target_len < (b).pos + (b).target_len) \
                                            ? 1 : 0 )

// sort by strand, min-seqid, min-position, max score
#define __tmap_map_sam_sort_coord_score_lt(a, b) (  ((a).strand < (b).strand) \
                                            || ( (a).strand == (b).strand && (a).seqid < (b).seqid) \
                                            || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos < (b).pos ) \
                                            || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos == (b).pos && (a).score > (b).score) \
                                            ? 1 : 0 )

TMAP_SORT_INIT(tmap_map_sam_sort_coord, tmap_map_sam_t, __tmap_map_sam_sort_coord_lt)
TMAP_SORT_INIT(tmap_map_sam_sort_coord_end, tmap_map_sam_t, __tmap_map_sam_sort_coord_end_lt)
TMAP_SORT_INIT(tmap_map_sam_sort_coord_score, tmap_map_sam_t, __tmap_map_sam_sort_coord_score_lt)

void
tmap_map_sam_malloc_aux(tmap_map_sam_t *s, int32_t algo_id)
{
  switch(s->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      s->aux.map1_aux = tmap_calloc(1, sizeof(tmap_map_map1_aux_t), "s->aux.map1_aux");
      break;
    case TMAP_MAP_ALGO_MAP2:
      s->aux.map2_aux = tmap_calloc(1, sizeof(tmap_map_map2_aux_t), "s->aux.map2_aux");
      break;
    case TMAP_MAP_ALGO_MAP3:
      s->aux.map3_aux = tmap_calloc(1, sizeof(tmap_map_map3_aux_t), "s->aux.map3_aux");
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      s->aux.map_vsw_aux = tmap_calloc(1, sizeof(tmap_map_map_vsw_aux_t), "s->aux.map_vsw_aux");
      break;
    default:
      break;
  }
}

inline void
tmap_map_sam_destroy_aux(tmap_map_sam_t *s)
{
  switch(s->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      free(s->aux.map1_aux);
      s->aux.map1_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP2:
      free(s->aux.map2_aux);
      s->aux.map2_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP3:
      free(s->aux.map3_aux);
      s->aux.map3_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      free(s->aux.map_vsw_aux);
      s->aux.map_vsw_aux = NULL;
      break;
    default:
      break;
  }
}

void
tmap_map_sam_destroy(tmap_map_sam_t *s)
{
  tmap_map_sam_destroy_aux(s);
  free(s->cigar);
  s->cigar = NULL;
  s->n_cigar = 0;
}

tmap_map_sams_t *
tmap_map_sams_init(tmap_map_sams_t *prev)
{
  tmap_map_sams_t *sams = tmap_calloc(1, sizeof(tmap_map_sams_t), "sams");
  sams->sams = NULL;
  sams->n = 0;
  if(NULL != prev) sams->max = prev->max;
  return sams;
}

void
tmap_map_sams_realloc(tmap_map_sams_t *s, int32_t n)
{
  int32_t i;
  if(n == s->n) return; 
  for(i=n;i<s->n;i++) {
      tmap_map_sam_destroy(&s->sams[i]);
  }
  s->sams = tmap_realloc(s->sams, sizeof(tmap_map_sam_t) * n, "s->sams");
  for(i=s->n;i<n;i++) {
      // nullify
      s->sams[i].algo_id = TMAP_MAP_ALGO_NONE;
      s->sams[i].n_cigar = 0;
      s->sams[i].cigar = NULL;
      s->sams[i].aux.map1_aux = NULL;
      s->sams[i].aux.map2_aux = NULL;
      s->sams[i].aux.map3_aux = NULL;
      s->sams[i].aux.map_vsw_aux = NULL;
      s->sams[i].ascore = INT32_MIN;
  }
  s->n = n;
}

void
tmap_map_sams_destroy(tmap_map_sams_t *s)
{
  int32_t i;
  if(NULL == s) return;
  for(i=0;i<s->n;i++) {
      tmap_map_sam_destroy(&s->sams[i]);
  }
  free(s->sams);
  free(s);
}

static inline void
tmap_map_sam_copy(tmap_map_sam_t *dest, tmap_map_sam_t *src)
{
  int32_t i;
  // shallow copy
  (*dest) = (*src);
  // aux data
  tmap_map_sam_malloc_aux(dest, src->algo_id);
  switch(src->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      (*dest->aux.map1_aux) = (*src->aux.map1_aux);
      break;
    case TMAP_MAP_ALGO_MAP2:
      (*dest->aux.map2_aux) = (*src->aux.map2_aux);
      break;
    case TMAP_MAP_ALGO_MAP3:
      (*dest->aux.map3_aux) = (*src->aux.map3_aux);
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      (*dest->aux.map_vsw_aux) = (*src->aux.map_vsw_aux);
      break;
    default:
      break;
  }
  // cigar
  if(0 < src->n_cigar && NULL != src->cigar) {
      dest->cigar = tmap_malloc(sizeof(uint32_t) * (1 + dest->n_cigar), "dest->cigar");
      for(i=0;i<dest->n_cigar;i++) {
          dest->cigar[i] = src->cigar[i];
      }
  }
}

void
tmap_map_sams_merge(tmap_map_sams_t *dest, tmap_map_sams_t *src) 
{
  int32_t i, j;
  if(NULL == src || 0 == src->n) return;

  j = dest->n;
  tmap_map_sams_realloc(dest, dest->n + src->n);
  for(i=0;i<src->n;i++,j++) {
      tmap_map_sam_copy(&dest->sams[j], &src->sams[i]);
  }
}

tmap_map_sams_t *
tmap_map_sams_clone(tmap_map_sams_t *src)
{
  int32_t i;
  tmap_map_sams_t *dest = NULL;
  
  // init
  dest = tmap_map_sams_init(src);
  if(0 == src->n) return dest;

  // realloc
  tmap_map_sams_realloc(dest, src->n);
  // copy over data
  for(i=0;i<src->n;i++) {
      tmap_map_sam_copy(&dest->sams[i], &src->sams[i]);
  }

  return dest;
}

void
tmap_map_sam_copy_and_nullify(tmap_map_sam_t *dest, tmap_map_sam_t *src)
{
  (*dest) = (*src);
  src->n_cigar = 0;
  src->cigar = NULL;
  switch(src->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      src->aux.map1_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP2:
      src->aux.map2_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP3:
      src->aux.map3_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      src->aux.map_vsw_aux = NULL;
      break;
    default:
      break;
  }
}

static void
tmap_map_sam_print(tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_map_sam_t *sam, int32_t sam_sff_tags, 
                   int32_t nh, int32_t aln_num, int32_t end_num, int32_t mate_unmapped, tmap_map_sam_t *mate)
{
  uint32_t mate_strand, mate_seqid, mate_pos, mate_tlen;
  // mate info
  mate_strand = mate_seqid = mate_pos = mate_tlen = 0;
  if(NULL != mate && 0 == mate_unmapped) {
      mate_strand = mate->strand;
      mate_seqid = mate->seqid;
      mate_pos = mate->pos;
  }
  if(NULL == sam) { // unmapped
      tmap_sam_print_unmapped(tmap_file_stdout, seq, sam_sff_tags, refseq,
                              end_num, mate_unmapped, 0,
                              mate_strand, mate_seqid, mate_pos);
  }
  else {
      // Note: samtools does not like this value
      if(INT32_MIN == sam->score_subo) {
          sam->score_subo++;
      }
      // Set the mate tlen
      if(NULL != mate && 0 == mate_unmapped) {
          // NB: assumes 5'->3' ordering of the fragments
          if(mate->pos < sam->pos) {
              mate_tlen = sam->pos + sam->target_len - mate->pos;
          }
          else {
              mate_tlen = mate->pos + mate->target_len - sam->pos;
          }
          // NB: first fragment is always positive, the rest are always negative
          if(1 == end_num) { // first end
              if(mate_tlen < 0) mate_tlen = -mate_tlen;
          }
          else { // second end
              if(0 < mate_tlen) mate_tlen = -mate_tlen;
          }
      }
      switch(sam->algo_id) {
        case TMAP_MAP_ALGO_MAP1:
          tmap_sam_print_mapped(tmap_file_stdout, seq, sam_sff_tags, refseq, 
                                sam->strand, sam->seqid, sam->pos, aln_num,
                                end_num, mate_unmapped, 0,
                                mate_strand, mate_seqid, mate_pos, mate_tlen,
                                sam->mapq, sam->cigar, sam->n_cigar,
                                sam->score, sam->ascore, nh, sam->algo_id, sam->algo_stage, "");
          break;
        case TMAP_MAP_ALGO_MAP2:
          if(0 < sam->aux.map2_aux->XI) {
              tmap_sam_print_mapped(tmap_file_stdout, seq, sam_sff_tags, refseq, 
                                    sam->strand, sam->seqid, sam->pos, aln_num,
                                    end_num, mate_unmapped, 0,
                                    mate_strand, mate_seqid, mate_pos, mate_tlen,
                                    sam->mapq, sam->cigar, sam->n_cigar,
                                    sam->score, sam->ascore, nh, sam->algo_id, sam->algo_stage, 
                                    "\tXS:i:%d\tXT:i:%d\t\tXF:i:%d\tXE:i:%d\tXI:i:%d",
                                    sam->score_subo,
                                    sam->n_seeds,
                                    sam->aux.map2_aux->XF, sam->aux.map2_aux->XE, 
                                    sam->aux.map2_aux->XI);
          }
          else {
              tmap_sam_print_mapped(tmap_file_stdout, seq, sam_sff_tags, refseq, 
                                    sam->strand, sam->seqid, sam->pos, aln_num,
                                    end_num, mate_unmapped, 0,
                                    mate_strand, mate_seqid, mate_pos, mate_tlen,
                                    sam->mapq, sam->cigar, sam->n_cigar,
                                    sam->score, sam->ascore, nh, sam->algo_id, sam->algo_stage, 
                                    "\tXS:i:%d\tXT:i:%d\tXF:i:%d\tXE:i:%d",
                                    sam->score_subo,
                                    sam->n_seeds,
                                    sam->aux.map2_aux->XF, sam->aux.map2_aux->XE);
          }
          break;
        case TMAP_MAP_ALGO_MAP3:
          tmap_sam_print_mapped(tmap_file_stdout, seq, sam_sff_tags, refseq, 
                                sam->strand, sam->seqid, sam->pos, aln_num,
                                end_num, mate_unmapped, 0,
                                mate_strand, mate_seqid, mate_pos, mate_tlen,
                                sam->mapq, sam->cigar, sam->n_cigar,
                                sam->score, sam->ascore, nh, sam->algo_id, sam->algo_stage, 
                                "\tXS:i:%d\tXT:i:%d",
                                sam->score_subo,
                                sam->n_seeds);
          break;
        case TMAP_MAP_ALGO_MAPVSW:
          tmap_sam_print_mapped(tmap_file_stdout, seq, sam_sff_tags, refseq, 
                                sam->strand, sam->seqid, sam->pos, aln_num,
                                end_num, mate_unmapped, 0,
                                mate_strand, mate_seqid, mate_pos, mate_tlen,
                                sam->mapq, sam->cigar, sam->n_cigar,
                                sam->score, sam->ascore, nh, sam->algo_id, sam->algo_stage, 
                                "\tXS:i:%d",
                                sam->score_subo);
          break;
      }
  }
}

void 
tmap_map_sams_print(tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_map_sams_t *sams, int32_t end_num,
                    tmap_map_sams_t *mates, int32_t sam_sff_tags) 
{
  int32_t i;
  tmap_map_sam_t *mate = NULL;
  int32_t mate_unmapped = 0;
  if(NULL != mates) {
      if(0 < mates->n) {
          // assumes the mates are sorted by their alignment score
          mate = &mates->sams[0];
      }
      else {
          mate_unmapped = 1;
      }
  }
  if(0 < sams->n) {
      for(i=0;i<sams->n;i++) {
          tmap_map_sam_print(seq, refseq, &sams->sams[i], sam_sff_tags, sams->max, i, end_num, mate_unmapped, mate);
      }
  }
  else {
      tmap_map_sam_print(seq, refseq, NULL, sam_sff_tags, sams->max, 0, end_num, mate_unmapped, mate);
  }
}

void
tmap_map_sams_filter1(tmap_map_sams_t *sams, int32_t aln_output_mode, int32_t algo_id)
{
  int32_t i, j, k;
  int32_t n_best = 0;
  int32_t best_score, cur_score;
  int32_t best_subo;

  if(sams->n <= 1) {
      return;
  }

  for(i=j=0;i<sams->n;i++) {
      if(TMAP_MAP_ALGO_NONE == algo_id
         || sams->sams[i].algo_id == algo_id) {
          j++;
      }
  }
  if(j <= 1) {
      return;
  }

  best_score = best_subo = INT32_MIN;
  n_best = 0;
  for(i=0;i<sams->n;i++) {
      if(TMAP_MAP_ALGO_NONE == algo_id
         || sams->sams[i].algo_id == algo_id) {
          cur_score = sams->sams[i].score;
          if(best_score < cur_score) {
              if(0 < n_best) {
                  best_subo = best_score;
              }
              best_score = cur_score;
              n_best = 1;
          }
          else if(!(cur_score < best_score)) { // equal
              best_subo = best_score; // more than one mapping
              n_best++;
          }
          else if(best_subo < cur_score) {
              best_subo = cur_score;
          }
          // check sub-optimal
          if(TMAP_MAP_ALGO_MAP2 == sams->sams[i].algo_id
             || TMAP_MAP_ALGO_MAP3 == sams->sams[i].algo_id) {
              cur_score = sams->sams[i].score_subo;
              if(best_subo < cur_score) {
                  best_subo = cur_score;
              }
          }

      }
  }

  // adjust mapping qualities
  if(1 < n_best) {
      for(i=0;i<sams->n;i++) {
          if(TMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              sams->sams[i].mapq = 0;
          }
      }
  }
  else {
      for(i=0;i<sams->n;i++) {
          if(TMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              cur_score = sams->sams[i].score;
              if(cur_score < best_score) { // not the best
                  sams->sams[i].mapq = 0;
              }
          }
      }
  }

  // adjust suboptimal
  if(TMAP_MAP_ALGO_NONE == algo_id) {
      for(i=0;i<sams->n;i++) {
          sams->sams[i].score_subo = best_subo;
      }
  }

  if(TMAP_MAP_OPT_ALN_MODE_ALL == aln_output_mode) {
      return;
  }

  // copy to the front
  if(n_best < sams->n) {
      for(i=j=0;i<sams->n;i++) {
          if(TMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              cur_score = sams->sams[i].score;
              if(cur_score < best_score) { // not the best
                  tmap_map_sam_destroy(&sams->sams[i]);
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                  }
                  j++;
              }
          }
          else {
              if(j < i) { // copy if we are not on the same index
                  tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
              }
              j++;
          }
      }
      // reallocate
      tmap_map_sams_realloc(sams, j);
  }

  if(TMAP_MAP_OPT_ALN_MODE_UNIQ_BEST == aln_output_mode) {
      if(1 < n_best) { // there can only be one
          if(TMAP_MAP_ALGO_NONE == algo_id) {
              tmap_map_sams_realloc(sams, 0);
          }
          else {
              // get rid of all of them
              for(i=j=0;i<sams->n;i++) {
                  if(sams->sams[i].algo_id == algo_id) {
                      tmap_map_sam_destroy(&sams->sams[i]);
                  }
                  else {
                      if(j < i) { // copy if we are not on the same index
                          tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                      }
                      j++;
                  }
              }
              tmap_map_sams_realloc(sams, j);
          }
      }
  }
  else if(TMAP_MAP_OPT_ALN_MODE_RAND_BEST == aln_output_mode) { // get a random
      int32_t r = (int32_t)(drand48() * n_best);

      // keep the rth one
      if(TMAP_MAP_ALGO_NONE == algo_id) {
          if(0 != r) {
              tmap_map_sam_destroy(&sams->sams[0]);
              tmap_map_sam_copy_and_nullify(&sams->sams[0], &sams->sams[r]);
          }
          // reallocate
          tmap_map_sams_realloc(sams, 1);
      }
      else {
          // keep the rth one
          for(i=j=k=0;i<sams->n;i++) {
              if(sams->sams[i].algo_id == algo_id) {
                  if(k == r) { // keep
                      if(j < i) { // copy if we are not on the same index
                          tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                      }
                      j++;
                  }
                  else { // free
                      tmap_map_sam_destroy(&sams->sams[i]);
                  }
                  k++;
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                  }
                  j++;
              }
          }
          tmap_map_sams_realloc(sams, j);
      }
  }
  else if(TMAP_MAP_OPT_ALN_MODE_ALL_BEST == aln_output_mode) {
      // do nothing
  }
  else {
      tmap_error("bug encountered", Exit, OutOfRange);
  }
}

void
tmap_map_sams_filter2(tmap_map_sams_t *sams, int32_t score_thr, int32_t mapq_thr)
{
  int32_t i, j;

  // filter based on score and mapping quality
  for(i=j=0;i<sams->n;i++) {
      if(sams->sams[i].score < score_thr || sams->sams[i].mapq < mapq_thr) {
          tmap_map_sam_destroy(&sams->sams[i]);
      }
      else {
          if(j < i) { // copy if we are not on the same index
              tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
          }
          j++;
      }
  }
  tmap_map_sams_realloc(sams, j);
}

void
tmap_map_sams_filter(tmap_map_sams_t *sams, int32_t aln_output_mode)
{
  tmap_map_sams_filter1(sams, aln_output_mode, TMAP_MAP_ALGO_NONE);
}

void
tmap_map_util_remove_duplicates(tmap_map_sams_t *sams, int32_t dup_window)
{
  int32_t i, next_i, j, k, end, best_score_i, best_score_n, best_score_subo;

  if(dup_window < 0 || sams->n <= 0) {
      return;
  }

  // sort
  // NB: since tmap_map_util_sw_gen_score only sets the end position of the
  // alignment, use that
  tmap_sort_introsort(tmap_map_sam_sort_coord_end, sams->n, sams->sams);

  // remove duplicates within a window
  for(i=j=0;i<sams->n;) {

      // get the change
      end = best_score_i = i;
      best_score_n = 0;
      best_score_subo = INT32_MIN;
      while(end+1 < sams->n) {
          if(sams->sams[end].seqid == sams->sams[end+1].seqid
             && sams->sams[end].strand == sams->sams[end+1].strand
             && fabs((sams->sams[end].pos + sams->sams[end].target_len) - (sams->sams[end+1].pos + sams->sams[end+1].target_len)) <= dup_window) {
              // track the best scoring
              if(sams->sams[best_score_i].score == sams->sams[end+1].score) {
                  best_score_i = end+1;
                  best_score_n++;
              }
              else if(sams->sams[best_score_i].score < sams->sams[end+1].score) {
                  best_score_i = end+1;
                  best_score_n = 1;
              }
              if(best_score_subo < sams->sams[end+1].score_subo) {
                  best_score_subo = sams->sams[end+1].score_subo;
              }
              end++;
          }
          else {
              break;
          }
      }
      next_i = end+1;

      // randomize the best scoring      
      if(1 < best_score_n) {
          k = (int32_t)(best_score_n * drand48()); // make this zero-based 
          best_score_n = 0; // make this one-based
          end = i;
          while(best_score_n <= k) { // this assumes we know there are at least "best_score
              if(sams->sams[best_score_i].score == sams->sams[end].score) {
                  best_score_i = end;
                  best_score_n++;
              }
              end++;
          }
      }

      // copy over the best
      if(j != best_score_i) {
          // destroy
          tmap_map_sam_destroy(&sams->sams[j]);
          // nullify
          tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[best_score_i]);
      }

      // copy over sub-optimal score
      sams->sams[j].score_subo = best_score_subo;

      // next
      i = next_i;
      j++;
  }

  // resize
  tmap_map_sams_realloc(sams, j);
}

inline void
tmap_map_util_mapq(tmap_map_sams_t *sams, int32_t seq_len, tmap_map_opt_t *opt)
{
  int32_t i, best_i;
  int32_t n_best = 0, n_best_subo = 0;
  int32_t best_score, cur_score, best_subo, best_subo2;
  int32_t mapq;
  int32_t stage = -1;
  int32_t algo_id = TMAP_MAP_ALGO_NONE; 
  int32_t best_repetitive = 0;

  // estimate mapping quality TODO: this needs to be refined
  best_i = 0;
  best_score = INT32_MIN;
  best_subo = best_subo2 = opt->score_thr;
  n_best = n_best_subo = 0;
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      if(best_score < cur_score) {
          // save sub-optimal
          best_subo = best_score;
          n_best_subo = n_best;
          // update
          best_score = cur_score;
          n_best = 1;
          best_i = i;
          stage = (algo_id == TMAP_MAP_ALGO_NONE) ? sams->sams[i].algo_stage-1 : -1;
          algo_id = (algo_id == TMAP_MAP_ALGO_NONE) ? sams->sams[i].algo_id : -1;
          if(sams->sams[i].algo_id == TMAP_MAP_ALGO_MAP2 && 1 == sams->sams[i].aux.map2_aux->flag) {
              best_repetitive = 1;
          }
          else {
              best_repetitive = 0;
          }
      }
      else if(cur_score == best_score) { // qual
          if(sams->sams[i].algo_id == TMAP_MAP_ALGO_MAP2 && 1 == sams->sams[i].aux.map2_aux->flag) {
              best_repetitive = 1;
          }
          n_best++;
      }
      else {
          if(best_subo < cur_score) {
              best_subo = cur_score;
              n_best_subo = 1;
          }
          else if(best_subo == cur_score) {
              n_best_subo++;
          }
      }
      // get the best subo-optimal score
      cur_score = sams->sams[i].score_subo;
      if(INT32_MIN == cur_score) {
          // ignore
      }
      else if(best_subo < cur_score) {
          best_subo2 = cur_score;
      }
  }
  if(best_subo < best_subo2) best_subo = best_subo2;
  if(1 < n_best || best_score <= best_subo || 0 < best_repetitive) {
      mapq = 0;
  }
  else {
      if(0 == n_best_subo) {
          n_best_subo = 1;
          best_subo = opt->score_thr; 
      }
      /*
         fprintf(stderr, "n_best=%d n_best_subo=%d\n",
         n_best, n_best_subo);
         fprintf(stderr, "best_score=%d best_subo=%d\n",
         best_score, best_subo);
         */
      // Note: this is the old calculationg, based on BWA-long
      //mapq = (int32_t)((n_best / (1.0 * n_best_subo)) * (best_score - best_subo) * (250.0 / best_score + 0.03 / opt->score_match) + .499);
      //
      double sf = 0.4; // initial scaling factor.  Note: 250 * sf is the maximum mapping quality.
      sf *= 250.0 / ((double)opt->score_match * seq_len); // scale based on the best possible alignment score
      sf *= (n_best / (1.0 * n_best_subo)); // scale based on number of sub-optimal mappings
      sf *= (double)(best_score - best_subo + 1 ); // scale based on distance to the sub-optimal mapping
      //sf *= (seq_len < 10) ? 1.0 : log10(seq_len); // scale based on longer reads having more information content
      mapq = (int32_t)(sf + 0.99999);
      if(mapq > 250) mapq = 250;
      if(mapq <= 0) mapq = 1;
  }
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      if(cur_score == best_score) {
          sams->sams[i].mapq = mapq;
      }
      else {
          sams->sams[i].mapq = 0;
      }
  }
}

#define __map_util_gen_ap(par, opt) do { \
    int32_t i; \
    for(i=0;i<25;i++) { \
        (par).matrix[i] = -(opt)->pen_mm; \
    } \
    for(i=0;i<4;i++) { \
        (par).matrix[i*5+i] = (opt)->score_match; \
    } \
    (par).gap_open = (opt)->pen_gapo; (par).gap_ext = (opt)->pen_gape; \
    (par).gap_end = (opt)->pen_gape; \
    (par).row = 5; \
    (par).band_width = (opt)->bw; \
} while(0)

// ACGTNBDHKMRSVWYN
// This gives a mask for match iupac characters versus A/C/G/T/N
//  A  C  G  T   N  B  D  H   K  M  R  S   V  W  Y  N
static int32_t matrix_iupac_mask[89] = {
    1, 0, 0, 0,  1, 0, 1, 1,  0, 1, 1, 0,  1, 1, 0, 1, // A
    0, 1, 0, 0,  1, 1, 0, 1,  0, 1, 0, 1,  1, 0, 1, 1, // C
    0, 0, 1, 0,  1, 1, 1, 0,  1, 0, 1, 1,  1, 0, 0, 1, // G
    0, 0, 0, 1,  1, 1, 1, 1,  1, 0, 0, 0,  0, 1, 1, 1, // T 
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1  // N
};

#define __map_util_gen_ap_iupac(par, opt) do { \
    int32_t i; \
    for(i=0;i<80;i++) { \
        if(0 < matrix_iupac_mask[i]) (par).matrix[i] = (opt)->score_match; \
        else (par).matrix[i] = -(opt)->pen_mm; \
    } \
    (par).gap_open = (opt)->pen_gapo; (par).gap_ext = (opt)->pen_gape; \
    (par).gap_end = (opt)->pen_gape; \
    (par).row = 16; \
    (par).band_width = (opt)->bw; \
} while(0)

tmap_map_sams_t *
tmap_map_util_sw_gen_score(tmap_refseq_t *refseq,
                 tmap_map_sams_t *sams, 
                 tmap_seq_t *seq,
                 tmap_map_opt_t *opt)
{
  int32_t i;
  int32_t start, end;
  tmap_map_sams_t *sams_tmp = NULL;
  tmap_seq_t *seqs[2] = {NULL, NULL};
  int32_t seq_len=0, tlen, target_mem=0;
  uint8_t *target=NULL;
  int32_t best_subo;
  tmap_vsw_query_t *vsw_query[2] = {NULL, NULL};
  tmap_vsw_opt_t *vsw_opt = NULL;
  uint32_t start_pos, end_pos;
  int32_t overflow, softclip_start, softclip_end;

  if(0 == sams->n) {
      return sams;
  }

  // the final mappings will go here 
  sams_tmp = tmap_map_sams_init(sams);
  tmap_map_sams_realloc(sams_tmp, sams->n);

  // sort by strand/chr/pos/score
  tmap_sort_introsort(tmap_map_sam_sort_coord, sams->n, sams->sams);
  softclip_start = (TMAP_MAP_OPT_SOFT_CLIP_LEFT == opt->softclip_type || TMAP_MAP_OPT_SOFT_CLIP_ALL == opt->softclip_type) ? 1 : 0;
  softclip_end = (TMAP_MAP_OPT_SOFT_CLIP_RIGHT == opt->softclip_type || TMAP_MAP_OPT_SOFT_CLIP_ALL == opt->softclip_type) ? 1 : 0;

  // initialize opt
  vsw_opt = tmap_vsw_opt_init(opt->score_match, opt->pen_mm, opt->pen_gapo, opt->pen_gape, opt->score_thr);

  // init seqs
  seq_len = tmap_seq_get_bases(seq)->l;
  seqs[0] = tmap_seq_clone(seq); // clone
  seqs[1] = tmap_seq_clone(seq); // clone
  tmap_seq_reverse_compliment(seqs[1]);
  tmap_seq_to_int(seqs[0]);
  tmap_seq_to_int(seqs[1]);

  vsw_query[0] = tmap_vsw_query_init((uint8_t*)tmap_seq_get_bases(seqs[0])->s, seq_len, 
                                     seq_len, softclip_start, softclip_end, vsw_opt);
  vsw_query[1] = tmap_vsw_query_init((uint8_t*)tmap_seq_get_bases(seqs[1])->s, seq_len, 
                                     seq_len, softclip_start, softclip_end, vsw_opt);

  i = start = end = 0;
  best_subo = INT32_MIN;
  start_pos = end_pos = 0;
  while(end < sams->n) {
      uint8_t strand, *query=NULL;
      uint32_t qlen;
      tmap_map_sam_t tmp_sam;

      // get the strand/start/end positions
      strand = sams->sams[end].strand;
      if(start == end) {
          start_pos = sams->sams[start].pos + 1; 
          end_pos = sams->sams[start].pos + sams->sams[start].target_len; 
      }

      // sub-optimal score
      if(best_subo < sams->sams[end].score_subo) {
          best_subo = sams->sams[end].score_subo;
      }

      // check if the hits can be banded
      if(end + 1 < sams->n) {
          if(sams->sams[end].strand == sams->sams[end+1].strand 
             && sams->sams[end].seqid == sams->sams[end+1].seqid
             && sams->sams[end+1].pos - sams->sams[end].pos <= opt->max_seed_band) {
              end++;
              if(end_pos < sams->sams[end].pos + sams->sams[end].target_len) {
                  end_pos = sams->sams[end].pos + sams->sams[end].target_len; // one-based
              }
              continue; // there may be more to add
          }
      }

      // choose a random one within the window
      if(start == end) {
          tmp_sam = sams->sams[start];
      }
      else {
          int32_t r = (int32_t)(drand48() * (end - start + 1));
          r += start;
          tmp_sam = sams->sams[r];
      }
      
      // update the query sequence
      query = (uint8_t*)tmap_seq_get_bases(seqs[strand])->s;
      qlen = tmap_seq_get_bases(seqs[strand])->l;

      // add in band width
      // one-based
      if(start_pos < opt->bw) {
          start_pos = 1;
      }
      else {
          start_pos -= opt->bw - 1;
      }
      end_pos += opt->bw - 1;
      if(refseq->annos[sams->sams[end].seqid].len < end_pos) {
          end_pos = refseq->annos[sams->sams[end].seqid].len; // one-based
      }

      // get the target sequence
      tlen = end_pos - start_pos + 1;
      if(target_mem < tlen) { // more memory?
          target_mem = tlen;
          tmap_roundup32(target_mem);
          target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }
      // NB: IUPAC codes are turned into mismatches
      if(NULL == tmap_refseq_subseq2(refseq, sams->sams[end].seqid+1, start_pos, end_pos, target, 1, NULL)) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }


      // Debugging
#ifdef TMAP_VSW_DEBUG
      int j;
      fprintf(stderr, "seqid:%u start_pos=%u end_pos=%u strand=%d\n", sams->sams[end].seqid+1, start_pos, end_pos, strand);
      for(j=0;j<qlen;j++) {
          fputc("ACGTN"[query[j]], stderr);
      }
      fputc('\n', stderr);
      for(j=0;j<tlen;j++) {
          fputc("ACGTN"[target[j]], stderr);
      }
      fputc('\n', stderr);
#endif

      // initialize the bounds
      tmp_sam.query_start = tmp_sam.query_end = 0;
      tmp_sam.target_start = tmp_sam.target_end = 0;

      // NB: if match/mismatch penalties are on the opposite strands, we may
      // have wrong scores
      if(0 == strand) {
          tmp_sam.score = tmap_vsw_sse2(vsw_query[strand], query, qlen,
                                        target, tlen, 
                                        softclip_start, softclip_end,
                                        vsw_opt, 
                                        &tmp_sam.score_fwd, &tmp_sam.score_rev,
                                        &tmp_sam.query_start, &tmp_sam.query_end,
                                        &tmp_sam.target_start, &tmp_sam.target_end,
                                        &overflow, opt->score_thr, 0);
      }
      else {
          tmp_sam.score = tmap_vsw_sse2(vsw_query[strand], query, qlen,
                                        target, tlen, 
                                        softclip_end, softclip_start,
                                        vsw_opt, 
                                        &tmp_sam.score_fwd, &tmp_sam.score_rev,
                                        &tmp_sam.query_start, &tmp_sam.query_end,
                                        &tmp_sam.target_start, &tmp_sam.target_end,
                                        &overflow, opt->score_thr, 0);
      }
      if(1 == overflow) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }

      if(opt->score_thr <= tmp_sam.score) {
          tmap_map_sam_t *s = &sams_tmp->sams[i];
          // shallow copy previous data 
          (*s) = tmp_sam; 

          // nullify the cigar
          s->n_cigar = 0;
          s->cigar = NULL;

          // adjust target length and position NB: query length is implicitly
          // stored in s->query_end (consider on the next pass)
          s->pos = start_pos - 1; // zero-based
          s->target_len = s->target_end + 1;
          /*
             fprintf(stderr, "strand=%d\n", strand);
             fprintf(stderr, "%d-%d %d-%d %d\n",
             s->query_start,
             s->query_end,
             s->target_start,
             s->target_end,
             s->target_len);
             */

          // # of seeds
          s->n_seeds = (end - start + 1);

          // update aux data
          tmap_map_sam_malloc_aux(s, s->algo_id);
          switch(s->algo_id) {
            case TMAP_MAP_ALGO_MAP1:
              (*s->aux.map1_aux) = (*tmp_sam.aux.map1_aux);
              break;
            case TMAP_MAP_ALGO_MAP2:
              (*s->aux.map2_aux) = (*tmp_sam.aux.map2_aux);
              break;
            case TMAP_MAP_ALGO_MAP3:
              (*s->aux.map3_aux) = (*tmp_sam.aux.map3_aux);
              break;
            case TMAP_MAP_ALGO_MAPVSW:
              (*s->aux.map_vsw_aux) = (*tmp_sam.aux.map_vsw_aux);
              break;
            default:
              tmap_error("bug encountered", Exit, OutOfRange);
              break;
          }

          i++;
      }
      else {
          tmp_sam.score = INT32_MIN;
      }

      // update start/end
      end++;
      start = end;
  }

  // realloc
  tmap_map_sams_realloc(sams_tmp, i);

  // update the sub-optimal
  for(i=0;i<sams_tmp->n;i++) {
      sams_tmp->sams[i].score_subo = best_subo;
  }

  // free memory
  for(i=0;i<2;i++) {
      if(NULL != seqs[i]) {
          tmap_seq_destroy(seqs[i]);
      }
  }
  tmap_map_sams_destroy(sams);
  free(target);
  tmap_vsw_opt_destroy(vsw_opt);
  tmap_vsw_query_destroy(vsw_query[0]);
  tmap_vsw_query_destroy(vsw_query[1]);

  return sams_tmp;
}

static void tmap_map_util_keytrim(uint8_t *query, int32_t qlen,
                                  uint8_t *target, int32_t tlen,
                                  int8_t strand, uint8_t key_base, 
                                  tmap_map_sam_t *s)
{
  int32_t j, k, l;
  int32_t op, op_len;
  //NB: we may only need to look at the first cigar
  op = op_len = 0;
  if(0 == strand) { // forward
      for(j = k = l = 0; j < s->n_cigar; j++) {
          // get the cigar
          op = TMAP_SW_CIGAR_OP(s->cigar[j]);
          op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[j]);
          if(op == BAM_CDEL) break; // looking for mismatch/insertion
          if(op == BAM_CSOFT_CLIP) break; // already start trimming

          while(0 < op_len) {
              if(query[k] != key_base) break;
              if(op == BAM_CINS) {
              }
              else if(op == BAM_CMATCH && target[l] != query[k]) {
                  l++;
              }
              else {
                  break;
              }
              op_len--;
              k++; // since we can only have match/mismatch/insertion
          }
          if(0 == k) { 
              // no trimming
              break;
          }
          else if(0 < op_len) {
              if(j == 0) {
                  // reallocate
                  s->n_cigar++;
                  s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*s->n_cigar, "s->cigar");
                  for(k=s->n_cigar-1;0<k;k--) { // shift up
                      s->cigar[k] = s->cigar[k-1];
                  }
                  s->cigar[0] = 0;
                  j++; // reflect the shift in j 
              }
              // NB: 0 < j
              // add to the leading soft-clip
              TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[j]) - op_len + TMAP_SW_CIGAR_LENGTH(s->cigar[0]));
              // reduce the cigar length
              TMAP_SW_CIGAR_STORE(s->cigar[j], op, op_len);
              break;
          }
          else { // NB: the full op was removed
              if(0 == j) {
                  TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[0]));
              }
              else {
                  // add to the leading soft-clip
                  TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[j]) + TMAP_SW_CIGAR_LENGTH(s->cigar[0]));
                  // shift down and overwrite the current cigar
                  for(k=j;k<s->n_cigar-1;k++) {
                      s->cigar[k] = s->cigar[k+1];
                  }
                  s->n_cigar--;
                  j--; // since j is incremented later
              }
          }
      }
      // update the position based on the number of reference bases we
      // skipped
      s->pos += l;
  }
  else { // reverse
      for(j = s->n_cigar-1, k = qlen-1, l = tlen-1; 0 <= j; j--) {
          // get the cigar
          op = TMAP_SW_CIGAR_OP(s->cigar[j]);
          op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[j]);
          if(op == BAM_CDEL) break; // looking for mismatch/insertion
          if(op == BAM_CSOFT_CLIP) break; // already start trimming
          while(0 < op_len) {
              if(query[k] != (3 - key_base)) break;
              if(op == BAM_CINS) {
              }
              else if(op == BAM_CMATCH && target[l] != query[k]) {
                  l--;
              }
              else {
                  break;
              }
              op_len--;
              k--; // since we can only have match/mismatch/insertion
          }
          if(qlen-1 == k) { 
              // no trimming
              break;
          }
          else if(0 < op_len) {
              if(j == s->n_cigar-1) {
                  // reallocate
                  s->n_cigar++;
                  s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*s->n_cigar, "s->cigar");
                  s->cigar[s->n_cigar-1] = 0;
              }
              // add to the ending soft-clip
              TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar-1], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[j]) - op_len + TMAP_SW_CIGAR_LENGTH(s->cigar[s->n_cigar-1]));
              // reduce the cigar length
              TMAP_SW_CIGAR_STORE(s->cigar[j], op, op_len);
              break;
          }
          else { // NB: the full op was removed
              if(j == s->n_cigar-1) {
                  TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar-1], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[s->n_cigar-1]));
              }
              else {
                  // add to the ending soft-clip
                  TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar-1], BAM_CSOFT_CLIP, TMAP_SW_CIGAR_LENGTH(s->cigar[s->n_cigar-1]));
                  // shift down and overwrite the current cigar
                  for(k=j;k<s->n_cigar-1;k++) {
                      s->cigar[k] = s->cigar[k+1];
                  }
                  s->n_cigar--;
                  j++; // since j is decremented later
              }
          }
      }
  }
}

tmap_map_sams_t *
tmap_map_util_sw_gen_cigar(tmap_refseq_t *refseq,
                 tmap_map_sams_t *sams, 
                 tmap_seq_t *seq,
                 tmap_map_opt_t *opt)
{
  int32_t i, j, matrix[25], matrix_iupac[80];
  int32_t start, end;
  tmap_map_sams_t *sams_tmp = NULL;
  tmap_sw_param_t par, par_iupac;
  tmap_sw_path_t *path = NULL;
  int32_t path_len, path_mem=0;
  tmap_seq_t *seqs[2] = {NULL, NULL};
  int32_t seq_len=0, tlen, target_mem=0;
  uint8_t *target=NULL;
  tmap_vsw_query_t *vsw_query[2] = {NULL, NULL};
  tmap_vsw_opt_t *vsw_opt = NULL;
  uint32_t start_pos, end_pos;
  int32_t overflow, softclip_start, softclip_end;
  uint8_t key_base = 0;
  int32_t iupac_init = 0;

  if(0 == sams->n) {
      return sams;
  }

  // the final mappings will go here 
  sams_tmp = tmap_map_sams_init(sams);
  tmap_map_sams_realloc(sams_tmp, sams->n);

  // scoring matrix
  par.matrix = matrix;
  __map_util_gen_ap(par, opt); 

  // sort by strand/chr/pos/score
  tmap_sort_introsort(tmap_map_sam_sort_coord, sams->n, sams->sams);
  softclip_start = (TMAP_MAP_OPT_SOFT_CLIP_LEFT == opt->softclip_type || TMAP_MAP_OPT_SOFT_CLIP_ALL == opt->softclip_type) ? 1 : 0;
  softclip_end = (TMAP_MAP_OPT_SOFT_CLIP_RIGHT == opt->softclip_type || TMAP_MAP_OPT_SOFT_CLIP_ALL == opt->softclip_type) ? 1 : 0;
  
  // initialize opt
  vsw_opt = tmap_vsw_opt_init(opt->score_match, opt->pen_mm, opt->pen_gapo, opt->pen_gape, opt->score_thr);

  // init seqs
  seq_len = tmap_seq_get_bases(seq)->l;
  seqs[0] = tmap_seq_clone(seq); // clone
  seqs[1] = tmap_seq_clone(seq); // clone
  tmap_seq_reverse_compliment(seqs[1]);
  tmap_seq_to_int(seqs[0]);
  tmap_seq_to_int(seqs[1]);

  vsw_query[0] = tmap_vsw_query_init((uint8_t*)tmap_seq_get_bases(seqs[0])->s, seq_len, 
                                     seq_len, softclip_start, softclip_end, vsw_opt);
  vsw_query[1] = tmap_vsw_query_init((uint8_t*)tmap_seq_get_bases(seqs[1])->s, seq_len, 
                                     seq_len, softclip_start, softclip_end, vsw_opt);
      
  if(1 == opt->softclip_key) {
      // get the last base
      if(NULL != opt->key_seq) {
          key_base = tmap_nt_char_to_int[(int)opt->key_seq[strlen(opt->key_seq)-1]]; 
      }
      else {
          if(TMAP_SEQ_TYPE_SFF != seq->type) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }
          key_base = tmap_nt_char_to_int[(int)seq->data.sff->gheader->key->s[seq->data.sff->gheader->key->l-1]];
      }
  }
      
  i = start = end = 0;
  start_pos = end_pos = 0;
  while(end < sams->n) {
      uint8_t strand, *query=NULL, *query_rc=NULL, *tmp_target;
      uint32_t qlen;
      tmap_map_sam_t tmp_sam;
      tmap_map_sam_t *s = NULL;
      int32_t query_start, query_end;
      int32_t conv = 0;
      
      // do not band when generating the cigar
      tmp_sam = sams->sams[end];

      // get the strand/start/end positions
      strand = tmp_sam.strand;
      start_pos = tmp_sam.pos + 1;
      end_pos = tmp_sam.pos + tmp_sam.target_len;

      // retrieve the reverse complimented query sequence
      query_rc  = (uint8_t*)tmap_seq_get_bases(seqs[1-strand])->s;
      query_rc += seq_len - tmp_sam.query_end - 1; // offset query
      qlen = tmp_sam.query_end + 1; // adjust based on the query end

      // get the target sequence
      tlen = tmp_sam.target_end + 1; // adjust based on the target end
      if(target_mem < tlen) { // more memory?
          target_mem = tlen;
          tmap_roundup32(target_mem);
          target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }
      // NB: IUPAC codes are turned into mismatches
      if(NULL == tmap_refseq_subseq2(refseq, sams->sams[end].seqid+1, start_pos, end_pos, target, 1, &conv)) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }
      
      /**
       * Step 1: find the start of the alignment
       */

      // NB: if match/mismatch penalties are on the opposite strands, we may
      // have wrong scores
      if(0 == strand) {
          tmp_sam.score = tmap_vsw_sse2(vsw_query[strand], query_rc, qlen,
                                        target, tlen, 
                                        softclip_start, softclip_end,
                                        vsw_opt, 
                                        &tmp_sam.score_fwd, &tmp_sam.score_rev,
                                        &tmp_sam.query_start, &tmp_sam.query_end,
                                        &tmp_sam.target_start, &tmp_sam.target_end,
                                        &overflow, opt->score_thr, 1);
      }
      else {
          tmp_sam.score = tmap_vsw_sse2(vsw_query[strand], query_rc, qlen,
                                        target, tlen, 
                                        softclip_end, softclip_start,
                                        vsw_opt, 
                                        &tmp_sam.score_fwd, &tmp_sam.score_rev,
                                        &tmp_sam.query_start, &tmp_sam.query_end,
                                        &tmp_sam.target_start, &tmp_sam.target_end,
                                        &overflow, opt->score_thr, 1);
      }
      if(1 == overflow) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }

      /*
      fprintf(stderr, "tmp_sam.query_start=%d\n", tmp_sam.query_start);
      fprintf(stderr, "tmp_sam.query_end=%d\n", tmp_sam.query_end);
      fprintf(stderr, "tmp_sam.target_start=%d\n", tmp_sam.target_start);
      fprintf(stderr, "tmp_sam.target_end=%d\n", tmp_sam.target_end);
      fprintf(stderr, "tmp_sam.pos=%u\n", tmp_sam.pos);
      */

      // adjust the query and target based off of the start of the alignment
      query = (uint8_t*)tmap_seq_get_bases(seqs[strand])->s;
      query += tmp_sam.query_start; // offset query
      qlen = tmp_sam.query_end - tmp_sam.query_start + 1; // update query length
      tmp_target = target;
      target += tmp_sam.target_start;
      tlen = tmp_sam.target_end - tmp_sam.target_start + 1;
      tmp_sam.pos += tmp_sam.target_start; // keep it zero based

      /**
       * Step 2: generate the cigar
       */

      if(0 < conv) { // NB: there were IUPAC bases
          // init iupac parameters, if necessary
          if(0 == iupac_init) {
              par_iupac.matrix = matrix_iupac;
              __map_util_gen_ap_iupac(par_iupac, opt);
              iupac_init = 1;
          }
          // Get the new target
          // NB: IUPAC codes are turned into mismatches
          start_pos += tmp_sam.target_start;
          if(NULL == tmap_refseq_subseq2(refseq, sams->sams[end].seqid+1, start_pos, end_pos, target, 0, NULL)) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }
      }

      // path memory
      if(path_mem <= tlen + seq_len) { // lengthen the path
          path_mem = tlen + seq_len;
          tmap_roundup32(path_mem);
          path = tmap_realloc(path, sizeof(tmap_sw_path_t)*path_mem, "path");
      }

      // Debugging
      /*
      for(j=0;j<qlen;j++) {
          fputc("ACGTN"[query[j]], stderr);
      }
      fputc('\n', stderr);
      for(j=0;j<tlen;j++) {
          fputc(tmap_iupac_int_to_char[target[j]], stderr);
      }
      fputc('\n', stderr);
      */


      // path memory
      if(path_mem <= tlen + qlen) { // lengthen the path
          path_mem = tlen + qlen;
          tmap_roundup32(path_mem);
          path = tmap_realloc(path, sizeof(tmap_sw_path_t)*path_mem, "path");
      }

      query_start = tmp_sam.query_start;
      query_end = tmp_sam.query_end;
      s = &sams_tmp->sams[i];
      
      // shallow copy previous data 
      (*s) = tmp_sam; 

      // Smith Waterman with banding
      // NB: we store the score from the banded version, which does not allow
      // ins then del, or del then ins.  The vectorized version does.
      // NB: left genomic indel justification is facilitated by always using the
      // forward strand target/query combination.
      // NB: iupac bases may also increase the score
      if(0 < conv) { // NB: there were IUPAC bases
          s->score = tmap_sw_global_banded_core(target, tlen, query, qlen, &par_iupac,
                                                tmp_sam.score_fwd, path, &path_len, 0); 
      }
      else {
          s->score = tmap_sw_global_banded_core(target, tlen, query, qlen, &par,
                                                tmp_sam.score_fwd, path, &path_len, 0); 
      }

      s->pos = s->pos + (path[path_len-1].i-1); // zero-based 
      if(path[path_len-1].ctype == TMAP_SW_FROM_I) {
          s->pos++;
      }
      s->cigar = tmap_sw_path2cigar(path, path_len, &s->n_cigar);
      if(0 == s->n_cigar) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }

      s->target_len = 0;
      for(j=0;j<s->n_cigar;j++) {
          switch(TMAP_SW_CIGAR_OP(s->cigar[j])) {
            case BAM_CMATCH:
            case BAM_CDEL:
              s->target_len += TMAP_SW_CIGAR_LENGTH(s->cigar[j]);
              break;
            default:
              break;
          }
      }

      // add soft clipping 
      if(0 < query_start) {
          // soft clip the front of the read
          s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1+s->n_cigar), "s->cigar");
          for(j=s->n_cigar-1;0<=j;j--) { // shift up
              s->cigar[j+1] = s->cigar[j];
          }
          TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, query_start);
          s->n_cigar++;
      }
      if(query_end < seq_len-1) {
          // soft clip the end of the read
          s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1+s->n_cigar), "s->cigar");
          TMAP_SW_CIGAR_STORE(s->cigar[s->n_cigar], BAM_CSOFT_CLIP, seq_len - query_end - 1);
          s->n_cigar++;
      }

      // update aux data
      tmap_map_sam_malloc_aux(s, s->algo_id);
      switch(s->algo_id) {
        case TMAP_MAP_ALGO_MAP1:
          (*s->aux.map1_aux) = (*tmp_sam.aux.map1_aux);
          break;
        case TMAP_MAP_ALGO_MAP2:
          (*s->aux.map2_aux) = (*tmp_sam.aux.map2_aux);
          break;
        case TMAP_MAP_ALGO_MAP3:
          (*s->aux.map3_aux) = (*tmp_sam.aux.map3_aux);
          break;
        case TMAP_MAP_ALGO_MAPVSW:
          (*s->aux.map_vsw_aux) = (*tmp_sam.aux.map_vsw_aux);
          break;
        default:
          tmap_error("bug encountered", Exit, OutOfRange);
          break;
      }

      // key trim the data
      if(1 == opt->softclip_key) {
          tmap_map_util_keytrim(query, qlen, target, tlen, strand, key_base, s);
      }

      i++;

      // update start/end
      end++;
      start = end;
      target = tmp_target;
  }
      
  // realloc
  tmap_map_sams_realloc(sams_tmp, i);

  // free memory
  for(i=0;i<2;i++) {
      if(NULL != seqs[i]) {
          tmap_seq_destroy(seqs[i]);
      }
  }
  tmap_map_sams_destroy(sams);
  free(path);
  free(target);
  tmap_vsw_query_destroy(vsw_query[0]);
  tmap_vsw_query_destroy(vsw_query[1]);
  tmap_vsw_opt_destroy(vsw_opt);

  return sams_tmp;
}

#define _tmap_map_util_sw_aux_path_adjust(_p, _ql, _tl) do { \
    (_p).i = _tl - (_p).i + 1; \
    (_p).j = _ql - (_p).j + 1; \
} while(0) \

int32_t 
tmap_map_util_sw_aux(tmap_map_sam_t *sam,
                 uint8_t *target, int32_t target_length,
                 uint8_t *query, int32_t query_length,
                 uint32_t seqid, uint32_t pos,
                 tmap_sw_param_t *par, tmap_sw_path_t *path, int32_t *path_len,
                 int32_t score_thr, int32_t softclip_type, int32_t strand)
{
  int32_t i, score, score_subo;

  if(1 == strand) { // reverse, but do not compliment, since we want consistent behavior on the nucleotides
      tmap_reverse_int(query, query_length);
      tmap_reverse_int(target, target_length);
  }

  /*
  fprintf(stderr, "strand=%d\n", strand);
  for(i=0;i<query_length;i++)
    fputc("ACGTN"[query[i]], stderr);
  fputc('\n', stderr);
  for(i=0;i<target_length;i++)
    fputc("ACGTN"[target[i]], stderr);
  fputc('\n', stderr);
  */

  switch(softclip_type) {
    case TMAP_MAP_OPT_SOFT_CLIP_ALL:
      //fprintf(stderr, "TMAP_MAP_OPT_SOFT_CLIP_ALL\n");
      score = tmap_sw_clipping_core(target, target_length, query, query_length, par, 1, 1, path, path_len, strand);
      break;
    case TMAP_MAP_OPT_SOFT_CLIP_LEFT:
      score = tmap_sw_clipping_core(target, target_length, query, query_length, par, 1, 0, path, path_len, strand);
      break;
    case TMAP_MAP_OPT_SOFT_CLIP_RIGHT:
      score = tmap_sw_clipping_core(target, target_length, query, query_length, par, 0, 1, path, path_len, strand);
      break;
    case TMAP_MAP_OPT_SOFT_CLIP_NONE:
    default:
      score = tmap_sw_clipping_core(target, target_length, query, query_length, par, 0, 0, path, path_len, strand);
      break;
  }
  score_subo = sam->score_subo;

  if(0 < (*path_len) && score_thr < score) {

      if(1 == strand) { // reverse compliment 
          // adjust path
          // TODO

          // reverse the path and adjust its values
          for(i=0;i<(*path_len)>>1;i++) {
              // reverse
              tmap_sw_path_t tmp= path[i];
              path[i] = path[(*path_len)-i-1];
              path[(*path_len)-i-1] = tmp;
              // adjust 
              _tmap_map_util_sw_aux_path_adjust(path[i], query_length, target_length);
              _tmap_map_util_sw_aux_path_adjust(path[(*path_len)-i-1], query_length, target_length);
          }
          if(1 == ((*path_len) & 1)) {
              _tmap_map_util_sw_aux_path_adjust(path[i], query_length, target_length);
          }
      }

      sam->strand = strand;
      sam->seqid = seqid;
      sam->pos = pos + (path[(*path_len)-1].i-1); // zero-based 
      if(0 == strand && path[(*path_len)-1].ctype == TMAP_SW_FROM_I) {
          sam->pos++;
      }
      sam->score = score;
      sam->score_subo = score_subo;
      sam->cigar = tmap_sw_path2cigar(path, (*path_len), &sam->n_cigar);

      if(0 == sam->n_cigar) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }

      // add soft clipping 
      if(1 < path[(*path_len)-1].j) {
          // soft clip the front of the read
          sam->cigar = tmap_realloc(sam->cigar, sizeof(uint32_t)*(1+sam->n_cigar), "sam->cigar");
          for(i=sam->n_cigar-1;0<=i;i--) { // shift up
              sam->cigar[i+1] = sam->cigar[i];
          }
          TMAP_SW_CIGAR_STORE(sam->cigar[0], BAM_CSOFT_CLIP, path[(*path_len)-1].j-1);
          sam->n_cigar++;
      }
      if(path[0].j < query_length) {
          // soft clip the end of the read
          sam->cigar = tmap_realloc(sam->cigar, sizeof(uint32_t)*(1+sam->n_cigar), "sam->cigar");
          TMAP_SW_CIGAR_STORE(sam->cigar[sam->n_cigar], BAM_CSOFT_CLIP, query_length - path[0].j);
          sam->n_cigar++;
      }
  }
  else {
      sam->score = INT32_MIN;
      sam->cigar = NULL;
      sam->n_cigar = 0;
      (*path_len) = 0;
  }
  if(1 == strand) { // reverse back
      tmap_reverse_int(query, query_length);
      tmap_reverse_int(target, target_length);
  }
  return (0 == sam->n_cigar) ? 0 : 1;
}

// TODO: make sure the "longest" read alignment is found
void
tmap_map_util_fsw(tmap_fsw_flowseq_t *fseq, tmap_seq_t *seq, 
                  uint8_t *flow_order, int32_t flow_order_len,
                  uint8_t *key_seq, int32_t key_seq_len,
                  tmap_map_sams_t *sams, tmap_refseq_t *refseq,
                  int32_t bw, int32_t softclip_type, int32_t score_thr,
                  int32_t score_match, int32_t pen_mm, int32_t pen_gapo, 
                  int32_t pen_gape, int32_t fscore)
{
  int32_t i, j, k, l;
  uint8_t *target = NULL;
  int32_t target_mem = 0, target_len = 0;

  tmap_fsw_path_t *path = NULL;
  int32_t path_mem = 0, path_len = 0;
  tmap_fsw_param_t param;
  int32_t matrix[25];
  int32_t start_softclip_len = 0;

  if(0 == sams->n) return;

  // generate the alignment parameters
  param.matrix = matrix;
  param.band_width = 0;
  param.offset = TMAP_MAP_OPT_FSW_OFFSET; // this sets the hp difference
  __tmap_fsw_gen_ap1(param, score_match, pen_mm, pen_gapo, pen_gape, fscore);

  // get flow sequence 
  fseq = tmap_fsw_flowseq_from_seq(fseq, seq, flow_order, flow_order_len, key_seq, key_seq_len);

  // go through each hit
  for(i=0;i<sams->n;i++) {
      tmap_map_sam_t *s = &sams->sams[i];
      uint32_t ref_start, ref_end;

      // get the reference end position
      // NB: soft-clipping at the start may cause ref_start to be moved
      param.band_width = 0;
      start_softclip_len = 0;
      ref_start = ref_end = s->pos + 1;
      for(j=0;j<s->n_cigar;j++) {
          int32_t op, op_len;

          op = TMAP_SW_CIGAR_OP(s->cigar[j]);
          op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[j]);

          switch(op) {
            case BAM_CMATCH:
              ref_end += op_len;
              break;
            case BAM_CDEL:
              if(param.band_width < op_len) param.band_width += op_len;
              ref_end += op_len;
              break;
            case BAM_CINS:
              if(param.band_width < op_len) param.band_width += op_len;
              break;
            case BAM_CSOFT_CLIP:
              if(0 == j) {
                  if(ref_start <= op_len) {
                      start_softclip_len = ref_start;
                      ref_start = 1;
                  }
                  else {
                      start_softclip_len = op_len;
                      ref_start = ref_start - op_len;
                  }
              }
              else ref_end += op_len;
              break;
            default:
              // ignore
              break;
          }
      }
      ref_end--; // NB: since we want this to be one-based

      // check bounds
      if(ref_start < 1) ref_start = 1;
      if(refseq->annos[s->seqid].len < ref_end) {
          ref_end = refseq->annos[s->seqid].len;
      }
      else if(ref_end < 1) ref_end = 1;

      // get the target sequence
      target_len = ref_end - ref_start + 1;
      if(target_mem < target_len) {
          target_mem = target_len;
          tmap_roundup32(target_mem);
          target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }
      target_len = tmap_refseq_subseq(refseq, ref_start + refseq->annos[s->seqid].offset, target_len, target);
      /*
      // NB: IUPAC codes are turned into mismatches
      if(NULL == tmap_refseq_subseq2(refseq, sams->sams[end].seqid+1, start_pos, end_pos, target, 1, NULL)) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }
      */

      if(1 == s->strand) { // reverse compliment
          tmap_reverse_compliment_int(target, target_len);
      }

      // add to the band width
      param.band_width += 2 * bw;

      // make sure we have enough memory for the path
      while(path_mem <= target_len + fseq->num_flows) { // lengthen the path
          path_mem = target_len + fseq->num_flows + 1;
          tmap_roundup32(path_mem);
          path = tmap_realloc(path, sizeof(tmap_fsw_path_t)*path_mem, "path");
      }

      /*
      fprintf(stderr, "strand=%d\n", s->strand);
      fprintf(stderr, "ref_start=%d ref_end=%d\n", ref_start, ref_end);
      fprintf(stderr, "base_calls:\n");
      for(j=0;j<fseq->num_flows;j++) {
          for(k=0;k<fseq->base_calls[j];k++) {
              fputc("ACGTN"[fseq->flow_order[j % fseq->flow_order_len]], stderr);
          }
      }
      fputc('\n', stderr);
      fprintf(stderr, "target:\n");
      for(j=0;j<target_len;j++) {
          fputc("ACGTN"[target[j]], stderr);
      }
      fputc('\n', stderr);
      for(j=0;j<fseq->flow_order_len;j++) {
          fputc("ACGTN"[fseq->flow_order[j]], stderr);
      }
      fputc('\n', stderr);
      */

      // re-align
      s->ascore = s->score;
      path_len = path_mem;
      //fprintf(stderr, "old score=%d\n", s->score);
      switch(softclip_type) {
        case TMAP_MAP_OPT_SOFT_CLIP_ALL:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq, &param, 
                                            1, 1, path, &path_len);
          break;
        case TMAP_MAP_OPT_SOFT_CLIP_LEFT:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq, &param, 
                                            1, 0, path, &path_len);
          break;
        case TMAP_MAP_OPT_SOFT_CLIP_RIGHT:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq, &param, 
                                            0, 1, path, &path_len);
          break;
        case TMAP_MAP_OPT_SOFT_CLIP_NONE:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq, &param, 
                                            0, 0, path, &path_len);
          break;
        default:
          tmap_error("soft clipping type was not recognized", Exit, OutOfRange);
          break;
      }
      s->score_subo = INT32_MIN;
      //fprintf(stderr, "new score=%d path_len=%d\n", s->score, path_len);

      if(0 < path_len) { // update

          /*
          for(j=0;j<path_len;j++) {
              fprintf(stderr, "j=%d path[j].i=%d path[j].j=%d path[j].type=%d\n", j, path[j].i, path[j].j, path[j].ctype);
          }
          */

          // score
          s->score = (int32_t)((s->score + 99.99)/100.0); 

          // position
          s->pos = (ref_start-1);
          // NB: must be careful of leading insertions and strandedness
          if(0 == s->strand) {
              if(0 <= path[path_len-1].j) { 
                  s->pos += (path[path_len-1].j);
              }
          }
          else {
              if(path[0].j < target_len) {
                  s->pos += target_len - path[0].j - 1;
              }
          }
          if(refseq->len < s->pos) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }

          // new cigar
          free(s->cigar);
          s->cigar = tmap_fsw_path2cigar(path, path_len, &s->n_cigar, 1);

          // soft-clipping
          if(0 < path[path_len-1].i) { // skipped beginning flows
              // get the number of bases to clip
              for(j=k=0;j<path[path_len-1].i;j++) {
                  k += fseq->base_calls[j];
              }
              if(0 < k) { // bases should be soft-clipped
                  s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
                  for(l=s->n_cigar-1;0<=l;l--) {
                      s->cigar[l+1] = s->cigar[l];
                  }
                  TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, k);
                  s->n_cigar++;
              }
          }

          // soft-clipping
          if(path[0].i+1 < fseq->num_flows) { // skipped ending flows
              // get the number of bases to clip 
              for(j=path[0].i+1,k=0;j<fseq->num_flows;j++) {
                  k += fseq->base_calls[j];
              }
              if(0 < k) { // bases should be soft-clipped
                  s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
                  s->cigar[s->n_cigar] = (k << 4) | 4;
                  s->n_cigar++;
              }
          }

          // reverse the cigar if on the reverse strand
          if(1 == s->strand) {
              for(i=0;i<s->n_cigar>>1;i++) {
                  uint32_t c;
                  c = s->cigar[i];
                  s->cigar[i] = s->cigar[s->n_cigar-i-1];
                  s->cigar[s->n_cigar-i-1] = c;
              }
          }
      
          s->target_len = 0;
          for(j=0;j<s->n_cigar;j++) {
              switch(TMAP_SW_CIGAR_OP(s->cigar[j])) {
                case BAM_CMATCH:
                case BAM_CDEL:
                  s->target_len += TMAP_SW_CIGAR_LENGTH(s->cigar[j]);
                  break;
                default:
                  break;
              }
          }
      }
  }
  // free
  free(target);
  free(path);
}
