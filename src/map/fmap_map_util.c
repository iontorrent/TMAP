#include <stdlib.h>
#include "../util/fmap_alloc.h"
#include "../util/fmap_error.h"
#include "../util/fmap_sam.h"
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_sa.h"
#include "../sw/fmap_sw.h"
#include "../sw/fmap_fsw.h"
#include "fmap_map_util.h"

void
fmap_map_sam_malloc_aux(fmap_map_sam_t *s, int32_t algo_id)
{
  switch(s->algo_id) {
    case FMAP_MAP_ALGO_MAP1:
      s->aux.map1_aux = fmap_calloc(1, sizeof(fmap_map_map1_aux_t), "s->aux.map1_aux");
      break;
    case FMAP_MAP_ALGO_MAP2:
      s->aux.map2_aux = fmap_calloc(1, sizeof(fmap_map_map2_aux_t), "s->aux.map2_aux");
      break;
    case FMAP_MAP_ALGO_MAP3:
      s->aux.map3_aux = fmap_calloc(1, sizeof(fmap_map_map3_aux_t), "s->aux.map3_aux");
      break;
    default:
      break;
  }
}

inline void
fmap_map_sam_destroy_aux(fmap_map_sam_t *s)
{
  switch(s->algo_id) {
    case FMAP_MAP_ALGO_MAP1:
      free(s->aux.map1_aux);
      s->aux.map1_aux = NULL;
      break;
    case FMAP_MAP_ALGO_MAP2:
      free(s->aux.map2_aux);
      s->aux.map2_aux = NULL;
      break;
    case FMAP_MAP_ALGO_MAP3:
      free(s->aux.map3_aux);
      s->aux.map3_aux = NULL;
      break;
    default:
      break;
  }
}

void
fmap_map_sam_destroy(fmap_map_sam_t *s)
{
  fmap_map_sam_destroy_aux(s);
  free(s->cigar);
  s->cigar = NULL;
  s->n_cigar = 0;
}

fmap_map_sams_t *
fmap_map_sams_init()
{
  return fmap_calloc(1, sizeof(fmap_map_sams_t), "return");
}

void
fmap_map_sams_realloc(fmap_map_sams_t *s, int32_t n)
{
  int32_t i;
  if(n == s->n) return; 
  for(i=n;i<s->n;i++) {
      fmap_map_sam_destroy_aux(&s->sams[i]);
  }
  s->sams = fmap_realloc(s->sams, sizeof(fmap_map_sam_t) * n, "s->sams");
  s->n = n;
}

void
fmap_map_sams_destroy(fmap_map_sams_t *s)
{
  int32_t i;
  if(NULL == s) return;
  for(i=0;i<s->n;i++) {
      fmap_map_sam_destroy_aux(&s->sams[i]);
  }
  free(s->sams);
  free(s);
}

static void
fmap_map_sam_print(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_map_sam_t *sam, int32_t sam_sff_tags)
{
  if(NULL == sam) { // unmapped
      fmap_sam_print_unmapped(fmap_file_stdout, seq, sam_sff_tags);
  }
  else {
      switch(sam->algo_id) {
        case FMAP_MAP_ALGO_MAP1:
          fmap_sam_print_mapped(fmap_file_stdout, seq, sam_sff_tags, refseq, 
                                sam->strand, sam->seqid, sam->pos,
                                sam->mapq, sam->cigar, sam->n_cigar,
                                sam->score, sam->algo_id, sam->algo_stage, 
                                "\tXM:i:%d\tXO:i:%d\tXG:i:%d",
                                sam->aux.map1_aux->n_mm,
                                sam->aux.map1_aux->n_gapo,
                                sam->aux.map1_aux->n_gape);
          break;
        case FMAP_MAP_ALGO_MAP2:
          if(0 < sam->aux.map2_aux->XI) {
              fmap_sam_print_mapped(fmap_file_stdout, seq, sam_sff_tags, refseq, 
                                    sam->strand, sam->seqid, sam->pos,
                                    sam->mapq, sam->cigar, sam->n_cigar,
                                    sam->score, sam->algo_id, sam->algo_stage, 
                                    "\tXS:i:%d\tXF:i:%d\tXE:i:%d\tXI:i:%d",
                                    sam->score_subo,
                                    sam->aux.map2_aux->XF, sam->aux.map2_aux->XE, 
                                    sam->aux.map2_aux->XI);
          }
          else {
              fmap_sam_print_mapped(fmap_file_stdout, seq, sam_sff_tags, refseq, 
                                    sam->strand, sam->seqid, sam->pos,
                                    sam->mapq, sam->cigar, sam->n_cigar,
                                    sam->score, sam->algo_id, sam->algo_stage, 
                                    "\tXS:i:%d\tXF:i:%d\tXE:i:%d",
                                    sam->score_subo,
                                    sam->aux.map2_aux->XF, sam->aux.map2_aux->XE);
          }
          break;
        case FMAP_MAP_ALGO_MAP3:
          fmap_sam_print_mapped(fmap_file_stdout, seq, sam_sff_tags, refseq, 
                                sam->strand, sam->seqid, sam->pos,
                                sam->mapq, sam->cigar, sam->n_cigar,
                                sam->score, sam->algo_id, sam->algo_stage, 
                                "\tXE:i:%d",
                                sam->aux.map3_aux->n_seeds);
          break;
      }
  }
}

void 
fmap_map_sams_print(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_map_sams_t *sams, int32_t sam_sff_tags) 
{
  int32_t i;
  if(0 < sams->n) {
      for(i=0;i<sams->n;i++) {
          fmap_map_sam_print(seq, refseq, &sams->sams[i], sam_sff_tags);
      }
  }
  else {
      fmap_map_sam_print(seq, refseq, NULL, sam_sff_tags);
  }
}

void
fmap_map_sam_copy_and_nullify(fmap_map_sam_t *dest, fmap_map_sam_t *src)
{
  (*dest) = (*src);
  src->cigar = NULL;
  switch(src->algo_id) {
    case FMAP_MAP_ALGO_MAP1:
      src->aux.map1_aux = NULL;
      break;
    case FMAP_MAP_ALGO_MAP2:
      src->aux.map2_aux = NULL;
      break;
    case FMAP_MAP_ALGO_MAP3:
      src->aux.map3_aux = NULL;
      break;
    default:
      break;
  }
}

void
fmap_map_sams_filter1(fmap_map_sams_t *sams, int32_t aln_output_mode, int32_t algo_id)
{
  int32_t i, j, k;
  int32_t n_best = 0;
  int32_t best_score, cur_score;

  if(sams->n <= 1) {
      return;
  }

  best_score = INT32_MIN;
  n_best = 0;
  for(i=0;i<sams->n;i++) {
      if(FMAP_MAP_ALGO_NONE == algo_id
         || sams->sams[i].algo_id == algo_id) {
          cur_score = sams->sams[i].score;
          if(best_score < cur_score) {
              best_score = cur_score;
              n_best = 1;
          }
          else if(!(cur_score < best_score)) { // equal
              n_best++;
          }
      }
  }

  // adjust mapping qualities
  if(1 < n_best) {
      for(i=0;i<sams->n;i++) {
          if(FMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              sams->sams[i].mapq = 0;
          }
      }
  }
  else {
      for(i=0;i<sams->n;i++) {
          if(FMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              cur_score = sams->sams[i].score;
              if(cur_score < best_score) { // not the best
                  sams->sams[i].mapq = 0;
              }
          }
      }
  }

  if(FMAP_MAP_UTIL_ALN_MODE_ALL == aln_output_mode) {
      return;
  }

  // copy to the front
  if(n_best < sams->n) {
      for(i=j=0;i<sams->n;i++) {
          if(FMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              cur_score = sams->sams[i].score;
              if(cur_score < best_score) { // not the best
                  fmap_map_sam_destroy(&sams->sams[i]);
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                  }
                  j++;
              }
          }
          else {
              if(j < i) { // copy if we are not on the same index
                  fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
              }
              j++;
          }
      }
      // reallocate
      fmap_map_sams_realloc(sams, j);
  }

  if(FMAP_MAP_UTIL_ALN_MODE_UNIQ_BEST == aln_output_mode) {
      if(1 < n_best) { // there can only be one
          if(FMAP_MAP_ALGO_NONE == algo_id) {
              fmap_map_sams_realloc(sams, 0);
          }
          else {
              // get rid of all of them
              for(i=j=0;i<sams->n;i++) {
                  if(sams->sams[i].algo_id == algo_id) {
                      fmap_map_sam_destroy(&sams->sams[i]);
                  }
                  else {
                      if(j < i) { // copy if we are not on the same index
                          fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                      }
                      j++;
                  }
              }
              fmap_map_sams_realloc(sams, j);
          }
      }
  }
  else if(FMAP_MAP_UTIL_ALN_MODE_RAND_BEST == aln_output_mode) { // get a random
      int32_t r = (int32_t)(drand48() * n_best);

      // keep the rth one
      if(FMAP_MAP_ALGO_NONE == algo_id) {
          if(0 != r) {
              fmap_map_sam_copy_and_nullify(&sams->sams[0], &sams->sams[r]);
          }
          // reallocate
          fmap_map_sams_realloc(sams, 1);
      }
      else {
          // keep the rth one
          for(i=j=k=0;i<sams->n;i++) {
              if(sams->sams[i].algo_id == algo_id) {
                  if(k == r) { // keep
                      if(j < i) { // copy if we are not on the same index
                          fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                      }
                      j++;
                  }
                  else { // free
                      fmap_map_sam_destroy(&sams->sams[i]);
                  }
                  k++;
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      fmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                  }
                  j++;
              }
          }
          fmap_map_sams_realloc(sams, j);
      }
  }
  else if(FMAP_MAP_UTIL_ALN_MODE_ALL_BEST == aln_output_mode) {
      // do nothing
  }
  else {
      fmap_error("bug encountered", Exit, OutOfRange);
  }
}

void
fmap_map_sams_filter(fmap_map_sams_t *sams, int32_t aln_output_mode)
{
  fmap_map_sams_filter1(sams, aln_output_mode, FMAP_MAP_ALGO_NONE);
}

void
fmap_map_util_map1_adjust_score(fmap_map_sams_t *sams, int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape)
{
  int32_t i, j, n_match;
  for(i=0;i<sams->n;i++) {
      fmap_map_sam_t *sam = &sams->sams[i];
      n_match = 0 - sam->aux.map1_aux->n_mm;
      for(j=0;j<sam->n_cigar;j++) {
          switch(FMAP_SW_CIGAR_OP(sam->cigar[j])) {
            case BAM_CMATCH:
              n_match += FMAP_SW_CIGAR_LENGTH(sam->cigar[j]);
              break;
            default:
              break;
          }
      }

      // update the score
      sam->score = n_match * score_match;
      sam->score -= sam->aux.map1_aux->n_mm * pen_mm;
      sam->score -= sam->aux.map1_aux->n_gapo * pen_gapo;
      sam->score -= sam->aux.map1_aux->n_gape * pen_gape;
  }
}

void
fmap_map_util_fsw(fmap_sff_t *sff, 
                  fmap_map_sams_t *sams, fmap_refseq_t *refseq,
                  int32_t bw, int32_t aln_global, int32_t score_thr,
                  int32_t score_match, int32_t pen_mm, int32_t pen_gapo, 
                  int32_t pen_gape, int32_t fscore)
{
  int32_t i, j, k, l;
  uint8_t *target = NULL;
  int32_t target_mem = 0, target_len = 0;

  fmap_fsw_flowseq_t *fseq[2] = {NULL, NULL};
  fmap_fsw_path_t *path = NULL;
  int32_t path_mem = 0, path_len = 0;
  fmap_fsw_param_t param;
  int32_t matrix[25];

  // generate the alignment parameters
  param.matrix = matrix;
  param.band_width = 0;
  param.offset = FMAP_MAP_UTIL_FSW_OFFSET; // this sets the hp difference
  __fmap_fsw_gen_ap1(param, score_match, pen_mm, pen_gapo, pen_gape, fscore);

  // go through each hit
  for(i=0;i<sams->n;i++) {
      fmap_map_sam_t *s = &sams->sams[i];
      uint32_t ref_start, ref_end, pacpos;

      // get flow sequence if necessary
      if(NULL == fseq[s->strand]) {
          fseq[s->strand] = fmap_fsw_sff_to_flowseq(sff);
          if(1 == s->strand) fmap_fsw_flowseq_reverse_compliment(fseq[s->strand]);
      }

      param.band_width = 0;
      ref_start = ref_end = s->pos;
      for(j=0;j<s->n_cigar;j++) {
          int32_t op, op_len;

          op = FMAP_SW_CIGAR_OP(s->cigar[j]);
          op_len = FMAP_SW_CIGAR_LENGTH(s->cigar[j]);

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
                      ref_start = 1;
                  }
                  else {
                      ref_start = ref_start - op_len + 1;
                  }
              }
              else ref_end += op_len;
              break;
            default:
              // ignore
              break;
          }
      }

      // check bounds
      if(ref_start < 1) ref_start = 1;
      if(refseq->annos[s->seqid].len < ref_end) {
          ref_end = refseq->annos[s->seqid].len;
      }

      // get the target sequence
      target_len = ref_end - ref_start + 1;
      if(target_mem < target_len) {
          target_mem = target_len;
          fmap_roundup32(target_mem);
          target = fmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }
      for(pacpos=ref_start;pacpos<=ref_end;pacpos++) {
          // add contig offset and make zero-based
          target[pacpos-ref_start] =
            fmap_refseq_seq_i(refseq, pacpos + refseq->annos[s->seqid].offset-1);
      }

      // add to the band width
      param.band_width += 2 * bw;

      // make sure we have enough memory for the path
      while(path_mem <= target_len + fseq[s->strand]->num_flows) { // lengthen the path
          path_mem = target_len + fseq[s->strand]->num_flows;
          fmap_roundup32(path_mem);
          target = fmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }

      // re-align
      if(0 == aln_global) {
          s->score = fmap_fsw_local_core(target, target_len, fseq[s->strand], &param, path, &path_len, score_thr, &s->score_subo);
      }
      else {
          s->score = fmap_fsw_fitting_core(target, target_len, fseq[s->strand], &param, path, &path_len);
          s->score_subo = INT32_MIN;
      }


      if(path_len < 0) { // update
          s->pos = (ref_start-1) + (path[path_len-1].j-1);
          free(s->cigar);
          s->cigar = fmap_fsw_path2cigar(path, path_len, &s->n_cigar, 1);

          if(0 < path[path_len-1].i) { // skipped beginning flows
              // get the number of bases to clip
              for(j=k=0;j<path[path_len-1].i;j++) {
                  k += fseq[s->strand]->base_calls[j];
              }
              if(0 < k) { // bases should be soft-clipped
                  s->cigar = fmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
                  for(l=s->n_cigar-1;0<=l;l--) {
                      s->cigar[l+1] = s->cigar[l];
                  }
                  FMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, k);
                  s->n_cigar++;
              }
          }

          if(path[0].i < fseq[s->strand]->num_flows) { // skipped ending flows
              // get the number of bases to clip 
              for(j=path[0].i,k=0;j<fseq[s->strand]->num_flows;j++) {
                  k += fseq[s->strand]->base_calls[j];
              }
              if(0 < k) { // bases should be soft-clipped
                  s->cigar = fmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
                  s->cigar[s->n_cigar] = (k << 4) | 4;
                  s->n_cigar++;
              }
          }
      }
  }
  // free
  if(NULL != fseq[0]) fmap_fsw_flowseq_destroy(fseq[0]);
  if(NULL != fseq[1]) fmap_fsw_flowseq_destroy(fseq[1]);
  free(target);
  free(path);
}
