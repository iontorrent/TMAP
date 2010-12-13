#include <stdlib.h>
#include "../util/fmap_alloc.h"
#include "../util/fmap_error.h"
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_sa.h"
#include "../sw/fmap_sw.h"
#include "../sw/fmap_fsw.h"
#include "fmap_map2_aux.h"
#include "fmap_map3.h"
#include "fmap_map_util.h"

static fmap_map_util_fsw_aln_t *
fmap_map_util_fsw_aln_init(int32_t n)
{
  fmap_map_util_fsw_aln_t *aln;
  aln = fmap_calloc(1, sizeof(fmap_map_util_fsw_aln_t), "aln");
  aln->n = n;
  if(0 < n) {
      aln->hits = fmap_calloc(n, sizeof(fmap_map_util_hit_t), "aln->hits");
  }
  return aln;
}

static void
fmap_map_util_fsw_aln_destroy(fmap_map_util_fsw_aln_t *x)
{
  int32_t i;
  for(i=0;i<x->n;i++) {
      free(x->hits[i].cigar);
  }
  free(x->hits);
  free(x);
}

static inline void
fmap_map_util_fsw(fmap_sff_t *sff, fmap_fsw_param_t *par, fmap_map_util_fsw_aln_t *aln, fmap_refseq_t *refseq,
                  int32_t bw, int32_t aln_global, int32_t score_thr)
{
  int32_t i, j, k, l;
  uint8_t *target = NULL;
  int32_t target_mem = 0, target_len = 0;

  fmap_fsw_flowseq_t *fseq[2] = {NULL, NULL};
  fmap_fsw_path_t *path = NULL;
  int32_t path_mem = 0, path_len = 0;
  fmap_fsw_param_t param;
  
  param = (*par);

  // go through each hit
  for(i=0;i<aln->n;i++) {
      fmap_map_util_hit_t *h = &aln->hits[i];
      uint32_t ref_start, ref_end, pacpos;

      // get flow sequence if necessary
      if(NULL == fseq[h->strand]) {
          fseq[h->strand] = fmap_fsw_sff_to_flowseq(sff);
          if(1 == h->strand) fmap_fsw_flowseq_reverse_compliment(fseq[h->strand]);
      }

      param.band_width = 0;
      ref_start = ref_end = h->pos;
      for(j=0;j<h->n_cigar;j++) {
          int32_t op, op_len;

          op = FMAP_SW_CIGAR_OP(h->cigar[j]);
          op_len = FMAP_SW_CIGAR_LENGTH(h->cigar[j]);

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
      if(refseq->annos[h->seqid].len < ref_end) {
          ref_end = refseq->annos[h->seqid].len;
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
            fmap_refseq_seq_i(refseq, pacpos + refseq->annos[h->seqid].offset-1);
      }

      // add to the band width
      param.band_width += 2 * bw;

      // make sure we have enough memory for the path
      while(path_mem <= target_len + fseq[h->strand]->num_flows) { // lengthen the path
          path_mem = target_len + fseq[h->strand]->num_flows;
          fmap_roundup32(path_mem);
          target = fmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }

      // re-align
      if(0 == aln_global) {
          h->score = fmap_fsw_local_core(target, target_len, fseq[h->strand], &param, path, &path_len, score_thr, &h->score_subo);
      }
      else {
          h->score = fmap_fsw_fitting_core(target, target_len, fseq[h->strand], &param, path, &path_len);
          h->score_subo = INT32_MIN;
      }


      if(path_len < 0) { // update
          h->pos = (ref_start-1) + (path[path_len-1].j-1);
          free(h->cigar);
          h->cigar = fmap_fsw_path2cigar(path, path_len, &h->n_cigar, 1);

          if(0 < path[path_len-1].i) { // skipped beginning flows
              // get the number of bases to clip
              for(j=k=0;j<path[path_len-1].i;j++) {
                  k += fseq[h->strand]->base_calls[j];
              }
              if(0 < k) { // bases should be soft-clipped
                  h->cigar = fmap_realloc(h->cigar, sizeof(uint32_t)*(1 + h->n_cigar), "h->cigar");
                  for(l=h->n_cigar-1;0<=l;l--) {
                      h->cigar[l+1] = h->cigar[l];
                  }
                  FMAP_SW_CIGAR_STORE(h->cigar[0], BAM_CSOFT_CLIP, k);
                  h->n_cigar++;
              }
          }

          if(path[0].i < fseq[h->strand]->num_flows) { // skipped ending flows
              // get the number of bases to clip 
              for(j=path[0].i,k=0;j<fseq[h->strand]->num_flows;j++) {
                  k += fseq[h->strand]->base_calls[j];
              }
              if(0 < k) { // bases should be soft-clipped
                  h->cigar = fmap_realloc(h->cigar, sizeof(uint32_t)*(1 + h->n_cigar), "h->cigar");
                  h->cigar[h->n_cigar] = (k << 4) | 4;
                  h->n_cigar++;
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

void
fmap_map_util_map1_adjust_score(fmap_map1_aln_t *aln, int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape)
{
  int32_t i, j, n_match;
  for(i=0;i<aln->n;i++) {
      fmap_map1_hit_t *h = &aln->hits[i];
      n_match = 0 - h->n_mm;
      for(j=0;j<h->n_cigar;j++) {
          switch(FMAP_SW_CIGAR_OP(h->cigar[j])) {
            case BAM_CMATCH:
              n_match += FMAP_SW_CIGAR_LENGTH(h->cigar[j]);
              break;
            default:
              break;
          }
      }

      // update the score
      h->score = n_match * score_match;
      h->score = h->n_mm * pen_mm;
      h->score = h->n_gapo * pen_mm;
      h->score = h->n_gape * pen_mm;
  }
}

/*
inline void
fmap_map_util_map1_fsw(fmap_sff_t *sff, fmap_map1_aln_t *aln, fmap_refseq_t *refseq, fmap_map1_opt_t *opt)
{
  int32_t i;
  fmap_map_util_fsw_aln_t *x;
  fmap_fsw_param_t par;
  int32_t matrix[25];

  x = fmap_map_util_fsw_aln_init(sam->num_entries);

  // generate the alignment parameters
  par.matrix = matrix;
  par.band_width = 0;
  par.offset = FMAP_MAP_UTIL_FSW_OFFSET; // this sets the hp difference
  __fmap_fsw_gen_ap(par, opt);

  // copy aln into x
  // TODO
  for(i=0;i<aln->n;i++) {
      x->hits[i].strand = aln->hits[i].strand;
      x->hits[i].seqid = aln->hits[i].seqid;
      x->hits[i].pos = aln->hits[i].pos;
      x->hits[i].score = aln->hits[i].AS;
      x->hits[i].score_subo = aln->hits[i].XS;
      x->hits[i].n_cigar = aln->hits[i].n_cigar;
      x->hits[i].cigar = aln->hits[i].cigar;
  }

  // re-align
  fmap_map_util_fsw(sff, &par, x, refseq, opt->bw, opt->aln_global, opt->score_thr);

  // copy x into aln
  for(i=0;i<aln->n;i++) {
      aln->hits[i].strand = x->hits[i].strand;
      aln->hits[i].seqid = x->hits[i].seqid;
      aln->hits[i].pos = x->hits[i].pos;
      aln->hits[i].AS = x->hits[i].score;
      aln->hits[i].XS = x->hits[i].score_subo;
      aln->hits[i].n_cigar = x->hits[i].n_cigar;
      aln->hits[i].cigar = x->hits[i].cigar;
      x->hits[i].cigar = NULL; // do not destroy this later
  }

  // destroy
  fmap_map_util_fsw_aln_destroy(x);
}
*/

inline void
fmap_map_util_map2_fsw(fmap_sff_t *sff, fmap_map2_sam_t *sam, fmap_refseq_t *refseq, fmap_map2_opt_t *opt)
{
  int32_t i;
  fmap_map_util_fsw_aln_t *x;
  fmap_fsw_param_t par;
  int32_t matrix[25];

  x = fmap_map_util_fsw_aln_init(sam->num_entries);

  // generate the alignment parameters
  par.matrix = matrix;
  par.band_width = 0;
  par.offset = FMAP_MAP_UTIL_FSW_OFFSET; // this sets the hp difference
  __fmap_fsw_gen_ap(par, opt);

  // copy sam into x
  for(i=0;i<sam->num_entries;i++) {
      x->hits[i].strand = sam->entries[i].strand;
      x->hits[i].seqid = sam->entries[i].seqid;
      x->hits[i].pos = sam->entries[i].pos;
      x->hits[i].score = sam->entries[i].AS;
      x->hits[i].score_subo = sam->entries[i].XS;
      x->hits[i].n_cigar = sam->entries[i].n_cigar;
      x->hits[i].cigar = sam->entries[i].cigar;
  }

  // re-align
  fmap_map_util_fsw(sff, &par, x, refseq, opt->bw, opt->aln_global, opt->score_thr);

  // copy x into sam
  for(i=0;i<sam->num_entries;i++) {
      sam->entries[i].strand = x->hits[i].strand;
      sam->entries[i].seqid = x->hits[i].seqid;
      sam->entries[i].pos = x->hits[i].pos;
      sam->entries[i].AS = x->hits[i].score;
      sam->entries[i].XS = x->hits[i].score_subo;
      sam->entries[i].n_cigar = x->hits[i].n_cigar;
      sam->entries[i].cigar = x->hits[i].cigar;
      x->hits[i].cigar = NULL; // do not destroy this later
  }

  // destroy
  fmap_map_util_fsw_aln_destroy(x);
}

inline void
fmap_map_util_map3_fsw(fmap_sff_t *sff, fmap_map3_aln_t *aln, fmap_refseq_t *refseq, fmap_map3_opt_t *opt)
{
  int32_t i;
  fmap_map_util_fsw_aln_t *x;
  fmap_fsw_param_t par;
  int32_t matrix[25];
  
  x = fmap_map_util_fsw_aln_init(aln->n);

  // generate the alignment parameters
  par.matrix = matrix;
  par.band_width = 0;
  par.offset = FMAP_MAP_UTIL_FSW_OFFSET; // this sets the hp difference
  __fmap_fsw_gen_ap(par, opt);

  // copy aln into x
  for(i=0;i<aln->n;i++) {
      x->hits[i].strand = aln->hits[i].strand;
      x->hits[i].seqid = aln->hits[i].seqid;
      x->hits[i].pos = aln->hits[i].pos;
      x->hits[i].score = aln->hits[i].score;
      x->hits[i].score_subo = aln->hits[i].score_subo;
      x->hits[i].n_cigar = aln->hits[i].n_cigar;
      x->hits[i].cigar = aln->hits[i].cigar;
  }

  // re-align
  fmap_map_util_fsw(sff, &par, x, refseq, opt->bw, opt->aln_global, opt->score_thr);

  // copy x into aln
  for(i=0;i<aln->n;i++) {
      aln->hits[i].strand = x->hits[i].strand;
      aln->hits[i].seqid = x->hits[i].seqid;
      aln->hits[i].pos = x->hits[i].pos;
      aln->hits[i].score = x->hits[i].score;
      aln->hits[i].score_subo = x->hits[i].score_subo;
      aln->hits[i].n_cigar = x->hits[i].n_cigar;
      aln->hits[i].cigar = x->hits[i].cigar;
      x->hits[i].cigar = NULL; // do not destroy this later
  }

  // destroy
  fmap_map_util_fsw_aln_destroy(x);
}
