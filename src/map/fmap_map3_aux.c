#include <stdlib.h>
#include "../util/fmap_alloc.h"
#include "../util/fmap_sort.h"
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_sa.h"
#include "../index/fmap_bwt_match.h"
#include "../sw/fmap_sw.h"
#include "fmap_map3.h"
#include "fmap_map3_aux.h"

#define __fmap_map3_hit_sort_lt(a, b) ( ((a).seqid < (b).seqid \
                                         || ( (a).seqid == (b).seqid && (a).pos < (b).pos ) \
                                         || ( (a).seqid == (b).seqid && (a).pos == (b).pos && (a).start < (b).start )) \
                                       ? 1 : 0)

FMAP_SORT_INIT(fmap_map3_aux_hit_t, fmap_map3_aux_hit_t, __fmap_map3_hit_sort_lt)

     // Note: this does not set the band width
#define __map3_gen_ap(par, opt) do { \
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
} while(0)

fmap_map3_aln_t *
fmap_map3_aln_init()
{
  return fmap_calloc(1, sizeof(fmap_map3_aln_t), "aln");
}

void
fmap_map3_aln_realloc(fmap_map3_aln_t *aln, int32_t n)
{
  int32_t i;

  for(i=n;i<aln->n;i++) {
      free(aln->hits[i].cigar);
  }
  aln->hits = fmap_realloc(aln->hits, sizeof(fmap_map3_hit_t)*n, "aln->hits");
  for(i=aln->n;i<n;i++) { // overly paranoid?
      aln->hits[i].strand = 0;
      aln->hits[i].seqid = aln->hits[i].pos = 0;
      aln->hits[i].score = 0;
      aln->hits[i].n_seeds = 0;
      aln->hits[i].n_cigar = 0;
      aln->hits[i].cigar = NULL;
  }
  aln->n = n;
}

void
fmap_map3_aln_destroy(fmap_map3_aln_t *aln)
{
  int32_t i;
  for(i=0;i<aln->n;i++) {
      free(aln->hits[i].cigar);
  }
  free(aln->hits);
  free(aln);
}

static inline void 
fmap_map3_aux_seed_add(fmap_map3_aux_seed_t **seeds,
                       int32_t *n_seeds,
                       int32_t *m_seeds,
                       int32_t k,
                       int32_t l,
                       int32_t start, 
                       int32_t offset)
{
  if((*m_seeds) <= (*n_seeds)) {
      (*m_seeds) = (0 == (*m_seeds)) ? 64 : ((*m_seeds) << 1);
      (*seeds) = fmap_realloc((*seeds), sizeof(fmap_map3_aux_seed_t)*(*m_seeds), "(*seeds)");
  }
  (*seeds)[(*n_seeds)].k = k;
  (*seeds)[(*n_seeds)].l = l;
  (*seeds)[(*n_seeds)].start = start;
  (*seeds)[(*n_seeds)].offset = offset;
  (*n_seeds)++;
}

static inline void
fmap_map3_aux_core_seed_helper(uint8_t *query,
                               int32_t query_length,
                               int32_t offset,
                               uint8_t *flow_order,
                               uint8_t flow_i,
                               fmap_refseq_t *refseq,
                               fmap_bwt_t *bwt,
                               fmap_sa_t *sa,
                               fmap_map3_opt_t *opt,
                               fmap_map3_aux_seed_t **seeds,
                               int32_t *n_seeds,
                               int32_t *m_seeds)
{ 
  int32_t i, k;
  int32_t n_bases;;
  fmap_bwt_match_occ_t prev_sa, cur_sa, next_sa, tmp_sa;

  // DEBUGGING
  if(query_length <= offset) return;
  if(flow_order[flow_i] != query[offset]) {
      fmap_error("bug encountered", Exit, OutOfRange);
  }

  // initialize prev
  prev_sa.k = 0;
  prev_sa.l = bwt->seq_len;
  prev_sa.hi = 0;
  prev_sa.offset = 0;

  i = offset;
  while(i < query_length) {
      // reached the seed length
      if(opt->seed_length < i - offset) { 
          break;
      }

      // get the homopolymer length
      n_bases = 0;
      if(flow_order[flow_i] == query[i]) { // non-empty flow
          n_bases = 1;
          while(i + n_bases < query_length && query[i] == query[i+n_bases]) {
              n_bases++;
          }
      }


      // move through the homopolymer, trying deletions if possible
      next_sa = prev_sa;
      for(k=0;k<n_bases;k++) {
          // reached the seed length
          if(opt->seed_length < i - offset + k) { 
              break;
          }

          // only delete if there are bases available and we are not deleting
          // the entire first flow 
          int32_t bases_to_align = opt->seed_length - (i - offset + k);
          int32_t bases_left = query_length - i - n_bases;
          if(0 < n_bases - k // bases to delete
             && n_bases - k <= opt->hp_diff // not too many to delete
             && (i != offset || 0 != k) // do not delete the entire flow
             && bases_to_align <= bases_left) { // enough bases 
              // match exactly from here onwards
              tmp_sa = next_sa;
              if(0 < fmap_bwt_match_exact_alt(bwt, bases_to_align, query + i + n_bases, &tmp_sa)
                 && (tmp_sa.l - tmp_sa.k + 1) <= opt->max_seed_hits) {
                  fmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, tmp_sa.k, tmp_sa.l, offset, n_bases-k);
              }
          }

          // move past this base in the hp
          tmp_sa = next_sa;
          fmap_bwt_match_2occ(bwt, &tmp_sa, flow_order[flow_i], &next_sa);
          if(next_sa.l < next_sa.k) { // no match, return
              return;
          }
      }

      // insert hp bases
      if(i + n_bases < offset + opt->seed_length) { // not the last flow
          cur_sa = next_sa; // already considered the 'n_bases' of this flow
          // insert
          for(k=1;k<=opt->hp_diff;k++) { // # of bases to insert
              fmap_bwt_match_2occ(bwt, &cur_sa, flow_order[flow_i], &tmp_sa);
              if(tmp_sa.l < tmp_sa.k) { // no match, do not continue
                  break;
              }
              // match exactly from here onwards
              if(0 < fmap_bwt_match_exact_alt(bwt, opt->seed_length - (i - offset) - n_bases, query + i + n_bases, &tmp_sa)
                 && (tmp_sa.l - tmp_sa.k + 1) <= opt->max_seed_hits) {
                  fmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, tmp_sa.k, tmp_sa.l, offset, k);
              }
              // move to the next
              cur_sa = tmp_sa;
          }
      }

      // next flow
      flow_i = (1+flow_i) & 3;
      i += n_bases;
      prev_sa = next_sa;
  }

  if(i - offset < opt->seed_length) fmap_error("bug encountered", Exit, OutOfRange);

  // add in the seed with no hp indels
  if((next_sa.l - next_sa.k + 1) <= opt->max_seed_hits) {
      fmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, next_sa.k, next_sa.l, offset, 0);
  }
}

static inline void
fmap_map3_aux_core_seed(uint8_t *query,
                        int32_t query_length,
                        uint8_t *flow_order,
                        fmap_refseq_t *refseq,
                        fmap_bwt_t *bwt,
                        fmap_sa_t *sa,
                        fmap_map3_opt_t *opt,
                        fmap_map3_aux_seed_t **seeds,
                        int32_t *n_seeds,
                        int32_t *m_seeds)
{
  int32_t i, j, flow_i;

  if(0 < opt->hp_diff) {
      i=flow_i=0;
      while(i<query_length - opt->seed_length + 1) {
          // move to the next flow
          j=0;
          while(query[i] != flow_order[flow_i]) {
              flow_i = (flow_i + 1) & 3;
              // sanity check
              j++;
              if(4 <= j) fmap_error("bug encountered", Exit, OutOfRange);
          }

          // add seeds
          fmap_map3_aux_core_seed_helper(query, query_length, i, flow_order, flow_i,
                                         refseq, bwt, sa, opt, seeds, n_seeds, m_seeds);

          // skip over this hp
          i++;
          while(i<query_length) {
              if(query[i] != query[i-1]) {
                  break;
              }
              i++;
          }
      }
  }
  else {
      fmap_bwt_match_occ_t cur_sa;
      for(i=0;i<query_length-opt->seed_length+1;i++) {
          if(0 < fmap_bwt_match_exact(bwt, opt->seed_length, query + i, &cur_sa)
                 && (cur_sa.l - cur_sa.k + 1) <= opt->max_seed_hits) {
              fmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, cur_sa.k, cur_sa.l, i, 0);
          }
      }
  }
}

// TODO: memory pools?
fmap_map3_aln_t *
fmap_map3_aux_core(fmap_seq_t *seq[2], 
                   uint8_t *flow[2],
                   fmap_refseq_t *refseq,
                   fmap_bwt_t *bwt,
                   fmap_sa_t *sa,
                   fmap_map3_opt_t *opt)
{
  int32_t i, j;

  int32_t seq_len[2];
  fmap_string_t *bases;
  uint8_t *query;

  uint32_t k, pacpos;
  uint32_t ref_start, ref_end;
  uint8_t *target = NULL;
  int32_t target_mem = 0, target_len;
  int32_t matrix[25];
  fmap_sw_param_t par;
  fmap_sw_path_t *path = NULL;
  int32_t path_len, path_mem=0, score, score_subo;
  int32_t bw = 0;

  fmap_map3_aux_seed_t *seeds[2];
  int32_t m_seeds[2], n_seeds[2];

  fmap_map3_aux_hit_t *hits[2];
  int32_t m_hits[2], n_hits[2];

  fmap_map3_aln_t *aln = NULL;

  // scoring matrix
  par.matrix = matrix;
  __map3_gen_ap(par, opt);

  aln = fmap_map3_aln_init();

  // band width
  bw = (opt->bw + 1) / 2;

  // seeds
  for(i=0;i<2;i++) { // forward/reverse-compliment
      bases = fmap_seq_get_bases(seq[i]);
      seq_len[i] = bases->l;
      query = (uint8_t*)bases->s;

      // pre-allocate mmemory
      n_seeds[i] = 0;
      m_seeds[i] = seq_len[i] - opt->seed_length + 1; // maximum number of seeds possible
      seeds[i] = fmap_malloc(m_seeds[i]*sizeof(fmap_map3_aux_seed_t), "seeds[i]");

      fmap_map3_aux_core_seed(query, seq_len[i], flow[i],
                              refseq, bwt, sa, opt, &seeds[i], &n_seeds[i], &m_seeds[i]);
  }

  // convert seeds to chr/pos
  for(i=0;i<2;i++) { // forward/reverse-compliment
      // allocate mmemory
      m_hits[i] = 0; // store the number of hits that are needed
      for(j=0;j<n_seeds[i];j++) {
          m_hits[i] += seeds[i][j].l - seeds[i][j].k + 1;
      }
      hits[i] = fmap_malloc(m_hits[i]*sizeof(fmap_map3_aux_hit_t), "hits[i]");
      n_hits[i] = 0;

      for(j=0;j<n_seeds[i];j++) { // go through all seeds
          for(k=seeds[i][j].k;k<=seeds[i][j].l;k++) { // through all occurrences
              pacpos = fmap_sa_pac_pos(sa, bwt, k);
              if(bwt->seq_len < pacpos + seeds[i][j].start + 1) { // before the beginning of the reference sequence
                  pacpos = 0;
              }
              else {
                  pacpos = bwt->seq_len - pacpos - seeds[i][j].start - 1;
              }
              if(0 < fmap_refseq_pac2real(refseq, pacpos, 1,
                                          &hits[i][n_hits[i]].seqid, &hits[i][n_hits[i]].pos)) {
                  if(hits[i][n_hits[i]].pos < opt->seed_length + seeds[i][j].offset - 1 ) {
                      hits[i][n_hits[i]].pos = 0;
                  }
                  else {
                      hits[i][n_hits[i]].pos -= opt->seed_length + seeds[i][j].offset - 1;
                  }
                  hits[i][n_hits[i]].start = seeds[i][j].start;
                  n_hits[i]++;
              }
          }
      }
  }

  // free the seeds
  for(i=0;i<2;i++) {
      free(seeds[i]);
      seeds[i]=NULL;
  }

  // sort hits
  for(i=0;i<2;i++) {
      fmap_sort_introsort(fmap_map3_aux_hit_t, n_hits[i], hits[i]);
  } 

  // band hits
  for(i=0;i<2;i++) {
      bases = fmap_seq_get_bases(seq[i]);
      seq_len[i] = bases->l;
      query = (uint8_t*)bases->s;
  
      int32_t start, end;
      start = end = 0;
      while(end < n_hits[i]) {
          if(end+1 < n_hits[i]) { 
              // check if the next interval can be banded
              if(hits[i][end].seqid == hits[i][end+1].seqid
                 && hits[i][end+1].pos - hits[i][end].pos <= opt->max_seed_band) {
                  end++;
                  continue; // there may be more to add
              }
          }

          // get the start of the target range
          if(hits[i][start].pos < bw) {
              ref_start = 1;
          }
          else {
              ref_start = hits[i][start].pos - bw + 1;
              // check bounds
              if(ref_start < 1) {
                  ref_start = 1;
              }
          }
          // get the end of the target range
          ref_end = hits[i][end].pos + seq_len[i] + bw - 1;
          if(refseq->annos[hits[i][end].seqid].len < ref_end) {
              // this assumes that the seed matched correctly (do not run
              // off the end)
              ref_end = refseq->annos[hits[i][end].seqid].len;
          }

          // get the target sequence
          target_len = ref_end - ref_start + 1;
          if(target_mem < target_len) { // more memory?
              target_mem = target_len;
              fmap_roundup32(target_mem);
              target = fmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
          }
          for(pacpos=ref_start;pacpos<=ref_end;pacpos++) {

              // add contig offset and make zero based
              target[pacpos-ref_start] = fmap_refseq_seq_i(refseq, pacpos + refseq->annos[hits[i][end].seqid].offset-1);
          }

          // get the band width
          par.band_width = hits[i][end].pos - hits[i][start].pos;
          par.band_width += 2 * bw; // add bases to the window

          if(path_mem <= target_len + seq_len[i]) { // lengthen the path
              path_mem = target_len + seq_len[i];
              fmap_roundup32(path_mem);
              path = fmap_realloc(path, sizeof(fmap_sw_path_t)*path_mem, "path");
          }

          // threshold the score by assuming that one seed's worth of
          // matches occurs in the alignment
          if(0 == opt->aln_global) {
              score = fmap_sw_local_core(target, target_len, query, seq_len[i], &par, path, &path_len, opt->score_thr, &score_subo);
          }
          else {
              score = fmap_sw_fitting_core(target, target_len, query, seq_len[i], &par, path, &path_len);
              score_subo = INT32_MIN;
          }

          if(0 < path_len) {
              fmap_map3_hit_t *hit;

              aln->n++;
              aln->hits = fmap_realloc(aln->hits, aln->n*sizeof(fmap_map3_hit_t), "aln->hits");

              // save the hit
              hit = &aln->hits[aln->n-1]; // for easy of writing code
              hit->strand = i;
              hit->seqid = hits[i][start].seqid; 
              hit->pos = (ref_start-1) + (path[path_len-1].i-1); // zero-based 
              hit->score = score;
              hit->score_subo = score_subo;
              hit->n_seeds = ((1 << 15) < end - start + 1) ? (1 << 15) : (end - start + 1);
              hit->cigar = fmap_sw_path2cigar(path, path_len, &hit->n_cigar);

              // add soft clipping after local alignment
              if(1 < path[path_len-1].j) {
                  // soft clip the front of the read
                  hit->cigar = fmap_realloc(hit->cigar, sizeof(uint32_t)*(1+hit->n_cigar), "hit->cigar");
                  for(j=hit->n_cigar-1;0<=j;j--) { // shift up
                      hit->cigar[j+1] = hit->cigar[j];
                  }
                  hit->cigar[0] = ((path[path_len-1].j-1) << 4) | 4; 
                  hit->n_cigar++;
              }
              if(path[0].j < seq_len[i]) { // 
                  // soft clip the end of the read
                  hit->cigar = fmap_realloc(hit->cigar, sizeof(uint32_t)*(1+hit->n_cigar), "hit->cigar");
                  hit->cigar[hit->n_cigar] = ((seq_len[i] - path[0].j) << 4) | 4; 
                  hit->n_cigar++;
              }
          }

          // update start/end
          end++;
          start = end;
      }
  }

  // free the hits
  for(i=0;i<2;i++) {
      free(hits[i]);
  }
  free(target);

  // free
  free(path);

  return aln;
}
