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
                                         || ( (a).seqid == (b).seqid && (a).pos == (b).pos && (a).offset < (b).offset )) \
                                       ? 1 : 0)

FMAP_SORT_INIT(fmap_map3_aux_hit_t, fmap_map3_aux_hit_t, __fmap_map3_hit_sort_lt)

     // Note: this does not set the band width
#define __map3_gen_ap(par, op) do { \
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

// TODO: memory pools?
fmap_map3_aln_t *
fmap_map3_aux_core(fmap_seq_t *seq[2], 
                   fmap_refseq_t *refseq,
                   fmap_bwt_t *bwt,
                   fmap_sa_t *sa,
                   fmap_map3_opt_t *opt)
{
  int32_t i, j;
  int32_t seq_len[2];
  fmap_string_t *bases;
  fmap_bwt_match_occ_t match_sa;
  uint8_t *query;
  uint32_t k, pacpos;
  uint32_t ref_start, ref_end;
  uint8_t *target = NULL;
  int32_t target_mem = 0, target_len;
  int32_t matrix[25];
  fmap_sw_param_t par;
  fmap_sw_path_t *path = NULL;
  int32_t path_len, score, score_subo;

  fmap_map3_aux_seed_t *seeds[2];
  int32_t m_seeds[2], n_seeds[2];

  fmap_map3_aux_hit_t *hits[2];
  int32_t m_hits[2], n_hits[2];

  fmap_map3_aln_t *aln = NULL;

  // scoring matrix
  par.matrix = matrix;
  __map3_gen_ap(par, opt);

  aln = fmap_map3_aln_init();

  // seeds
  for(i=0;i<2;i++) { // forward/reverse-compliment
      bases = fmap_seq_get_bases(seq[i]);
      seq_len[i] = bases->l;
      query = (uint8_t*)bases->s;

      // allocate mmemory
      n_seeds[i] = 0;
      m_seeds[i] = seq_len[i] - opt->seed_length + 1; // maximum number of seeds possible
      seeds[i] = fmap_malloc(m_seeds[i]*sizeof(fmap_map3_aux_seed_t), "seeds[i]");

      m_hits[i] = 0; // store the number of hits that are needed

      for(j=0;j<seq_len[i] - opt->seed_length + 1;j++) { // all valid offsets
          if(0 < fmap_bwt_match_exact(bwt, opt->seed_length, query + j, &match_sa)
             && (match_sa.l - match_sa.k + 1) <= opt->max_seed_hits) {
              seeds[i][n_seeds[i]].k = match_sa.k;
              seeds[i][n_seeds[i]].l = match_sa.l;
              seeds[i][n_seeds[i]].offset = j;
              n_seeds[i]++;
              m_hits[i] += match_sa.l - match_sa.k + 1;
          }
      }
  }

  // convert seeds to chr/pos
  for(i=0;i<2;i++) { // forward/reverse-compliment
      // allocate mmemory
      hits[i] = fmap_malloc(m_hits[i]*sizeof(fmap_map3_aux_hit_t), "hits[i]");
      n_hits[i] = 0;

      for(j=0;j<n_seeds[i];j++) { // go through all seeds
          for(k=seeds[i][j].k;k<=seeds[i][j].l;k++) { // through all occurrences
              pacpos = fmap_sa_pac_pos(sa, bwt, k);
              if(bwt->seq_len < pacpos + seeds[i][j].offset + 1) { // before the beginning of the reference sequence
                  pacpos = 0;
              }
              else {
                  pacpos = bwt->seq_len - pacpos - seeds[i][j].offset + 1;
              }
              if(0 < fmap_refseq_pac2real(refseq, pacpos, 1,
                                          &hits[i][n_hits[i]].seqid, &hits[i][n_hits[i]].pos)) {
                  hits[i][n_hits[i]].pos -= opt->seed_length;
                  hits[i][n_hits[i]].offset = seeds[i][j].offset;
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

          // get the target range
          if(hits[i][start].pos < opt->sw_offset) ref_start = 1;
          else ref_start = hits[i][start].pos - opt->sw_offset;
          ref_end = hits[i][end].pos + seq_len[i] + opt->sw_offset - 1;
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
          par.band_width += 2 * opt->sw_offset; // add bases to the window

          path = fmap_calloc(2*target_len, sizeof(fmap_sw_path_t), "path");

          // threshold the score by assuming that one seed's worth of
          // matches occurs in the alignment
          score = fmap_sw_local_core(target, target_len, query, seq_len[i], &par, path, &path_len, opt->seed_length * opt->score_match, &score_subo);

          // i - target
          // j - query

          if(0 < score) {
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

          // free
          free(path);

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

  return aln;
}
