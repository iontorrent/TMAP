/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include "../util/tmap_alloc.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_sa.h"
#include "../index/tmap_bwt_match.h"
#include "../sw/tmap_sw.h"
#include "tmap_map_util.h"
#include "tmap_map3.h"
#include "tmap_map3_aux.h"

static inline void 
tmap_map3_aux_seed_add(tmap_map3_aux_seed_t **seeds,
                       int32_t *n_seeds,
                       int32_t *m_seeds,
                       int32_t k,
                       int32_t l,
                       int32_t start, 
                       int32_t n_diff,
                       int32_t seed_step)
{
  /*
  if(offset < INT8_MIN || INT8_MAX < offset) {
      tmap_error("offset for hp enumeration was out of range", Warn, OutOfRange);
  }
  */
  if((*m_seeds) <= (*n_seeds)) {
      (*m_seeds) = (0 == (*m_seeds)) ? 64 : ((*m_seeds) << 1);
      (*seeds) = tmap_realloc((*seeds), sizeof(tmap_map3_aux_seed_t)*(*m_seeds), "(*seeds)");
  }
  (*seeds)[(*n_seeds)].k = k;
  (*seeds)[(*n_seeds)].l = l;
  (*seeds)[(*n_seeds)].start = start;
  (*seeds)[(*n_seeds)].n_diff = n_diff;
  (*seeds)[(*n_seeds)].seed_step = seed_step;
  (*n_seeds)++;
}

static inline void
tmap_map3_aux_core_seed_helper(uint8_t *query,
                               int32_t query_length,
                               int32_t offset,
                               uint8_t *flow_order,
                               uint8_t flow_i,
                               tmap_refseq_t *refseq,
                               tmap_bwt_t *bwt,
                               tmap_sa_t *sa,
                               tmap_map_opt_t *opt,
                               tmap_map3_aux_seed_t **seeds,
                               int32_t *n_seeds,
                               int32_t *m_seeds,
                               int32_t seed_length)
{ 
  int32_t i, k;
  int32_t n_bases;;
  tmap_bwt_match_occ_t prev_sa, cur_sa, next_sa, tmp_sa;

  if(query_length <= offset) return;
  if(flow_order[flow_i] != query[offset]) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }

  // initialize prev
  prev_sa.k = 0;
  prev_sa.l = bwt->seq_len;
  prev_sa.hi = 0;
  prev_sa.offset = 0;

  i = offset;
  while(i < query_length) {
      // reached the seed length
      if(seed_length < i - offset) { 
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
          if(seed_length < i - offset + k) { 
              break;
          }

          // only delete if there are bases available and we are not deleting
          // the entire first flow 
          int32_t bases_to_align = seed_length - (i - offset + k);
          int32_t bases_left = query_length - i - n_bases;
          if(0 < n_bases - k // bases to delete
             && n_bases - k <= opt->hp_diff // not too many to delete
             && (i != offset || 0 != k) // do not delete the entire flow
             && bases_to_align <= bases_left) { // enough bases 
              // match exactly from here onwards
              tmp_sa = next_sa;
              if(0 < tmap_bwt_match_exact_alt(bwt, bases_to_align, query + i + n_bases, &tmp_sa)
                 && (tmp_sa.l - tmp_sa.k + 1) <= opt->max_seed_hits) {
                  tmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, tmp_sa.k, tmp_sa.l, offset, n_bases-k, 0);
              }
          }

          // move past this base in the hp
          tmp_sa = next_sa;
          tmap_bwt_match_2occ(bwt, &tmp_sa, flow_order[flow_i], &next_sa);
          if(next_sa.l < next_sa.k) { // no match, return
              return;
          }
      }

      // insert hp bases
      if(i + n_bases < offset + seed_length) { // not the last flow
          cur_sa = next_sa; // already considered the 'n_bases' of this flow
          // insert
          for(k=1;k<=opt->hp_diff;k++) { // # of bases to insert
              tmap_bwt_match_2occ(bwt, &cur_sa, flow_order[flow_i], &tmp_sa);
              if(tmp_sa.l < tmp_sa.k) { // no match, do not continue
                  break;
              }
              // match exactly from here onwards
              if(0 < tmap_bwt_match_exact_alt(bwt, seed_length - (i - offset) - n_bases, query + i + n_bases, &tmp_sa)
                 && (tmp_sa.l - tmp_sa.k + 1) <= opt->max_seed_hits) {
                  tmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, tmp_sa.k, tmp_sa.l, offset, k, 0);
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

  if(i - offset < seed_length) tmap_error("bug encountered", Exit, OutOfRange);

  // add in the seed with no hp indels
  if((next_sa.l - next_sa.k + 1) <= opt->max_seed_hits) {
      tmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, next_sa.k, next_sa.l, offset, 0, 0);
  }
}

static inline void
tmap_map3_aux_core_seed(uint8_t *query,
                        int32_t query_length,
                        uint8_t *flow_order,
                        tmap_refseq_t *refseq,
                        tmap_bwt_t *bwt,
                        tmap_sa_t *sa,
                        tmap_map_opt_t *opt,
                        tmap_map3_aux_seed_t **seeds,
                        int32_t *n_seeds,
                        int32_t *m_seeds,
                        int32_t seed_length,
                        int32_t seed_step)
{
  int32_t i, j, flow_i;

  if(0 < opt->hp_diff) {
      i=flow_i=0;
      while(i<query_length - seed_length + 1) {
          // move to the next flow
          j=0;
          while(query[i] != flow_order[flow_i]) {
              flow_i = (flow_i + 1) & 3;
              // sanity check
              j++;
              if(4 <= j) tmap_error("bug encountered", Exit, OutOfRange);
          }

          // add seeds
          tmap_map3_aux_core_seed_helper(query, query_length, i, flow_order, flow_i,
                                         refseq, bwt, sa, opt, seeds, n_seeds, m_seeds,
                                         seed_length);

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
      int count;
      tmap_bwt_match_occ_t cur_sa;
      j = count = 0;
      for(i=query_length-seed_length;0<=i;i--) {
          if(0 < tmap_bwt_match_exact(bwt, seed_length, query + i, &cur_sa)) {
              count++;
              if((cur_sa.l - cur_sa.k + 1) <= opt->max_seed_hits) {
                  tmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, cur_sa.k, cur_sa.l, i, 0, 0);
                  j++;
              }
              else {
                  // seed stepping
                  if(0 < seed_step) {
                      int32_t k = i + seed_length;
                      int32_t n = 0;
                      while(k + seed_step < query_length && 0 < tmap_bwt_match_exact(bwt, seed_step, query + k, &cur_sa)) {
                          if((cur_sa.l - cur_sa.k + 1) <= opt->max_seed_hits) {
                              tmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, cur_sa.k, cur_sa.l, i, 0, n);
                          }
                          k += seed_step;
                          n++;
                      }
                  }
              }
          }
          else {
              // skip over if we came up short
              i -= (seed_length - cur_sa.offset); 
          }
      }
      // remove seeds if there were too many repetitive hits
      // NB: does not count seed steps
      if(j / (double)count < opt->hit_frac) {
          (*n_seeds) -= j;
      }
      /*
      for(i=0;i<query_length-seed_length+1;i++) {
          if(0 < tmap_bwt_match_exact(bwt, seed_length, query + i, &cur_sa)
                 && (cur_sa.l - cur_sa.k + 1) <= opt->max_seed_hits) {
              tmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, cur_sa.k, cur_sa.l, i, 0);
          }
      }
      */
  }
}

// TODO: memory pools?
tmap_map_sams_t *
tmap_map3_aux_core(tmap_seq_t *seq[2], 
                   uint8_t *flow_order,
                   int32_t flow_order_len,
                   tmap_refseq_t *refseq,
                   tmap_bwt_t *bwt,
                   tmap_sa_t *sa,
                   tmap_map_opt_t *opt)
{
  int32_t i, j, n, seed_length;
  int32_t seq_len[2];
  tmap_string_t *bases;
  uint8_t *query;
  uint8_t *flow[2];

  tmap_map3_aux_seed_t *seeds[2];
  int32_t m_seeds[2], n_seeds[2];

  tmap_map_sams_t *sams = NULL;

  if(0 < opt->hp_diff) {
      // set up the flow order to be used
      if(NULL == flow_order) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }
      else {
          flow[0] = tmap_malloc(sizeof(uint8_t)*flow_order_len, "flow[0]");
          flow[1] = tmap_malloc(sizeof(uint8_t)*flow_order_len, "flow[0]");
          for(i=0;i<flow_order_len;i++) {
              flow[0][i] = flow_order[i]; // forward
              flow[1][flow_order_len-1-i] = 3 - flow_order[i]; // reverse complimenmt
          }
      }
  }

  // init
  sams = tmap_map_sams_init(NULL);

  // update the seed length based on the read length
  seed_length = opt->seed_length;
  if(0 == opt->seed_length_set) {
      i = tmap_seq_get_bases(seq[0])->l;
      while(0 < i) {
          seed_length++;
          i >>= 1; // divide by two
      }
  }
  
  // check that the sequences are long enough
  for(i=0;i<2;i++) { // forward/reverse-compliment
      bases = tmap_seq_get_bases(seq[i]);
      seq_len[i] = bases->l;
      if(seq_len[i] - seed_length < 0) {
          return sams;
      }
  }

  // seeds
  for(i=n=0;i<2;i++) { // forward/reverse-compliment
      bases = tmap_seq_get_bases(seq[i]);
      seq_len[i] = bases->l;
      query = (uint8_t*)bases->s;

      n_seeds[i] = 0;
      // pre-allocate mmemory
      if(0 < opt->seed_step) {
          m_seeds[i] = 0;
          int32_t k;
          for(j=seed_length,k=0;j<=seq_len[i];j+=opt->seed_step,k++) {
              if(UINT8_MAX < k) {
                  tmap_error("seed step out of range", Warn, OutOfRange);
                  break;
              }
              m_seeds[i] += seq_len[i] - j + 1; // maximum number of seeds possible
          }
      }
      else {
          m_seeds[i] = seq_len[i] - seed_length + 1; // maximum number of seeds possible
      }
      seeds[i] = tmap_malloc(m_seeds[i]*sizeof(tmap_map3_aux_seed_t), "seeds[i]");

      // seed the alignment
      tmap_map3_aux_core_seed(query, seq_len[i], flow[i],
                              refseq, bwt, sa, opt, &seeds[i], &n_seeds[i], &m_seeds[i],
                              seed_length, opt->seed_step);

      // for SAM storage
      for(j=0;j<n_seeds[i];j++) {
          n += seeds[i][j].l - seeds[i][j].k + 1;
      }
  }

  // make enough room
  tmap_map_sams_realloc(sams, n);
  
  // convert seeds to chr/pos
  for(i=n=0;i<2;i++) { // forward/reverse-compliment
      for(j=0;j<n_seeds[i];j++) { // go through all seeds
          uint32_t seqid, pos;
          uint32_t k, pacpos;
          uint8_t seed_length_ext = seed_length + (seeds[i][j].seed_step * opt->seed_length);
          for(k=seeds[i][j].k;k<=seeds[i][j].l;k++) { // through all occurrences
              tmap_map_sam_t *s = NULL;
              pacpos = tmap_sa_pac_pos(sa, bwt, k);
              if(bwt->seq_len < pacpos + seeds[i][j].start + 1) { // before the beginning of the reference sequence
                  pacpos = 0;
              }
              else {
                  pacpos = bwt->seq_len - pacpos - seeds[i][j].start - 1;
              }
              if(0 < tmap_refseq_pac2real(refseq, pacpos, 1, &seqid, &pos)) {
                  uint32_t tmp_seqid, tmp_pos;
                  // we need to check if we crossed contig boundaries
                  if(0 < tmap_refseq_pac2real(refseq, pacpos + seeds[i][j].start, 1, &tmp_seqid, &tmp_pos) &&
                     seqid < tmp_seqid) {
                      seqid = tmp_seqid;
                      pos = 0;
                  }
                  else {
                      if(pos < seed_length_ext + seeds[i][j].n_diff - 1 ) {
                          pos = 0;
                      }
                      else {
                          pos -= seed_length_ext + seeds[i][j].n_diff - 1;
                      }
                  }

                  // save
                  s = &sams->sams[n];

                  // save the hit
                  s->algo_id = TMAP_MAP_ALGO_MAP3;
                  s->algo_stage = 0;
                  s->strand = i;
                  s->seqid = seqid;
                  s->pos = pos;
                  s->target_len = seq_len[i]; 
                  if(refseq->annos[seqid].len < s->target_len) {
                      s->target_len = refseq->annos[seqid].len;
                  }
                  s->score_subo = INT32_MIN;

                  // map3 aux data
                  tmap_map_sam_malloc_aux(s, TMAP_MAP_ALGO_MAP3);

                  n++;
              }
          }
      }
  }

  // resize
  tmap_map_sams_realloc(sams, n);

  // free the seeds
  for(i=0;i<2;i++) {
      free(seeds[i]);
      seeds[i]=NULL;
  }

  // free
  if(0 < opt->hp_diff) {
      free(flow[0]);
      free(flow[1]);
  }

  return sams;
}
