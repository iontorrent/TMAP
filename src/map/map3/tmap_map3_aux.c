/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include "../../util/tmap_alloc.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_sa.h"
#include "../../index/tmap_index.h"
#include "../../index/tmap_bwt_match.h"
#include "../../index/tmap_bwt_match_hash.h"
#include "../../sw/tmap_sw.h"
#include "../util/tmap_map_util.h"
#include "tmap_map3.h"
#include "tmap_map3_aux.h"

static inline void 
tmap_map3_aux_seed_add(tmap_map3_aux_seed_t **seeds,
                       int32_t *n_seeds,
                       int32_t *m_seeds,
                       tmap_bwt_int_t k,
                       tmap_bwt_int_t l,
                       int32_t start, 
                       int16_t seed_length)
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
  (*seeds)[(*n_seeds)].seed_length = seed_length;
  (*n_seeds)++;
}

static inline void
tmap_map3_aux_core_seed(uint8_t *query,
                        int32_t query_length,
                        tmap_refseq_t *refseq,
                        tmap_bwt_t *bwt,
                        tmap_sa_t *sa,
                        tmap_bwt_match_hash_t *hash,
                        tmap_map_opt_t *opt,
                        tmap_map3_aux_seed_t **seeds,
                        int32_t *n_seeds,
                        int32_t *m_seeds,
                        int32_t seed_length,
                        int32_t seed_step,
                        int32_t fwd_search)
{
  int32_t i, j;

  int k, count;
  tmap_bwt_match_occ_t cur_sa, prev_sa;
  j = count = 0;
  if(1 == fwd_search) {
      for(i=0;i<query_length-seed_length+1;i++) {
          if(0 < tmap_bwt_match_hash_exact(bwt, seed_length, query + i, &cur_sa, hash)) {
              count++;
              if((cur_sa.l - cur_sa.k + 1) <= opt->max_seed_hits) {
                  // extend further
                  prev_sa = cur_sa;
                  k = i + 1;
                  while(k < query_length - seed_length) {
                      tmap_bwt_match_hash_2occ(bwt, &prev_sa, query[k], &cur_sa, hash);
                      if(cur_sa.l < cur_sa.k) {
                          // use prev
                          cur_sa = prev_sa;
                          break;
                      }
                      else {
                          // keep going
                          prev_sa = cur_sa;
                          k++;
                      }
                  }
                  k--; // k is always one greater
                  tmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, cur_sa.k, cur_sa.l, k, seed_length + k - i);
                  j++;
                  // skip over 
                  if(0 < opt->skip_seed_frac) {
                      i += opt->skip_seed_frac * (seed_length + k - i - 1); // - 1 since i will be incremented
                  }
              }
              else {
                  // seed stepping
                  if(0 < seed_step) {
                      k = i + seed_length;
                      int32_t n = 0;
                      while(k + seed_step < query_length && 0 < tmap_bwt_match_hash_exact(bwt, seed_step, query + k, &cur_sa, hash)) {
                          if((cur_sa.l - cur_sa.k + 1) <= opt->max_seed_hits) {
                              tmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, cur_sa.k, cur_sa.l, i, seed_length + k - i);
                              j++;
                              // skip over 
                              if(0 < opt->skip_seed_frac) {
                                  i += opt->skip_seed_frac * (seed_length + k - i - 1); // - 1 since i will be incremented
                              }
                              break;
                          }
                          k += seed_step;
                          n++;
                      }
                  }
              }
          }
      }
  }
  else {
      for(i=query_length-seed_length;0<=i;i--) {
          if(0 < tmap_bwt_match_hash_exact(bwt, seed_length, query + i, &cur_sa, hash)) {
              count++;
              if((cur_sa.l - cur_sa.k + 1) <= opt->max_seed_hits) {
                  tmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, cur_sa.k, cur_sa.l, i, seed_length);
                  j++;
                  if(0 < opt->skip_seed_frac) {
                      i -= opt->skip_seed_frac * (seed_length - 1); // -1 since i will be incremented
                  }
              }
              else {
                  // seed stepping
                  if(0 < seed_step) {
                      int32_t k = i + seed_length;
                      int32_t n = 0;
                      while(k + seed_step < query_length && 0 < tmap_bwt_match_hash_exact_alt(bwt, seed_step, query + k, &cur_sa, hash)) {
                          if((cur_sa.l - cur_sa.k + 1) <= opt->max_seed_hits) {
                              tmap_map3_aux_seed_add(seeds, n_seeds, m_seeds, cur_sa.k, cur_sa.l, i, seed_length + k - i);
                              j++;
                              if(0 < opt->skip_seed_frac) {
                                  i -= opt->skip_seed_frac * (seed_length + k - i - 1); // -1 since i will be incremented
                              }
                              // break when the e find the first hit
                              break;
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
  }
  // remove seeds if there were too many repetitive hits
  // NB: does count seed steps
  if(j / (double)count < opt->hit_frac) {
      (*n_seeds) = 0;
      //(*n_seeds) -= j;
  }
}

// TODO: memory pools?
tmap_map_sams_t *
tmap_map3_aux_core(tmap_seq_t *seq, 
                   tmap_refseq_t *refseq,
                   tmap_bwt_t *bwt,
                   tmap_sa_t *sa,
                   tmap_bwt_match_hash_t *hash,
                   tmap_map_opt_t *opt)
{
  int32_t i, j, n, seed_length;
  int32_t seq_len;
  tmap_string_t *bases;
  uint8_t *query;
  tmap_map3_aux_seed_t *seeds;
  int32_t m_seeds, n_seeds;
  tmap_map_sams_t *sams = NULL;

  // init
  sams = tmap_map_sams_init(NULL);

  // update the seed length based on the read length
  seed_length = opt->seed_length;
  if(0 == opt->seed_length_set) {
      i = tmap_seq_get_bases_length(seq);
      while(0 < i) {
          seed_length++;
          i >>= 1; // divide by two
      }
  }
  
  // check that the sequence is long enough
  bases = tmap_seq_get_bases(seq);
  seq_len = bases->l;
  if(seq_len - seed_length < 0) {
      return sams;
  }

  // seeds
  bases = tmap_seq_get_bases(seq);
  query = (uint8_t*)bases->s;

  n_seeds = 0;
  // pre-allocate mmemory
  if(0 < opt->seed_step) {
      m_seeds = 0;
      int32_t k;
      for(j=seed_length,k=0;j<=seq_len;j+=opt->seed_step,k++) {
          if(UINT16_MAX < k) {
              /*
              fprintf(stderr, "j=%d seed_length=%d k=%d seq_len=%d opt->seed_step=%d\n",
                      j, seed_length, k, seq_len, opt->seed_step);
                      */
              tmap_error("seed step out of range", Warn, OutOfRange);
              return sams;
          }
          m_seeds += seq_len - j + 1; // maximum number of seeds possible
      }
  }
  else {
      m_seeds = seq_len - seed_length + 1; // maximum number of seeds possible
  }
  seeds = tmap_malloc(m_seeds*sizeof(tmap_map3_aux_seed_t), "seeds");

  // seed the alignment
  tmap_map3_aux_core_seed(query, seq_len, refseq, bwt, sa, 
                          hash, opt, &seeds, &n_seeds, &m_seeds,
                          seed_length, opt->seed_step, opt->fwd_search);
  if(0 == n_seeds) { // try the seeding in the opposite direction
      tmap_map3_aux_core_seed(query, seq_len, refseq, bwt, sa, 
                              hash, opt, &seeds, &n_seeds, &m_seeds,
                              seed_length, opt->seed_step, 1-opt->fwd_search);
  }

  // for SAM storage
  n = 0;
  for(j=0;j<n_seeds;j++) {
      n += seeds[j].l - seeds[j].k + 1;
  }

  // make enough room
  tmap_map_sams_realloc(sams, n);
  
  // convert seeds to chr/pos
  n = 0;
  for(j=0;j<n_seeds;j++) { // go through all seeds
      uint32_t seqid, pos, pos_adj;
      tmap_bwt_int_t k, pacpos;
      uint16_t seed_length_ext = seeds[j].seed_length;
      uint16_t start = seeds[j].start;
      uint8_t strand;
      for(k=seeds[j].k;k<=seeds[j].l;k++) { // through all occurrences
          tmap_map_sam_t *s = NULL;
          pacpos = bwt->seq_len - tmap_sa_pac_pos_hash(sa, bwt, k, hash);
          if(0 < tmap_refseq_pac2real(refseq, pacpos, 1, &seqid, &pos, &strand)) {
              // adjust based on offset
              pos_adj = start + seed_length_ext - 1; // NB: we only used the forward read, so adjustment is based off of this...
              /*
              fprintf(stderr, "seqid=%u pos=%u pos_adj=%u start=%u seed_length_ext=%u strand=%u n=%llu\n",
                      seqid, pos, pos_adj, start, seed_length_ext, strand, seeds[j].l - seeds[j].k + 1);
                      */
              // contig boundary
              if(pos <= pos_adj) pos = 0; 
              else pos -= pos_adj;

              // save
              s = &sams->sams[n];
              tmap_map_sam_init(s);

              // save the hit
              s->algo_id = TMAP_MAP_ALGO_MAP3;
              s->algo_stage = opt->algo_stage;
              s->strand = strand;
              s->seqid = seqid;
              s->pos = pos;
              s->target_len = seq_len; 
              if(refseq->annos[seqid].len < s->target_len) {
                  s->target_len = refseq->annos[seqid].len;
              }
              s->score_subo = INT32_MIN;

              // map3 aux data
              tmap_map_sam_malloc_aux(s);

              n++;
          }
      }
  }

  // resize
  tmap_map_sams_realloc(sams, n);

  // free the seeds
  free(seeds);
  seeds=NULL;

  return sams;
}
