/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_fibheap.h"
#include "../util/tmap_definitions.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_bwt_match.h"
#include "../index/tmap_sa.h"
#include "../sw/tmap_sw.h"
#include "tmap_map_util.h"
#include "tmap_map1.h"
#include "tmap_map1_aux.h"

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

// TODO redefine
#define aln_score(m,o,e,p) ((m)*(p)->pen_mm + (o)*((p)->pen_gapo + (p)->pen_gape) + (e)*(p)->pen_gape)

#define tmap_map1_aux_reverse_query(_query, _ql, _i) \
  for(_i=0;_i<(_ql>>1);_i++) { \
      uint8_t _tmp = _query[_i]; \
      _query[_i] = _query[_ql-_i-1]; \
      _query[_ql-_i-1] = _tmp; \
  }

static int 
tmap_map1_aux_stack_cmp(void *a, void *b)
{
  tmap_map1_aux_stack_entry_t *x = (tmap_map1_aux_stack_entry_t*)a;
  tmap_map1_aux_stack_entry_t *y = (tmap_map1_aux_stack_entry_t*)b;

  /* sort by score, then offset */
  if(x->score < y->score 
     || (x->score == y->score && x->offset > y->offset)) {
      return -1;
  }
  else if(x->score == y->score && x->offset == y->offset) {
      return 0;
  }
  else {
      return 1;
  }
}

tmap_map1_aux_stack_t *
tmap_map1_aux_stack_init()
{
  int32_t i;
  tmap_map1_aux_stack_t *stack = NULL;
  stack = tmap_calloc(1, sizeof(tmap_map1_aux_stack_t), "stack");
  stack->entry_pool_length = 1024; // TODO: move to a define block
  stack->entry_pool = tmap_malloc(stack->entry_pool_length*sizeof(tmap_map1_aux_stack_entry_t*), "stack->entry_pool");
  for(i=0;i<stack->entry_pool_length;i++) {
      stack->entry_pool[i] = tmap_malloc(sizeof(tmap_map1_aux_stack_entry_t), "stack->entry_pool[i]");
  }
  stack->heap = tmap_fibheap_makeheap(tmap_map1_aux_stack_cmp);
  //stack->heap = tmap_fibheap_makekeyheap();

  return stack;
}

void
tmap_map1_aux_stack_destroy(tmap_map1_aux_stack_t *stack)
{
  int32_t i;
  tmap_fibheap_deleteheap(stack->heap);
  for(i=0;i<stack->entry_pool_length;i++) {
      free(stack->entry_pool[i]);
  }
  free(stack->entry_pool);
  free(stack);
}

static void
tmap_map1_aux_stack_reset(tmap_map1_aux_stack_t *stack)
{
  // reset heap
  tmap_fibheap_deleteheap(stack->heap);
  stack->heap = tmap_fibheap_makeheap(tmap_map1_aux_stack_cmp);
  //stack->heap = tmap_fibheap_makekeyheap();
  // move to the beginning of the memory pool
  stack->entry_pool_i = 0;
  stack->best_score = INT32_MAX;
}

static inline void
tmap_map1_aux_stack_push(tmap_map1_aux_stack_t *stack, int32_t strand, 
                         int32_t offset,
                         tmap_bwt_match_occ_t *match_sa_prev,
                         int32_t n_mm, int32_t n_gapo, int32_t n_gape,
                         int32_t state, int32_t is_diff, 
                         tmap_map1_aux_stack_entry_t *prev_entry,
                         const tmap_map_opt_t *opt)
{
  tmap_map1_aux_stack_entry_t *entry = NULL;

  // check to see if we need more memory
  if(stack->entry_pool_length <= stack->entry_pool_i) { 
      int32_t i = stack->entry_pool_length;
      stack->entry_pool_length <<= 2;
      stack->entry_pool = tmap_realloc(stack->entry_pool, 
                                       sizeof(tmap_map1_aux_stack_entry_t*)*stack->entry_pool_length, 
                                       "stack->entry_pool");
      while(i<stack->entry_pool_length) {
          stack->entry_pool[i] = tmap_malloc(sizeof(tmap_map1_aux_stack_entry_t), "stack->entry_pool[i]");
          i++;
      }
  }

  entry = stack->entry_pool[stack->entry_pool_i];
  entry->score = aln_score(n_mm, n_gapo, n_gape, opt);
  entry->n_mm = n_mm;
  entry->n_gapo = n_gapo;
  entry->n_gape = n_gape;
  entry->state = state;
  entry->strand = strand;
  entry->match_sa = (*match_sa_prev); 
  entry->i = stack->entry_pool_i;
  entry->offset = offset;
  if(NULL == prev_entry) {
      entry->last_diff_offset = (1 == is_diff) ? offset : 0;
      entry->prev_i = -1;
  }
  else {
      entry->last_diff_offset = (1 == is_diff) ? offset : prev_entry->last_diff_offset;
      entry->prev_i = prev_entry->i;
  }

  if(stack->best_score > entry->score) stack->best_score = entry->score;

  tmap_fibheap_insert(stack->heap, entry);
  //tmap_fibheap_insertkey(stack->heap, entry->score, entry);

  stack->entry_pool_i++;
}

static inline tmap_map1_aux_stack_entry_t *
tmap_map1_aux_stack_pop(tmap_map1_aux_stack_t *stack)
{
  tmap_map1_aux_stack_entry_t *best, *next_best;

  best = tmap_fibheap_extractmin(stack->heap);
  next_best = tmap_fibheap_min(stack->heap);
  if(NULL != next_best) {
      stack->best_score  = next_best->score;
  }
  else {
      stack->best_score = INT32_MAX;
  }

  return best;
}

static inline void 
tmap_map1_aux_stack_shadow(int32_t x, int32_t len, uint32_t max, 
                           int32_t last_diff_offset, tmap_bwt_match_width_t *w)
{
  int32_t i, j;
  for(i=j=len-1;last_diff_offset<i;i--) {
      if(x < w[i].w) { 
          w[i].w -= x;
      }
      else if(x == w[i].w) {
          w[i].bid = 1;
          j++;
          w[i].w = max - j;
      }
      else {
          tmap_error("else happened", Exit, OutOfRange);
      }
  }
}

static inline int32_t
tmap_map1_aux_stack_size(tmap_map1_aux_stack_t *stack)
{
  return stack->heap->tmap_fibheap_n;
}

static inline int32_t int_log2(uint32_t v)
{
  int32_t c = 0;
  if(v & 0xffff0000u) { v >>= 16; c |= 16; }
  if(v & 0xff00) { v >>= 8; c |= 8; }
  if(v & 0xf0) { v >>= 4; c |= 4; }
  if(v & 0xc) { v >>= 2; c |= 2; }
  if(v & 0x2) c |= 1;
  return c;
}

static inline int
tmap_map1_aux_get_bam_state(int state)
{
  switch(state) {
    case STATE_M:
      return BAM_CMATCH;
    case STATE_I:
      return BAM_CINS;
    case STATE_D:
      return BAM_CDEL;
    default:
      // in theory this should not occur
      return BAM_CMATCH;
  }
}

#define __add_to_cigar() do { \
    if(0 < op_len) { \
        if(sam->n_cigar <= cigar_i) { \
            sam->n_cigar++; \
            sam->cigar = tmap_realloc(sam->cigar, sizeof(uint32_t)*sam->n_cigar, "sam->cigar"); \
        } \
        TMAP_SW_CIGAR_STORE(sam->cigar[cigar_i], tmap_map1_aux_get_bam_state(op), op_len); \
        cigar_i++; \
    } \
} while(0)

static tmap_map_sams_t *
tmap_map1_sam_to_real(tmap_map_sams_t *sams, tmap_string_t *bases[2], int32_t seed2_len,
                       tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2], tmap_map_opt_t *opt) 
{
  tmap_map_sams_t *sams_tmp = NULL;
  uint32_t i, j, k, l, m, n;
  int32_t matrix[25];
  uint8_t *target = NULL;
  int32_t target_length;
  tmap_sw_param_t par;
  tmap_sw_path_t *path = NULL;
  int32_t path_len, path_mem=0;
  int32_t seq_len; 
  uint8_t *query=NULL;

  seq_len = bases[0]->l;
  par.matrix = matrix;
  __gen_ap(par, opt);
  target_length = ((seq_len + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + seq_len;
  target = tmap_calloc(target_length, sizeof(uint8_t), "target");

  // max # of entries
  for(i=n=0;i<sams->n;i++) {
      n += sams->sams[i].pos - sams->sams[i].seqid + 1; // l - k + 1
  }

  // alloc
  sams_tmp = tmap_map_sams_init();
  tmap_map_sams_realloc(sams_tmp, n);
            
  // reverse (it has been complimented)
  tmap_string_reverse(bases[1]);

  // copy over
  for(i=j=0;i<sams->n;i++) {
      tmap_map_sam_t *sam;

      sam = &sams->sams[i];

      // go through SA interval
      for(k=sams->sams[i].seqid;k<=sams->sams[i].pos;k++) { // k -> l
          uint32_t pos = 0, seqid = 0, pacpos = 0, lt;
          int32_t aln_ref_l, added;
          uint8_t strand;

          lt = target_length;
          strand = sams->sams[i].strand;

          // query sequence
          query = (uint8_t*)bases[sam->strand]->s;
          // get the target sequence to which we will align
          if(0 == strand) { // forward
              pacpos = bwt[1-strand]->seq_len - tmap_sa_pac_pos(sa[1-strand], bwt[1-strand], k) - sam->aux.map1_aux->aln_ref; // pacpos is zero-based
          }
          else { // reverse
              pacpos = tmap_sa_pac_pos(sa[1-strand], bwt[1-strand], k) + sam->aux.map1_aux->aln_ref;
              pacpos = (pacpos < seq_len) ? 0 : (pacpos - seq_len); // guard against going off the edge of the reference
          }
          for(l = pacpos, m = 0; l < pacpos + lt && l < refseq->len; l++) {
              target[m++] = tmap_refseq_seq_i(refseq, l); 
          }
          lt = m;

          /*
          fprintf(stderr, "Q=");
          for(l=0;l<seq_len;l++) {
              fputc("ACGTN"[query[l]], stderr);
          }
          fputc('\n', stderr);
          fprintf(stderr, "T=");
          for(l=0;l<lt;l++) {
              fputc("ACGTN"[target[l]], stderr);
          }
          fputc('\n', stderr);
          */

          //fprintf(stderr, "i=%d j=%d k=%d lt=%d pacpos=%d seq_len=%d\n", i, j, k, lt, pacpos, seq_len);

          // get more memory if required
          if(path_mem <= lt + seq_len) { 
              path_mem = lt + seq_len;
              tmap_roundup32(path_mem);
              path = tmap_realloc(path, sizeof(tmap_sw_path_t)*path_mem, "path");
          }
          // align
          added = tmap_map_util_sw(&sams_tmp->sams[j],
                                   target, lt, 
                                   query, seq_len, 
                                   seqid, 0,
                                   &par, path, &path_len, 
                                   opt->score_thr, opt->softclip_type, strand);
          if(0 == path_len) {
              tmap_error("0 == path_len", Exit, OutOfRange);
          }

          // adjust pacpos
          pacpos += sams_tmp->sams[j].pos + 1; // now pacpos is one-based
          aln_ref_l = path[0].i - path[path_len-1].i + 1; 

          // save the hit
          if(1 == added && 0 < tmap_refseq_pac2real(refseq, pacpos, aln_ref_l, &seqid, &pos)) {
              // copy over previous parameters
              sams_tmp->sams[j].algo_id = TMAP_MAP_ALGO_MAP1;
              sams_tmp->sams[j].algo_stage = 0;

              sams_tmp->sams[j].pos = pos-1; // adjust
            
              // aux
              tmap_map_sam_malloc_aux(&sams_tmp->sams[j], TMAP_MAP_ALGO_MAP1);
              sams_tmp->sams[j].aux.map1_aux->n_mm = sam->aux.map1_aux->n_mm;
              sams_tmp->sams[j].aux.map1_aux->n_gapo = sam->aux.map1_aux->n_gapo;
              sams_tmp->sams[j].aux.map1_aux->n_gape = sam->aux.map1_aux->n_gape;
              // get the number of non-inserted bases 
              sams_tmp->sams[j].aux.map1_aux->aln_ref = 0;
              for(l=0;l<sams_tmp->sams[j].n_cigar;l++) {
                  switch(TMAP_SW_CIGAR_OP(sams_tmp->sams[j].cigar[l])) {
                    case BAM_CMATCH:
                    case BAM_CDEL:
                      sams_tmp->sams[j].aux.map1_aux->aln_ref += TMAP_SW_CIGAR_LENGTH(sams_tmp->sams[j].cigar[l]); 
                      break;
                    default:
                      break;
                  }
              }
              j++;
          }
      }
  }
  
  // reverse back
  tmap_string_reverse(bases[1]);

  // destroy
  tmap_map_sams_destroy(sams);

  // realloc
  tmap_map_sams_realloc(sams_tmp, j);

  // free memory
  free(path);
  free(target);
  
  return sams_tmp;
}

tmap_map_sams_t *
tmap_map1_aux_core(tmap_seq_t *seq[2], tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2],
                   tmap_bwt_match_width_t *width[2], tmap_bwt_match_width_t *seed_width[2], tmap_map_opt_t *opt,
                   tmap_map1_aux_stack_t *stack, int32_t seed2_len)
{
  int32_t max_mm = opt->max_mm, max_gapo = opt->max_gapo, max_gape = opt->max_gape, seed_max_mm = opt->seed_max_mm;
  int32_t best_score, next_best_score;
  int32_t best_cnt = 0;
  int32_t j, num_n = 0;
  int32_t max_edit_score;
  tmap_bwt_match_occ_t match_sa_start;
  tmap_string_t *bases[2]={NULL,NULL};
  tmap_map_sams_t *sams = NULL;
  
  sams = tmap_map_sams_init();

  best_score = next_best_score = aln_score(max_mm+1, max_gapo+1, max_gape+1, opt);

  max_edit_score = opt->pen_mm;
  if(max_edit_score < opt->pen_gapo + opt->pen_gape) max_edit_score = opt->pen_gapo;
  if(max_edit_score < opt->pen_gape) max_edit_score = opt->pen_gape;

  bases[0] = tmap_seq_get_bases(seq[0]);
  bases[1] = tmap_seq_get_bases(seq[1]);

  // check whether there are too many N
  for(j=num_n=0;j<bases[0]->l;j++) {
      if(3 < bases[0]->s[j]) {
          num_n++;
      }
  }
  if(max_mm < num_n) {
      return NULL;
  }

  match_sa_start.offset = 0;
  match_sa_start.hi = 0;
  match_sa_start.k = 0;
  match_sa_start.l = bwt[0]->seq_len;

  tmap_map1_aux_stack_reset(stack); // reset stack
  tmap_map1_aux_stack_push(stack, 0, 0, &match_sa_start, 0, 0, 0, STATE_M, 0, NULL, opt);
  tmap_map1_aux_stack_push(stack, 1, 0, &match_sa_start, 0, 0, 0, STATE_M, 0, NULL, opt);

  while(0 < tmap_map1_aux_stack_size(stack) && tmap_map1_aux_stack_size(stack) < opt->max_entries) {
      tmap_map1_aux_stack_entry_t *e = NULL;
      int32_t strand, len=-1; 
      int32_t n_mm, n_gapo, n_gape;
      int32_t n_seed_mm=0, offset;
      const uint8_t *str=NULL;
      int32_t sam_found;
      tmap_bwt_match_width_t *width_cur;
      const tmap_bwt_match_width_t *seed_width_cur = NULL;
      tmap_bwt_match_occ_t match_sa_cur, match_sa_next[4];

      e = tmap_map1_aux_stack_pop(stack); // get the best entry
      if(best_score + max_edit_score < e->score) break; // no need to continue

      strand = e->strand; // strand;
      match_sa_cur = e->match_sa; 
      n_mm = max_mm - e->n_mm;
      n_gapo = max_gapo - e->n_gapo;
      n_gape = max_gape - e->n_gape;
      if(n_mm < 0 || n_gapo < 0 || n_gape < 0) continue; // too many edits

      offset = e->offset; // zero-based
      str = (uint8_t*)bases[strand]->s;
      len = (0 < seed2_len) ? seed2_len : bases[strand]->l;
      //len = bases[strand]->l;
      width_cur = width[strand];

      if(NULL != seed_width && offset < opt->seed_length) { 
          seed_width_cur = seed_width[strand];
          n_seed_mm = seed_max_mm - e->n_mm; // consider only mismatches in the seed
      }
      else {
          seed_width_cur = NULL;
      }

      // if there are no more gaps, check if we are allowed mismatches
      if(0 == n_gapo && 0 == n_gape && offset < len && n_mm < width_cur[offset].bid) {
          continue;
      } 

      // check whether a sam is found
      sam_found = 0;
      if(len == offset) {
          sam_found = 1;
      }
      else if(0 == n_mm // no mismatches from any state
              && (e->state == STATE_M && 0 == n_gapo) // in STATE_M but no more gap opens
              && (e->state != STATE_M && 0 == n_gape)) { // in STATE_I/STATE_D but no more extensions
          if(0 < tmap_bwt_match_exact_alt(bwt[1-strand], offset, str, &match_sa_cur)) { // the alignment must match exactly to sam
              sam_found = 1;
          }
          else {
              continue; // no sam, skip
          }
      }

      if(1 == sam_found) { // alignment found
          int32_t score = aln_score(e->n_mm, e->n_gapo, e->n_gape, opt);
          int32_t do_add = 1;
          if(sams->n == 0) {
              best_score = score;
              best_cnt = 0;
          }
          if(score == best_score) {
              best_cnt += match_sa_cur.l - match_sa_cur.k + 1;
          }
          else {
              if(opt->max_best_cals < best_cnt) {
                  // ignore if too many "best" have been found
                  break;
              }
              if(score < next_best_score) {
                  next_best_score = score;
              }
              else if(next_best_score < score) {
                  // no need to examine further
                  break;
              }
          }
          if(0 < e->n_gapo) { // check if the same alignment has been found 
              for(j=0;j<sams->n;j++) {
                  if(sams->sams[j].seqid == match_sa_cur.k && sams->sams[j].pos == match_sa_cur.l) {
                      if(score < sams->sams[j].score) { // bug!
                          tmap_error("bug encountered", Exit, OutOfRange);
                      }
                      do_add = 0;
                      break;
                  }
              }
          }
          if(do_add) { // append
              uint32_t op, op_len, cigar_i;
              tmap_map_sam_t *sam = NULL;
              tmap_map1_aux_stack_entry_t *cur = e;
          
              tmap_map_sams_realloc(sams, sams->n+1);
              sam = &sams->sams[sams->n-1];

              sam->algo_id = TMAP_MAP_ALGO_MAP1;
              sam->algo_stage = 0;
              sam->score = cur->score;
              sam->strand = cur->strand;
              sam->seqid = cur->match_sa.k;
              sam->pos = cur->match_sa.l;
              
              // aux data
              tmap_map_sam_malloc_aux(sam, TMAP_MAP_ALGO_MAP1);
              sam->aux.map1_aux->n_mm = cur->n_mm;
              sam->aux.map1_aux->n_gapo = cur->n_gapo;
              sam->aux.map1_aux->n_gape = cur->n_gape;
              // aux data: reference length
              sam->aux.map1_aux->aln_ref = 0;
              cigar_i = 0;
              if(offset < len) { // we used 'tmap_bwt_match_exact_alt' 
                  op = STATE_M;
                  op_len = len - offset;
              }
              else {
                  op = -1;
                  op_len = 0;
              }
              while(0 < cur->offset) {
                  if(op != cur->state) {
                      if(STATE_M == op || STATE_D == op) {
                          sam->aux.map1_aux->aln_ref += op_len;
                      }
                      op = cur->state;
                      op_len = 1;
                  }
                  else {
                      op_len++;
                  }
                  cur = stack->entry_pool[cur->prev_i];
              }
              if(STATE_M == op || STATE_D == op) {
                  sam->aux.map1_aux->aln_ref += op_len;
              }

              // TODO: use the shadow ?
              //tmap_map1_aux_stack_shadow(l - k + 1, len, bwt[1-strand]->seq_len, e->last_diff_offset, width_cur);
          }
      }
      else {
          // use a bound for mismatches
          int32_t allow_diff = 1, allow_mm = 1;
          if(offset+1 < len) {
              if(n_mm-1 < width_cur[offset+1].bid) { 
                  allow_diff = 0;
              }
              else if(width_cur[offset].bid == n_mm-1
                      && width_cur[offset+1].bid == n_mm-1
                      && width_cur[offset+1].w == width_cur[offset].w) {
                  allow_mm = 0;
              }
              if(NULL != seed_width_cur && offset+1 < opt->seed_length) {
                  if(n_seed_mm < seed_width_cur[offset+1].bid) {
                      allow_diff = 0;
                  }
                  if(offset < len
                     && seed_width_cur[offset].bid == n_seed_mm-1
                     && seed_width_cur[offset+1].bid == n_seed_mm-1
                     && seed_width_cur[offset+1].w == seed_width_cur[offset].w) {
                      allow_mm = 0;
                  }
              }
          }

          // retrieve the next SA interval
          tmap_bwt_match_2occ4(bwt[1-strand], &e->match_sa, match_sa_next); 

          // insertions/deletions
          if(opt->indel_ends_bound <= offset && offset < len - opt->indel_ends_bound) { // do not add gaps round the ends
              if(STATE_M == e->state) { // gap open
                  if(0 < n_gapo) { // gap open is allowed
                      // insertion
                      // remember to use 'offset+1' and 'match_sa_cur' to skip over a read base
                      tmap_map1_aux_stack_push(stack, strand, offset+1, &match_sa_cur, e->n_mm, e->n_gapo + 1, e->n_gape, STATE_I, 1, e, opt);

                      // deletion
                      for(j = 0; j != 4; ++j) {
                          if(match_sa_next[j].k <= match_sa_next[j].l) {
                              //   remember that a gap deletion does not consume a
                              //   read base, so use 'offset'
                              tmap_map1_aux_stack_push(stack, strand, offset, &match_sa_next[j], e->n_mm, e->n_gapo + 1, e->n_gape, STATE_D, 1, e, opt);
                          }
                      }
                  }
              }
              else if(STATE_I == e->state) { // extension of an insertion
                  if(0 < n_gape) { // gap extension is allowed
                      // remember to use 'offset+1' and 'match_sa_cur' to skip over a read base
                      tmap_map1_aux_stack_push(stack, strand, offset+1, &match_sa_cur, e->n_mm, e->n_gapo, e->n_gape + 1, STATE_I, 1, e, opt);
                  }
              }
              else if(STATE_D == e->state) { // extension of a deletion
                  if(0 < n_gape 
                     && e->match_sa.l - e->match_sa.k + 1 < opt->max_cals_del) { // gap extension is allowed
                      for(j = 0; j != 4; ++j) {
                          if(match_sa_next[j].k <= match_sa_next[j].l) {
                              //   remember that a gap deletion does not consume a
                              //   read base, so use 'match_sa_next[j].offset-1'
                              tmap_map1_aux_stack_push(stack, strand, offset, &match_sa_next[j], e->n_mm, e->n_gapo, e->n_gape + 1, STATE_D, 1, e, opt);
                          }
                      }
                  }
              }
          }

          // mismatches
          if(1 == allow_mm && 1 == allow_diff && offset < len) { // mismatches allowed
              for(j=0;j<4;j++) {
                  int32_t c = (str[offset] + j) & 3;
                  int32_t is_mm = (0 < j || 3 < str[offset]);
                  if(match_sa_next[c].k <= match_sa_next[c].l) {
                      tmap_map1_aux_stack_push(stack, strand, offset+1, &match_sa_next[c], e->n_mm + is_mm, e->n_gapo, e->n_gape, STATE_M, is_mm, e, opt);
                  }
              }
          } 
          else if(str[offset] < 4) { // try exact match only
              int32_t c = str[offset] & 3;
              if(match_sa_next[c].k <= match_sa_next[c].l) {
                  tmap_map1_aux_stack_push(stack, strand, offset+1, &match_sa_next[c], e->n_mm, e->n_gapo, e->n_gape, STATE_M, 0, e, opt);
              }
          }
      }
  }

  return tmap_map1_sam_to_real(sams, bases, seed2_len, refseq, bwt, sa, opt);
}
