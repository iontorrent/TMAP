#include <stdlib.h>

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_seq.h"
#include "../util/fmap_fibheap.h"
#include "../util/fmap_definitions.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_bwt_match.h"
#include "../index/fmap_sa.h"
#include "fmap_map1.h"
#include "fmap_map1_aux.h"

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

// TODO redefine
#define aln_score(m,o,e,p) ((m)*(p)->pen_mm + (o)*(p)->pen_gapo + (e)*(p)->pen_gape)

fmap_map1_aux_stack_t *
fmap_map1_aux_stack_init()
{
  int32_t i;
  fmap_map1_aux_stack_t *stack = NULL;
  stack = fmap_calloc(1, sizeof(fmap_map1_aux_stack_t), "stack");
  stack->entry_pool_length = 1024; // TODO: move to a define block
  stack->entry_pool = fmap_malloc(stack->entry_pool_length*sizeof(fmap_map1_aux_stack_entry_t*), "stack->entry_pool");
  for(i=0;i<stack->entry_pool_length;i++) {
      stack->entry_pool[i] = fmap_malloc(sizeof(fmap_map1_aux_stack_entry_t), "stack->entry_pool[i]");
  }
  stack->heap = fmap_fibheap_makekeyheap();

  return stack;
}

void
fmap_map1_aux_stack_destroy(fmap_map1_aux_stack_t *stack)
{
  int32_t i;
  fmap_fibheap_deleteheap(stack->heap);
  for(i=0;i<stack->entry_pool_length;i++) {
      free(stack->entry_pool[i]);
  }
  free(stack->entry_pool);
  free(stack);
}

static void
fmap_map1_aux_stack_reset(fmap_map1_aux_stack_t *stack)
{
  // reset heap
  fmap_fibheap_deleteheap(stack->heap);
  stack->heap = fmap_fibheap_makekeyheap();
  // move to the beginning of the memory pool
  stack->entry_pool_i = 0;
  stack->best_score = INT32_MAX;
}

static inline void
fmap_map1_aux_stack_push(fmap_map1_aux_stack_t *stack, int32_t strand, 
                         fmap_bwt_match_occ_t *match_sa_prev,
                         int32_t n_mm, int32_t n_gapo, int32_t n_gape,
                         int32_t state, int32_t is_diff, 
                         fmap_map1_aux_stack_entry_t *prev_entry,
                         const fmap_map1_opt_t *opt)
{
  fmap_map1_aux_stack_entry_t *entry = NULL;

  // check to see if we need more memory
  if(stack->entry_pool_length <= stack->entry_pool_i) { 
      int32_t i = stack->entry_pool_length;
      stack->entry_pool_length <<= 2;
      stack->entry_pool = fmap_realloc(stack->entry_pool, 
                                       sizeof(fmap_map1_aux_stack_entry_t*)*stack->entry_pool_length, 
                                       "stack->entry_pool");
      while(i<stack->entry_pool_length) {
          stack->entry_pool[i] = fmap_malloc(sizeof(fmap_map1_aux_stack_entry_t), "stack->entry_pool[i]");
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
  if(NULL == prev_entry) {
      entry->last_diff_offset = (1 == is_diff) ? entry->match_sa.offset : 0;
      entry->prev_i = -1;
  }
  else {
      entry->last_diff_offset = (1 == is_diff) ? entry->match_sa.offset : prev_entry->last_diff_offset;
      entry->prev_i = prev_entry->i;
  }

  if(stack->best_score > entry->score) stack->best_score = entry->score;

  fmap_fibheap_insertkey(stack->heap, entry->score, entry);

  stack->entry_pool_i++;
}

static inline fmap_map1_aux_stack_entry_t *
fmap_map1_aux_stack_pop(fmap_map1_aux_stack_t *stack)
{
  fmap_map1_aux_stack_entry_t *best, *next_best;

  best = fmap_fibheap_extractmin(stack->heap);
  next_best = fmap_fibheap_min(stack->heap);
  if(NULL != next_best) {
      stack->best_score  = next_best->score;
  }
  else {
      stack->best_score = INT32_MAX;
  }

  return best;
}

static inline void 
fmap_map1_aux_stack_shadow(int32_t x, int32_t len, uint32_t max, 
                           int32_t last_diff_offset, fmap_bwt_match_width_t *w)
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
          fmap_error("else happened", Exit, OutOfRange);
      }
  }
}

static inline int32_t
fmap_map1_aux_stack_size(fmap_map1_aux_stack_t *stack)
{
  return stack->heap->fmap_fibheap_n;
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
fmap_map1_aux_get_bam_state(int state)
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
        if(aln->cigar_length <= cigar_i) { \
            aln->cigar_length++; \
            aln->cigar = fmap_realloc(aln->cigar, sizeof(uint32_t)*aln->cigar_length, "aln->cigar"); \
        } \
        aln->cigar[cigar_i] = (op_len << 4 ) | fmap_map1_aux_get_bam_state(op); \
        cigar_i++; \
    } \
} while(0)

fmap_map1_aln_t **
fmap_map1_aux_core(fmap_seq_t *seq[2], fmap_bwt_t *bwt,
                   fmap_bwt_match_width_t *width[2], fmap_bwt_match_width_t *seed_width[2], fmap_map1_opt_t *opt,
                   fmap_map1_aux_stack_t *stack)
{
  int32_t max_mm = opt->max_mm, max_gapo = opt->max_gapo, max_gape = opt->max_gape;
  int32_t best_score = aln_score(max_mm+1, max_gapo+1, max_gape+1, opt);
  int32_t best_cnt = 0;
  int32_t j, _j;
  int32_t min_edit_score;
  int32_t alns_num = 0; 
  fmap_bwt_match_occ_t match_sa_start;
  fmap_map1_aln_t **alns=NULL;

  // HERE
  int32_t hw = bwt->hash_width;
  bwt->hash_width = 0;

  if(0 == bwt->is_rev) fmap_error("0 == bwt->is_rev", Exit, OutOfRange);

  min_edit_score = opt->pen_mm;
  if(opt->pen_gapo < min_edit_score) min_edit_score = opt->pen_gapo;
  if(opt->pen_gape < min_edit_score) min_edit_score = opt->pen_gape;

  // check whether there are too many N
  for(j=_j=0;j<seq[0]->seq.l;j++) {
      if(3 < seq[0]->seq.s[j]) _j++;
  }
  if(max_mm < _j) {
      return NULL;
  }

  match_sa_start.offset = 0;
  match_sa_start.hi = 0;
  match_sa_start.k = 0;
  match_sa_start.l = bwt->seq_len;

  fmap_map1_aux_stack_reset(stack); // reset stack
  fmap_map1_aux_stack_push(stack, 0, &match_sa_start, 0, 0, 0, STATE_M, 0, NULL, opt);
  fmap_map1_aux_stack_push(stack, 1, &match_sa_start, 0, 0, 0, STATE_M, 0, NULL, opt);

  while(0 < fmap_map1_aux_stack_size(stack) && fmap_map1_aux_stack_size(stack) < opt->max_entries) {
      fmap_map1_aux_stack_entry_t *e = NULL;
      int32_t strand, len; 
      int32_t n_mm, n_gapo, n_gape;
      int32_t n_seed_mm;
      const uint8_t *str=NULL;
      int32_t hit_found;
      fmap_bwt_match_width_t *width_cur;
      const fmap_bwt_match_width_t *seed_width_cur = NULL;
      fmap_bwt_match_occ_t match_sa_cur, match_sa_tmp, match_sa_next[4];
      // int32_t i, m, m_seed = 0;

      e = fmap_map1_aux_stack_pop(stack); // get the best entry
      if(best_score + min_edit_score < e->score) break; // no need to continue
                  
      strand = e->strand; // strand;
      match_sa_cur = e->match_sa; 
      n_mm = opt->max_mm - e->n_mm;
      n_gapo = opt->max_gapo - e->n_gapo;
      n_gape = opt->max_gape - e->n_gape;
      if(n_mm < 0 || n_gapo < 0 || n_gape < 0) continue; // too many edits

      str = (uint8_t*)seq[strand]->seq.s; 
      len = seq[strand]->seq.l;
      width_cur = width[strand];
      seed_width_cur = (NULL == seed_width) ? NULL : seed_width[strand];

      /*
      fprintf(stderr, "strand=%d offset=%d k=%u l=%u hi=%u e->n_mm=%d e->n_gapo=%d e->n_gape=%d score=%d\n",
              strand, match_sa_cur.offset, match_sa_cur.k, match_sa_cur.l, match_sa_cur.hi, e->n_mm, e->n_gapo, e->n_gape, e->score);
              */

      if(NULL != seed_width_cur) { // apply seeding
          seed_width_cur = seed_width[strand];
          n_seed_mm = opt->seed_max_mm - e->n_mm; // ignore gaps TODO: is this correct?
      }
      // if there are no more gaps, check if we are allowed mismatches
      if(match_sa_cur.offset+1 < len && 0 == n_gapo && 0 == n_gape && n_mm < width_cur[match_sa_cur.offset].bid) {
          continue;
      } 

      // check whether a hit is found
      hit_found = 0;
      if(len == match_sa_cur.offset) {
          hit_found = 1;
      }
      else if(0 == n_mm // no mismatches from any state
              && (e->state == STATE_M && 0 == n_gapo) // in STATE_M but no more gap opens
              && (e->state != STATE_M && 0 == n_gape)) { // in STATE_I/STATE_D but no more extensions
          if(0 < fmap_bwt_match_exact_alt(bwt, match_sa_cur.offset, str, &match_sa_cur.k, &match_sa_cur.l)) { // the alignment must match exactly to hit
              hit_found = 1;
          }
          else {
              continue; // no hit, skip
          }
      }
      
      if(1 == hit_found) { // alignment found
          int32_t score = aln_score(e->n_mm, e->n_gapo, e->n_gape, opt);
          int32_t do_add = 1;
          //printf("#2 hits found: %d:(%u,%u)\n", e->n_mm+e->n_gapo, k, l);
          if(alns_num == 0) {
              best_score = score;
              best_cnt = 0;
          }
          if(score == best_score) {
              best_cnt += match_sa_cur.l - match_sa_cur.k + 1;
          }
          else if(opt->max_best_cals < best_cnt) {
              // ignore if too many "best" have been found
              break;
          }
          if(0 < e->n_gapo) { // check if the same alignment has been found 
              for(j=0;j<alns_num;j++) {
                  if(alns[j]->k == match_sa_cur.k && alns[j]->l) {
                      if(score < alns[j]->score) { // bug!
                          fmap_error("bug encountered", Exit, OutOfRange);
                      }
                      do_add = 0;
                      break;
                  }
              }
          }
          if(do_add) { // append
              uint32_t op, op_len, cigar_i;
              fmap_map1_aln_t *aln = NULL;
              fmap_map1_aux_stack_entry_t *cur = e;

              alns_num++;
              alns = fmap_realloc(alns, sizeof(fmap_map1_aln_t*)*alns_num, "alns");
              alns[alns_num-1] = fmap_calloc(1, sizeof(fmap_map1_aln_t), "alns[alns_num-1]");
              aln = alns[alns_num-1];

              aln->score = cur->score;
              aln->n_mm = cur->n_mm;
              aln->n_gapo = cur->n_gapo;
              aln->n_gape = cur->n_gape;
              aln->strand = cur->strand;
              aln->k = cur->match_sa.k;
              aln->l = cur->match_sa.l;

              // cigar
              aln->cigar_length = 1 + cur->n_gapo;
              aln->cigar = fmap_malloc(sizeof(uint32_t)*aln->cigar_length, "aln->cigar");
              cigar_i = 0;

              if(match_sa_cur.offset < len) { // we used 'fmap_bwt_match_exact_alt' 
                  op = STATE_M;
                  op_len = len - match_sa_cur.offset;
              }
              else {
                  op = -1;
                  op_len = 0;
              }

              while(0 < cur->match_sa.offset) {
                  if(op != cur->state) {
                      __add_to_cigar();
                      op = cur->state;
                      op_len = 1;
                  }
                  else {
                      op_len++;
                  }
                  cur = stack->entry_pool[cur->prev_i];
              }
              __add_to_cigar();

              if(cigar_i < aln->cigar_length) { // reallocate to fit
                  aln->cigar_length = cigar_i;
                  aln->cigar = fmap_realloc(aln->cigar, sizeof(uint32_t)*aln->cigar_length, "aln->cigar");
              }

              // reverse the cigar 
              for(cigar_i=0;cigar_i<(aln->cigar_length >> 1);cigar_i++) {
                  op = aln->cigar[aln->cigar_length-1-cigar_i];
                  aln->cigar[aln->cigar_length-1-cigar_i] = aln->cigar[cigar_i];
                  aln->cigar[cigar_i] = op;
              }

              // print
              /*
                 for(cigar_i=0;cigar_i<aln->cigar_length;cigar_i++) {
                 fprintf(stderr, "%d%c", (aln->cigar[cigar_i] >> 4), "MIDNSHP"[aln->cigar[cigar_i] & 0xf]);
                 }
                 fputc('\n', stderr);
                 */

              // TODO: use the shadow ?
              //fmap_map1_aux_stack_shadow(l - k + 1, len, bwt->seq_len, e->last_diff_offset, width_cur);
          }
          continue; // could use an else statement
      }

      // retrieve the next SA interval
      fmap_bwt_match_2occ4(bwt, &e->match_sa, match_sa_next); 
      for(j=0;j<4;j++) {
          match_sa_next[j].k += bwt->L2[j] + 1; 
          match_sa_next[j].l += bwt->L2[j]; 
          /*
          fprintf(stderr, "e->match_sa.k=%u e->match_sa.l=%u match_sa_next[%d].k=%u match_sa_next[%d].l=%u\n",
                  e->match_sa.k, e->match_sa.l, j, match_sa_next[j].k, j, match_sa_next[j].l);
                  */
      }

      // insertions/deletions
      if(opt->indel_ends_bound <= match_sa_cur.offset && match_sa_cur.offset <= len - opt->indel_ends_bound) { // do not add gaps round the ends
          if(STATE_M == e->state) { // gap open
              if(0 < n_gapo) { // gap open is allowed
                  // insertion
                  match_sa_tmp = match_sa_cur;
                  match_sa_tmp.offset++; // skip over a base
                  fmap_map1_aux_stack_push(stack, strand, &match_sa_tmp, e->n_mm, e->n_gapo + 1, e->n_gape, STATE_I, 1, e, opt);
                  
                  // deletion
                  for(j = 0; j != 4; ++j) {
                      if(match_sa_next[j].k <= match_sa_next[j].l) {
                          match_sa_tmp = match_sa_next[j];
                          //   remember that a gap deletion does not consume a
                          //   read base, so use 'match_sa_next[j].offset-1'
                          match_sa_tmp.offset--; 
                          fmap_map1_aux_stack_push(stack, strand, &match_sa_tmp, e->n_mm, e->n_gapo + 1, e->n_gape, STATE_D, 1, e, opt);
                      }
                  }
              }
          }
          else if(STATE_I == e->state) { // extension of an insertion
              if(0 < n_gape) { // gap extension is allowed
                  match_sa_tmp = match_sa_cur;
                  match_sa_tmp.offset++; // skip over a base
                  fmap_map1_aux_stack_push(stack, strand, &match_sa_tmp, e->n_mm, e->n_gapo, e->n_gape + 1, STATE_I, 1, e, opt);
              }
          }
          else if(STATE_D == e->state) { // extension of a deletion
              if(0 < n_gape 
                 && e->match_sa.l - e->match_sa.k + 1 < opt->max_cals_del) { // gap extension is allowed
                  for(j = 0; j != 4; ++j) {
                      if(match_sa_next[j].k <= match_sa_next[j].l) {
                          match_sa_tmp = match_sa_next[j];
                          //   remember that a gap deletion does not consume a
                          //   read base, so use 'match_sa_next[j].offset-1'
                          match_sa_tmp.offset--; 
                          fmap_map1_aux_stack_push(stack, strand, &match_sa_tmp, e->n_mm, e->n_gapo, e->n_gape + 1, STATE_D, 1, e, opt);
                      }
                  }
              }
          }
      }
      // mismatches
      if(match_sa_cur.offset < len && 0 < n_mm && width_cur[match_sa_cur.offset].bid <= n_mm) { // mismatches allowed
          for(j=0;j<4;j++) {
              int32_t c = (str[match_sa_cur.offset] + j) & 3;
              int32_t is_mm = (0 < j || 3 < str[match_sa_cur.offset]);
              if(match_sa_next[c].k <= match_sa_next[c].l) {
                  fmap_map1_aux_stack_push(stack, strand, &match_sa_next[c], e->n_mm + is_mm, e->n_gapo, e->n_gape, STATE_M, is_mm, e, opt);
              }
          }
      } 
      else if(str[match_sa_cur.offset] < 4) { // try exact match only
          int32_t c = str[match_sa_cur.offset] & 3;
          if(match_sa_next[c].k <= match_sa_next[c].l) {
              fmap_map1_aux_stack_push(stack, strand, &match_sa_next[c], e->n_mm, e->n_gapo, e->n_gape, STATE_M, 0, e, opt);
          }
      }
  }

  if(0 < alns_num) { // add trailing null
      alns = fmap_realloc(alns, sizeof(fmap_map1_aln_t*)*(1+alns_num), "alns");
      alns[alns_num] = NULL;
  }

  // HERE
  bwt->hash_width = hw;

  return alns;
}
