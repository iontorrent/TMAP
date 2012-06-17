/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>

#include "../../util/tmap_error.h"
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_fibheap.h"
#include "../../util/tmap_definitions.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_bwt_match.h"
#include "../../index/tmap_bwt_match_hash.h"
#include "../../index/tmap_sa.h"
#include "../../index/tmap_index.h"
#include "../../sw/tmap_sw.h"
#include "../util/tmap_map_util.h"
#include "tmap_map1.h"
#include "tmap_map1_aux.h"

#define TMAP_MAP1_AUX_STACK_INIT_SIZE 4096
#define TMAP_MAP1_AUX_STACK_TOO_BIG 262144 // TODO: are these good values?

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

#define aln_score(m,o,e,p) ((m)*((p)->pen_mm) + (o)*((p)->pen_gapo + (p)->pen_gape) + (e)*((p)->pen_gape))

#define tmap_map1_aux_reverse_query(_query, _ql, _i) \
  for(_i=0;_i<(_ql>>1);_i++) { \
      uint8_t _tmp = _query[_i]; \
      _query[_i] = _query[_ql-_i-1]; \
      _query[_ql-_i-1] = _tmp; \
  }


// Dummy tuple
typedef struct {
    tmap_bwt_int_t k, l;
} tmap_map1_aux_occ_t;

/*
static int 
tmap_map1_aux_stack_cmp(void *a, void *b)
{
  tmap_map1_aux_stack_entry_t *x = (tmap_map1_aux_stack_entry_t*)a;
  tmap_map1_aux_stack_entry_t *y = (tmap_map1_aux_stack_entry_t*)b;

  // sort by score, then offset 
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
*/

/*
static void
tmap_map1_aux_stack_entry_print(tmap_file_t *fp, tmap_map1_aux_stack_entry_t *e)
{
  fprintf(stderr, "score=%u, n_mm=%u, n_gapo=%d, n_gape=%d, state=%u, offset=%d, last_diff_offset=%d, k=%llu, l=%llu, i=%u, prev_i=%d\n",
  //tmap_file_fprintf(fp, "score=%u, n_mm=%u, n_gapo=%d, n_gape=%d, state=%u, offset=%d, last_diff_offset=%d, k=%llu, l=%llu, i=%u, prev_i=%d\n",
                    e->score, e->n_mm, e->n_gapo, e->n_gape, e->state, e->offset, e->last_diff_offset, 
                    (unsigned long long int)e->match_sa.k, (unsigned long long int)e->match_sa.l, e->i, e->prev_i);
}
*/

static void
tmap_map1_aux_stack_init_helper(tmap_map1_aux_stack_t *stack)
{
  int32_t i;

  // small memory pool
  stack->entry_pool_length = TMAP_MAP1_AUX_STACK_INIT_SIZE; 
  stack->entry_pool = tmap_calloc(stack->entry_pool_length, sizeof(tmap_map1_aux_stack_entry_t*), "stack->entry_pool");
  for(i=0;i<stack->entry_pool_length;i++) {
      stack->entry_pool[i] = tmap_calloc(1, sizeof(tmap_map1_aux_stack_entry_t), "stack->entry_pool[i]");
  }

  // nullify bins
  stack->n_bins = 0;
  stack->bins = NULL;

  // be paranoid
  stack->entry_pool_i = 0;
  stack->best_score = 0;   
  stack->n_entries = 0;
}

tmap_map1_aux_stack_t *
tmap_map1_aux_stack_init()
{
  tmap_map1_aux_stack_t *stack = NULL;
  stack = tmap_calloc(1, sizeof(tmap_map1_aux_stack_t), "stack");

  tmap_map1_aux_stack_init_helper(stack);

  return stack;
}

static void
tmap_map1_aux_stack_destroy_helper(tmap_map1_aux_stack_t *stack, int32_t full)
{
  int32_t i;
  for(i=0;i<stack->n_bins;i++) {
      free(stack->bins[i].entries);
  }
  free(stack->bins);
  for(i=0;i<stack->entry_pool_length;i++) {
      free(stack->entry_pool[i]);
  }
  free(stack->entry_pool);
  if(1 == full) {
      free(stack);
  }
}

void
tmap_map1_aux_stack_destroy(tmap_map1_aux_stack_t *stack)
{
  int32_t i;
  for(i=0;i<stack->n_bins;i++) {
      free(stack->bins[i].entries);
  }
  free(stack->bins);
  for(i=0;i<stack->entry_pool_length;i++) {
      free(stack->entry_pool[i]);
  }
  free(stack->entry_pool);
  free(stack);
}

static tmap_map1_aux_stack_t*
tmap_map1_aux_stack_reset(tmap_map1_aux_stack_t *stack,
                          int32_t max_mm, int32_t max_gapo, int32_t max_gape, 
                          const tmap_map_opt_t *opt)
{
  int32_t i;
  //int32_t i, j;
  int32_t n_bins_needed = 0;
  // move to the beginning of the memory pool
  stack->entry_pool_i = 0;
  stack->best_score = INT32_MAX;

  if(TMAP_MAP1_AUX_STACK_TOO_BIG < stack->entry_pool_length) {
      tmap_map1_aux_stack_destroy_helper(stack, 0);
      tmap_map1_aux_stack_init_helper(stack);
  }

  // clear the bins 
  for(i=0;i<stack->n_bins;i++) {
      /*
      for(j=0;j<stack->bins[i].n_entries;j++) {
          stack->bins[i].entries[j] = NULL;
      }
      */
      stack->bins[i].n_entries = 0;

  }
  // resize the bins if necessary
  n_bins_needed = aln_score(max_mm+1, max_gapo+1, max_gape+1, opt);
  if(stack->n_bins < n_bins_needed) {
      // realloc
      tmap_roundup32(n_bins_needed);
      stack->bins = tmap_realloc(stack->bins, sizeof(tmap_map1_aux_bin_t) * n_bins_needed, "stack->bins"); 
      // initialize
      for(i=stack->n_bins;i<n_bins_needed;i++) {
          stack->bins[i].n_entries = stack->bins[i].m_entries = 0;
          stack->bins[i].entries = NULL;
      }
      stack->n_bins = n_bins_needed;
  }
  stack->n_entries = 0;
  
  return stack;
}

static inline void
tmap_map1_aux_stack_push(tmap_map1_aux_stack_t *stack, 
                         int32_t offset,
                         tmap_bwt_match_occ_t *match_sa_prev,
                         int32_t n_mm, int32_t n_gapo, int32_t n_gape,
                         int32_t state, int32_t is_diff, 
                         tmap_map1_aux_stack_entry_t *prev_entry,
                         const tmap_map_opt_t *opt)
{
  int32_t i;
  int32_t n_bins_needed = 0;
  tmap_map1_aux_stack_entry_t *entry = NULL;
  tmap_map1_aux_bin_t *bin = NULL;

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
  entry->match_sa = (*match_sa_prev); 
  entry->i = stack->entry_pool_i;
  entry->offset = offset;
  if(NULL == prev_entry) {
      entry->last_diff_offset = offset;
      entry->prev_i = -1;
  }
  else {
      entry->last_diff_offset = (1 == is_diff) ? (offset) : prev_entry->last_diff_offset; 
      entry->prev_i = prev_entry->i;
  }

  if(stack->n_bins <= entry->score) {
      //tmap_bug();
      // resize the bins if necessary
      n_bins_needed = entry->score + 1;
      // realloc
      tmap_roundup32(n_bins_needed);
      stack->bins = tmap_realloc(stack->bins, sizeof(tmap_map1_aux_bin_t) * n_bins_needed, "stack->bins"); 
      // initialize
      for(i=stack->n_bins;i<n_bins_needed;i++) {
          stack->bins[i].n_entries = stack->bins[i].m_entries = 0;
          stack->bins[i].entries = NULL;
      }
      stack->n_bins = n_bins_needed;
  }
  if(stack->n_bins <= entry->score) {
      tmap_bug();
  }
  bin = &stack->bins[entry->score];
  
  // - remove duplicates
  // - most likely formed by tandem repeats or indels
  // - too computationally expensive, and not necessary
  /*
  for(i=0;i<bin->n_entries;i++) {
      if(bin->entries[i]->match_sa.k == entry->match_sa.k 
         && bin->entries[i]->match_sa.l == entry->match_sa.l 
         && bin->entries[i]->offset == entry->offset
         && bin->entries[i]->state == entry->state) {
          return;
      }
  }
  */
  
  // update best score
  if(stack->best_score > entry->score) stack->best_score = entry->score;

  if(bin->m_entries <= bin->n_entries) {
      bin->m_entries++;
      tmap_roundup32(bin->m_entries);
      bin->entries = tmap_realloc(bin->entries, sizeof(tmap_map1_aux_bin_t) * bin->m_entries, "bin->entries");
  }
  bin->entries[bin->n_entries] = entry;
  bin->n_entries++;

  stack->entry_pool_i++;
  stack->n_entries++;
}

static inline tmap_map1_aux_stack_entry_t *
tmap_map1_aux_stack_pop(tmap_map1_aux_stack_t *stack)
{
  int32_t i;
  tmap_map1_aux_bin_t *bin;
  tmap_map1_aux_stack_entry_t *best = NULL;

  if(0 == stack->n_entries) {
      return NULL;
  }
  
  // remove from the appropriate bin
  bin = &stack->bins[stack->best_score];
  if(0 == bin->n_entries) {
      tmap_bug();
  }
  best = bin->entries[bin->n_entries-1];
  bin->entries[bin->n_entries-1] = NULL;
  bin->n_entries--;
  stack->n_entries--;

  if(0 == stack->n_entries) {
      stack->best_score = INT32_MAX;
  }
  else if(0 == bin->n_entries) { // find the next best
      for(i=stack->best_score;i<stack->n_bins;i++) {
          if(0 < stack->bins[i].n_entries) {
              stack->best_score = i;
              break;
          }
      }
      if(i == stack->n_bins) {
          tmap_bug();
      }
  }

  return best;
}

static inline int32_t
tmap_map1_aux_stack_size(tmap_map1_aux_stack_t *stack)
{
  return stack->n_entries;
}

static inline void 
tmap_map1_aux_stack_shadow(tmap_bwt_int_t x, tmap_bwt_int_t max, 
                           int32_t last_diff_offset, tmap_bwt_match_width_t *w)
{
  int32_t i, j;
  for(i=j=0;i<last_diff_offset;i++) {
      /*
      fprintf(stderr, "B i=%d j=%d last_diff_offset=%d x=%u max=%u w[i].w=%u w[i].bid=%d\n",
              i, j, last_diff_offset, x, max, w[i].w, w[i].bid);
              */
      if(x < w[i].w) { 
          w[i].w -= x;
      }
      else {
          w[i].bid++;
          w[i].w = max - (++j);
      }
      /*
      fprintf(stderr, "A i=%d j=%d last_diff_offset=%d x=%u max=%u w[i].w=%u w[i].bid=%d\n",
              i, j, last_diff_offset, x, max, w[i].w, w[i].bid);
              */
  }
}

static tmap_map_sams_t *
tmap_map1_sam_to_real(tmap_map_sams_t *sams, tmap_map1_aux_occ_t *occs, tmap_string_t *bases, int32_t seed2_len,
                       tmap_refseq_t *refseq, tmap_bwt_t *bwt, tmap_sa_t *sa, tmap_bwt_match_hash_t *hash, tmap_map_opt_t *opt) 
{
  tmap_map_sams_t *sams_tmp = NULL;
  tmap_map_sam_t *sam_cur = NULL;
  uint32_t i, j, aln_ref;
  tmap_bwt_int_t k, n, num_all_sa;

  // max # of entries
  for(i=n=0;i<sams->n;i++) {
      n += occs[i].l - occs[i].k + 1;
  }
  num_all_sa = n;

  // bound the # of returned hits
  if(opt->max_best_cals < n) {
      n = opt->max_best_cals;
  }

  // alloc
  sams_tmp = tmap_map_sams_init(sams);
  tmap_map_sams_realloc(sams_tmp, n);
            
  // copy over
  for(i=j=0;i<sams->n && j<n;i++) {
      tmap_map_sam_t *sam;

      sam = &sams->sams[i];

      // go through SA interval
      for(k=occs[i].k;k<=occs[i].l;k++) {
          uint32_t pos = 0, seqid = 0;
          tmap_bwt_int_t pacpos = 0;
          uint8_t strand;

          sam_cur = &sams_tmp->sams[j];
          tmap_map_sam_init(sam_cur);

          strand = sams->sams[i].strand;
          aln_ref = sam->aux.map1_aux->aln_ref;
              
          pacpos = bwt->seq_len - tmap_sa_pac_pos_hash(sa, bwt, k, hash);
          pacpos = (pacpos < aln_ref) ? 0 : (pacpos - aln_ref);
          pacpos++; // make one-based
          
          // NB: addressing the symptom, not the problem
          // This happens when we are at the end of the reference on the reverse
          // strand
          if(pacpos <= refseq->len && refseq->len < pacpos + aln_ref - 1) {
              aln_ref = refseq->len - pacpos + 1;
          }
          else if(refseq->len * 2 < pacpos + aln_ref - 1) {
              aln_ref = (refseq->len * 2) - pacpos + 1;
          }
          
          // save the hit
          if(0 < tmap_refseq_pac2real(refseq, pacpos, aln_ref, &seqid, &pos, &strand)) {
              // copy over previous parameters
              sam_cur->algo_id = TMAP_MAP_ALGO_MAP1;
              sam_cur->algo_stage = opt->algo_stage;
              sam_cur->strand = strand;
              sam_cur->seqid = seqid;
              sam_cur->pos = pos-1; // adjust to zero-based
              sam_cur->target_len = aln_ref;
              sam_cur->score_subo = INT32_MIN;
              if(0 < opt->seed2_length && seed2_len < bases->l) { // adjust if we used a secondary seed
                  // adjust both the target length and position
                  if(0 == strand) { // forward
                      sam_cur->target_len += (bases->l - seed2_len);
                      // NB: sam_cur->pos is zero based
                      if(refseq->annos[sam_cur->seqid].len < sam_cur->pos + sam_cur->target_len) {
                          sam_cur->target_len = refseq->annos[sam_cur->seqid].len - sam_cur->pos;
                      }
                  }
                  else { // reverse
                      if(sam_cur->pos < (bases->l - seed2_len)) { // before the start of the chromosome
                          sam_cur->target_len += sam_cur->pos;
                          sam_cur->pos = 0;
                      }
                      else { // move to the end of the read
                          sam_cur->target_len += (bases->l - seed2_len);
                          sam_cur->pos -= (bases->l - seed2_len);
                      }
                  }
              }
              else if(sam_cur->target_len < bases->l) { // do not adjust, we used the full read
                  sam_cur->target_len = bases->l;
              }

              // aux
              tmap_map_sam_malloc_aux(sam_cur);
              sam_cur->aux.map1_aux->n_mm = sam->aux.map1_aux->n_mm;
              sam_cur->aux.map1_aux->n_gapo = sam->aux.map1_aux->n_gapo;
              sam_cur->aux.map1_aux->n_gape = sam->aux.map1_aux->n_gape;
              sam_cur->aux.map1_aux->aln_ref = 0;
              sam_cur->aux.map1_aux->num_all_sa = num_all_sa;
              j++;

              // only save the top n hits
              if(n <= j) {  // equivalent to opt->max_bset_cals <= j
                  break;
              }
          }
      }
  }
  
  // destroy
  tmap_map_sams_destroy(sams);

  // realloc
  tmap_map_sams_realloc(sams_tmp, j);

  //NB: do not use occs later
  free(occs);
  occs = NULL;
  
  return sams_tmp;
}

tmap_map_sams_t *
tmap_map1_aux_core(tmap_seq_t *seq, tmap_index_t *index,
                   tmap_bwt_match_hash_t *hash,
                   tmap_bwt_match_width_t *width, tmap_bwt_match_width_t *seed_width, tmap_map_opt_t *opt,
                   tmap_map1_aux_stack_t *stack, int32_t seed2_len)
{
  int32_t max_mm = opt->max_mm, max_gapo = opt->max_gapo, max_gape = opt->max_gape, seed_max_diff = opt->seed_max_diff;
  int32_t best_score, next_best_score;
  int32_t best_cnt = 0;
  int32_t i, j, num_n = 0;
  int32_t max_edit_score;
  tmap_bwt_match_occ_t match_sa_start;
  tmap_string_t *bases=NULL;
  tmap_map_sams_t *sams = NULL;
  int32_t max_diff, best_diff;
  tmap_bwt_int_t k, l;
  tmap_refseq_t *refseq = index->refseq;
  tmap_bwt_t *bwt = index->bwt;
  tmap_sa_t *sa = index->sa;
  tmap_map1_aux_occ_t *occs = NULL;


  max_edit_score = opt->pen_mm;
  //if(max_edit_score < opt->pen_gapo + opt->pen_gape) max_edit_score = opt->pen_gapo + opt->pen_gape;
  //if(max_edit_score < opt->pen_gape) max_edit_score = opt->pen_gape;

  bases = tmap_seq_get_bases(seq);
  /*
  fputc('\n', stderr);
  for(i=0;i<bases->l;i++) {
      fputc("ACGTN"[(int)bases->s[i]], stderr);
  }
  fputc('\n', stderr);
  */
  
  // the maximum # of differences
  if(bases->l <= TMAP_MAP_OPT_MAX_DIFF_READ_LENGTH) {
      best_diff = max_diff = opt->max_diff_table[bases->l];
  }
  else {
      best_diff = max_diff = opt->max_diff_table[TMAP_MAP_OPT_MAX_DIFF_READ_LENGTH];
  }
  
  // bound differenes by the maximum # of differences
  if(max_diff < max_mm) max_mm = max_diff;
  if(max_diff < max_gapo) max_gapo = max_diff;
  //if(max_diff < max_gape) max_gape = max_diff;
  
  best_score = next_best_score = aln_score(max_mm+1, max_gapo+1, max_gape+1, opt);

  // check whether there are too many N
  for(j=bases->l-seed2_len,num_n=0;j<bases->l;j++) {
      if(3 < bases->s[j]) {
          num_n++;
      }
  }
  if(max_mm < num_n || max_diff < num_n) {
      return tmap_map_sams_init(NULL);
  }

  // initialize
  sams = tmap_map_sams_init(NULL);
  occs = NULL;

  match_sa_start.offset = 0;
  match_sa_start.hi = 0;
  match_sa_start.k = 0;
  match_sa_start.l = bwt->seq_len;

  stack = tmap_map1_aux_stack_reset(stack, max_mm, max_gapo, max_gape, opt); // reset stack
  tmap_map1_aux_stack_push(stack, bases->l, &match_sa_start, 0, 0, 0, STATE_M, 0, NULL, opt);

  while(0 < tmap_map1_aux_stack_size(stack) && tmap_map1_aux_stack_size(stack) < opt->max_entries) {
      tmap_map1_aux_stack_entry_t *e = NULL;
      int32_t len=-1; 
      int32_t n_seed_mm=0, offset, width_cur_i;
      const uint8_t *str=NULL;
      int32_t sam_found, m;
      tmap_bwt_match_width_t *width_cur = NULL;
      const tmap_bwt_match_width_t *seed_width_cur = NULL;
      tmap_bwt_match_occ_t match_sa_cur, match_sa_next[4];
      
      // get the best entry
      e = tmap_map1_aux_stack_pop(stack); 

      // bound with best score
      if(best_score + max_edit_score < e->score) {
          break; // no need to continue
      }

      // some more information
      match_sa_cur = e->match_sa; 

      // check if we have too many edits
      m = max_diff - (e->n_mm + e->n_gapo + e->n_gape);
      if(m < 0) {
          continue; // too many edits
      }

      // get the rest of the information
      offset = e->offset; // zero-based
      str = (uint8_t*)bases->s;
      len = bases->l;
      width_cur = width;
      width_cur_i = seed2_len - (len - offset);

      if(NULL != seed_width) {
          seed_width_cur = seed_width;
          n_seed_mm = seed_max_diff - (e->n_mm + e->n_gapo + e->n_gape); // consider only mismatches in the seed
      }
      else {
          seed_width_cur = NULL;
      }
      if(0 < width_cur_i && m < width_cur[width_cur_i-1].bid) { // too many edits
          continue;
      }

      // check whether a sam is found
      sam_found = 0;
      if(len - seed2_len == offset) {
          sam_found = 1;
      }
      else if(max_mm == e->n_mm // no mismatches from any state
              && ((e->state == STATE_M && max_gapo == e->n_gapo) // in STATE_M but no more gap opens
                  || (e->state != STATE_M && max_gape == e->n_gape))) { // in STATE_I/STATE_D but no more extensions
          if(0 < tmap_bwt_match_hash_exact_alt_reverse(bwt, offset, str, &match_sa_cur, hash)) { // the alignment must match exactly to sam
              sam_found = 2;
          }
          else {
              continue; // no sam, skip
          }
      }

      if(0 < sam_found) { // alignment found
          // check for duplicates
          if(0 < sams->n) {
              for(i=0;i<sams->n;i++) {
                  // check contained
                  if(match_sa_cur.k <= occs[i].k
                     && occs[i].k <= match_sa_cur.l) { // MK <= SK <= ML
                      if(occs[i].l <= match_sa_cur.l) { // MK <= SK <= SL <= ML
                          // Want (SK - MK) + (ML - SL)
                          k = occs[i].k - match_sa_cur.k; // (SK - MK)
                          k += match_sa_cur.l - occs[i].l; // (ML - SL)
                          occs[i].l = match_sa_cur.l; // Make SL = ML
                      }
                      else { // MK <= SK <= ML <= SL
                          k = occs[i].k - match_sa_cur.k; // (SK - MK)
                      }
                      occs[i].k = match_sa_cur.k; // Make SK = MK
                      break;
                  }
                  else if(match_sa_cur.k <= occs[i].l
                          && occs[i].l <= match_sa_cur.l) { // MK <= SL <= ML
                      if(match_sa_cur.k <= occs[i].k) { // MK <= SK <= SL <= ML
                          // Want (SK - MK) + (ML - SL)
                          k = occs[i].k - match_sa_cur.k; // (SK - MK)
                          k += match_sa_cur.l - occs[i].l; // (ML - SL)
                          occs[i].k = match_sa_cur.k; // Make SK = MK
                      }
                      else { // SK <= MK <= SL <= ML
                          k = match_sa_cur.l - occs[i].l; // (ML - SL)
                      }
                      occs[i].l = match_sa_cur.l; // Make SL = ML
                      break;
                  }
              }
              if(i < sams->n) {
                  // shadow
                  if(0 < k) {
                      //tmap_map1_aux_stack_shadow(k, bwt->seq_len, e->last_diff_offset, width_cur);
                      width_cur_i = seed2_len - (len - e->last_diff_offset);
                      tmap_map1_aux_stack_shadow(k, seed2_len, width_cur_i, width_cur);
                  }
                  sam_found = 0;
                  continue;
              }
          }

          int32_t score = aln_score(e->n_mm, e->n_gapo, e->n_gape, opt);
          int32_t do_add = 1;
          if(sams->n == 0) {
              best_score = score;
              best_cnt = 0;
              best_diff = e->n_mm + e->n_gapo + e->n_gape;
          }
          if(score == best_score) {
              best_cnt += match_sa_cur.l - match_sa_cur.k + 1;
          }
          else {
              if(best_diff + 1 <= max_diff) {
                  max_diff = best_diff + 1;
              }
              if(score < next_best_score) {
                  next_best_score = score;
              }
              else if(next_best_score < score) {
                  // no need to examine further
                  break;
              }
          }
          if(do_add) { // append
              uint32_t op, op_len, cigar_i;
              tmap_map_sam_t *sam = NULL;
              tmap_map1_aux_stack_entry_t *cur = NULL;
  
              tmap_map_sams_realloc(sams, sams->n+1);
              occs = tmap_realloc(occs, sizeof(tmap_map1_aux_occ_t) * sams->n, "occs");

              sam = &sams->sams[sams->n-1];

              sam->algo_id = TMAP_MAP_ALGO_MAP1;
              sam->algo_stage = 0;
              sam->score = e->score;

              // aux data
              tmap_map_sam_malloc_aux(sam);
              k = occs[sams->n-1].k = match_sa_cur.k;
              l = occs[sams->n-1].l= match_sa_cur.l;
              sam->aux.map1_aux->n_mm = e->n_mm;
              sam->aux.map1_aux->n_gapo = e->n_gapo;
              sam->aux.map1_aux->n_gape = e->n_gape;

              // aux data: reference length
              cur = e;
              i = e->i;
              sam->aux.map1_aux->aln_ref = 0;
              cigar_i = 0;
              if(2 == sam_found) { // we used 'tmap_bwt_match_exact_alt_reverse' 
                  op = STATE_M;
                  op_len = offset;
              }
              else {
                  op = -1;
                  op_len = 0;
              }
              while(0 <= i) {
                  cur = stack->entry_pool[i];
                  if(len == cur->offset) break;
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
                  //fprintf(stderr, "cur->state=%c op_len=%d cur->prev_i=%d k=%u l=%u\n", "MIDS"[cur->state], op_len, cur->prev_i, cur->match_sa.k, cur->match_sa.l);
                  i = cur->prev_i;
              }
              if(STATE_M == op || STATE_D == op) {
                  sam->aux.map1_aux->aln_ref += op_len;
              }

              /*
              fprintf(stderr, "shadow 2 k=%u l=%u len=%d offset=%d last_diff_offset=%d\n",
                      k, l, len, offset, e->last_diff_offset);
              fprintf(stderr, "e->n_mm=%d e->n_gapo=%d e->n_gape=%d\n",
                      e->n_mm, e->n_gapo, e->n_gape);
              */
              //tmap_map1_aux_stack_shadow(l - k + 1, bwt->seq_len, e->last_diff_offset, width_cur);
              width_cur_i = seed2_len - (len - e->last_diff_offset);
              tmap_map1_aux_stack_shadow(l - k + 1, seed2_len, width_cur_i, width_cur);
              if(opt->max_best_cals < best_cnt) {
                  // ignore if too many "best" have been found
                  occs[sams->n-1].l -= (best_cnt - opt->max_best_cals); // only save the maximum
                  break;
              }
          }
      }
      else {
          int32_t allow_diff = 1, allow_mm = (e->n_mm < max_mm) ? 1 : 0;

          // decrement the offset
          offset--;

          // use a bound for mismatches
          if(0 < offset) {
              int32_t seed_width_cur_i = offset - (len - opt->seed_length);
              width_cur_i = seed2_len - (len - offset);
              if(0 < width_cur_i) {
                  if(m-1 < width_cur[width_cur_i-1].bid) { 
                      allow_diff = 0;
                  }
                  else if(width_cur[width_cur_i-1].bid == m-1
                          && width_cur[width_cur_i].bid == m-1
                          && width_cur[width_cur_i-1].w == width_cur[width_cur_i].w) {
                      allow_mm = 0;
                  }
              }
              if(0 < seed_width_cur_i) {
                  if(NULL != seed_width_cur && 0 < seed_width_cur_i) {
                      if(n_seed_mm-1 < seed_width_cur[seed_width_cur_i-1].bid) {
                          allow_diff = 0;
                      }
                      else if(seed_width_cur[seed_width_cur_i-1].bid == n_seed_mm-1
                              && seed_width_cur[seed_width_cur_i].bid == n_seed_mm-1
                              && seed_width_cur[seed_width_cur_i-1].w == seed_width_cur[seed_width_cur_i].w) {
                          allow_mm = 0;
                      }
                  }
              }
          }

          // retrieve the next SA interval
          tmap_bwt_match_hash_2occ4(bwt, &e->match_sa, match_sa_next, hash); 

          // insertions/deletions
          if(allow_diff 
             && opt->indel_ends_bound + e->n_gapo + e->n_gape <= offset
             && opt->indel_ends_bound + e->n_gapo + e->n_gape <= len - offset) { // check to add gaps
              if(STATE_M == e->state) { // gap open
                  if(e->n_gapo < max_gapo) { // gap open is allowed
                      // insertion
                      tmap_map1_aux_stack_push(stack, offset, &match_sa_cur, e->n_mm, e->n_gapo + 1, e->n_gape, STATE_I, 1, e, opt);

                      // deletion
                      for(j = 0; j != 4; ++j) {
                          if(match_sa_next[j].k <= match_sa_next[j].l) {
                              //   remember that a gap deletion does not consume a
                              //   read base, so use 'offset+1'
                              tmap_map1_aux_stack_push(stack, offset+1, &match_sa_next[j], e->n_mm, e->n_gapo + 1, e->n_gape, STATE_D, 1, e, opt);
                          }
                      }
                  }
              }
              else if(STATE_I == e->state) { // extension of an insertion
                  if(e->n_gape < max_gape) { // gap extension is allowed
                      tmap_map1_aux_stack_push(stack, offset, &match_sa_cur, e->n_mm, e->n_gapo, e->n_gape + 1, STATE_I, 1, e, opt);
                  }
              }
              else if(STATE_D == e->state) { // extension of a deletion
                  if(e->n_gape < max_gape) {
                      if(e->n_gape + e->n_gapo < max_diff || e->match_sa.l - e->match_sa.k + 1 < opt->max_cals_del) { // gap extension is allowed
                          for(j = 0; j != 4; ++j) {
                              if(match_sa_next[j].k <= match_sa_next[j].l) {
                                  //   remember that a gap deletion does not consume a
                                  //   read base, so use 'offset+1'
                                  tmap_map1_aux_stack_push(stack, offset+1, &match_sa_next[j], e->n_mm, e->n_gapo, e->n_gape + 1, STATE_D, 1, e, opt);
                              }
                          }
                      }
                  }
              }
          }

          // mismatches
          if(1 == allow_mm && 1 == allow_diff) { // mismatches allowed
              for(j=0;j<4;j++) {
                  int32_t c = (str[offset] + j) & 3;
                  int32_t is_mm = (0 < j || 3 < str[offset]);
                  if(match_sa_next[c].k <= match_sa_next[c].l) {
                      tmap_map1_aux_stack_push(stack, offset, &match_sa_next[c], e->n_mm + is_mm, e->n_gapo, e->n_gape, STATE_M, is_mm, e, opt);
                  }
              }
          } 
          else if(str[offset] < 4) { // try exact match only
              int32_t c = str[offset] & 3;
              if(match_sa_next[c].k <= match_sa_next[c].l) {
                  tmap_map1_aux_stack_push(stack, offset, &match_sa_next[c], e->n_mm, e->n_gapo, e->n_gape, STATE_M, 0, e, opt);
              }
          }
      }
  }

  return tmap_map1_sam_to_real(sams, occs, bases, seed2_len, refseq, bwt, sa, hash, opt);
}
