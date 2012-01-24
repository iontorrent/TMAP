/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/resource.h>
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_sort.h"
#include "../../util/tmap_vec.h"
#include "../../util/tmap_hash.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_bwtl.h"
#include "../../index/tmap_bwt_match.h"
#include "../../index/tmap_bwt_match_hash.h"
#include "../../index/tmap_sa.h"
#include "../../index/tmap_index.h"
#include "../../sw/tmap_sw.h"
#include "../util/tmap_map_util.h"
#include "tmap_map2.h"
#include "tmap_map2_aux.h"
#include "tmap_map2_core.h"

typedef struct {
    tmap_bwt_int_t k, l;
} tmap_map2_qintv_t;

#define qintv_eq(a, b) ((a).k == (b).k && (a).l == (b).l)
#define qintv_hash(a) ((a).k>>7^(a).l<<17)

// for hashing
TMAP_HASH_INIT(tmap_map2_qintv, tmap_map2_qintv_t, uint64_t, 1, qintv_hash, qintv_eq)
TMAP_HASH_MAP_INIT_INT64(64, uint64_t)

// for sorting generically
TMAP_SORT_INIT_GENERIC(int32_t)

// TODO: document
static const tmap_map2_cell_t tmap_map2_core_default_cell = { {0, 0, 0, 0}, TMAP_MAP2_MINUS_INF, TMAP_MAP2_MINUS_INF,
     TMAP_MAP2_MINUS_INF, 0, 0, 0, -1, -1, {-1, -1, -1, -1} };

/* --- BEGIN: utilities --- */
static tmap_hash_t(64) *
tmap_map2_core_connectivity(const tmap_bwtl_t *b)
{
  tmap_hash_t(64) *h;
  uint32_t k, l, cntk[4], cntl[4];
  uint64_t x;
  tmap_hash_iter_t iter;
  int32_t j, ret;
  tmap_vec_t(uint64_t) stack;

  tmap_vec_init(stack);
  h = tmap_hash_init(64);
  tmap_hash_resize(64, h, b->seq_len * 4);
  x = b->seq_len;
  tmap_vec_push(uint64_t, stack, x);
  while(tmap_vec_size(stack)) {
      x = tmap_vec_pop(stack);
      k = x>>32; l = (uint32_t)x;
      tmap_bwtl_2occ4(b, k-1, l, cntk, cntl);
      for(j = 0; j != 4; ++j) {
          k = b->L2[j] + cntk[j] + 1;
          l = b->L2[j] + cntl[j];
          if(k > l) continue;
          x = (uint64_t)k << 32 | l;
          iter = tmap_hash_put(64, h, x, &ret);
          if(ret) { // if not present
              tmap_hash_value(h, iter) = 1;
              tmap_vec_push(uint64_t, stack, x);
          } else ++tmap_hash_value(h, iter);
      }
  }
  tmap_vec_destroy(stack);
  return h;
}

// pick up top T matches at a node
static void 
tmap_map2_core_cut_tail(tmap_map2_entry_t *u, int32_t T, tmap_map2_entry_t *aux)
{
  int32_t i, *a, n, x;
  if(u->n <= T) return;
  if(aux->max < u->n) {
      aux->max = u->n;
      aux->array = tmap_realloc(aux->array, aux->max * sizeof(tmap_map2_cell_t), "aux->array");
  }
  a = (int*)aux->array;
  for(i = n = 0; i != u->n; ++i)
    if(u->array[i].match_sa.l && u->array[i].G > TMAP_MAP2_MINUS_INF)
      a[n++] = -u->array[i].G;
  if(n <= T) return;
  x = -tmap_sort_small(int32_t, n, a, T);
  n = 0;
  for(i = 0; i < u->n; ++i) {
      tmap_map2_cell_t *p = u->array + i;
      if(p->G == x) ++n;
      if(p->G < x || (p->G == x && n >= T)) {
          p->match_sa.k = p->match_sa.l = 0; p->G = TMAP_MAP2_MINUS_INF;
          p->match_sa.hi = p->match_sa.offset = 0;
          if(p->ppos >= 0) u->array[p->ppos].cpos[p->pj] = -1;
      }
  }
}
// remove duplicated cells
static inline void 
tmap_map2_core_remove_duplicate(tmap_map2_entry_t *u, tmap_hash_t(tmap_map2_qintv) *hash)
{
  int32_t i, ret, j;
  tmap_hash_iter_t k;
  tmap_map2_qintv_t key;
  tmap_hash_clear(tmap_map2_qintv, hash);
  for(i = 0; i != u->n; ++i) {
      tmap_map2_cell_t *p = u->array + i;
      if(p->match_sa.l== 0) continue;
      key.k = p->match_sa.k; key.l = p->match_sa.l;
      k = tmap_hash_put(tmap_map2_qintv, hash, key, &ret);
      j = -1;
      if(ret == 0) {
          if((uint32_t)tmap_hash_value(hash, k) >= p->G) j = i;
          else {
              j = tmap_hash_value(hash, k)>>32;
              tmap_hash_value(hash, k) = (uint64_t)i<<32 | p->G;
          }
      } else tmap_hash_value(hash, k) = (uint64_t)i<<32 | p->G;
      if(j >= 0) {
          p = u->array + j;
          p->match_sa.k = p->match_sa.l = 0; p->G = TMAP_MAP2_MINUS_INF;
          p->match_sa.hi = p->match_sa.offset = 0;
          if(p->ppos >= 0) u->array[p->ppos].cpos[p->pj] = -3;
      }
  }
}

// merge two entries
static void 
tmap_map2_core_merge_entry(const tmap_map_opt_t * __restrict opt, 
                           tmap_map2_entry_t *u, tmap_map2_entry_t *v)
{
  int32_t i;
  if(u->n + v->n >= u->max) {
      u->max = u->n + v->n;
      u->array = tmap_realloc(u->array, u->max * sizeof(tmap_map2_cell_t), "u->array");
  }
  for(i = 0; i != v->n; ++i) {
      tmap_map2_cell_t *p = v->array + i;
      if(p->ppos >= 0) p->ppos += u->n;
      if(p->cpos[0] >= 0) p->cpos[0] += u->n;
      if(p->cpos[1] >= 0) p->cpos[1] += u->n;
      if(p->cpos[2] >= 0) p->cpos[2] += u->n;
      if(p->cpos[3] >= 0) p->cpos[3] += u->n;
  }
  memcpy(u->array + u->n, v->array, v->n * sizeof(tmap_map2_cell_t));
  u->n += v->n;
}

static inline tmap_map2_cell_t *
tmap_map2_core_push_array_p(tmap_map2_entry_t *e)
{
  if(e->n == e->max) {
      e->max = e->max? e->max<<1 : 256;
      e->array = tmap_realloc(e->array, e->max * sizeof(tmap_map2_cell_t), "e->array");
  }
  return e->array + e->n;
}
/* --- END: utilities --- */

/* --- BEGIN: processing partial hits --- */
static void 
tmap_map2_core_save_hits(const tmap_bwtl_t *bwtl, int32_t thres, tmap_map2_hit_t *hits, tmap_map2_entry_t *u)
{
  int32_t i;
  uint32_t k; // for the target
  for(i = 0; i < u->n; ++i) {
      tmap_map2_cell_t *p = u->array + i;
      if(p->G < thres) continue;
      for(k = u->tk; k <= u->tl; ++k) {
          int32_t beg, end;
          tmap_map2_hit_t *q = NULL;
          beg = bwtl->sa[k]; end = beg + p->tlen - 1;
          if(p->G > hits[beg*2].G) {
              hits[beg*2+1] = hits[beg*2];
              q = hits + beg * 2;
          } else if(p->G > hits[beg*2+1].G) q = hits + beg * 2 + 1;
          if(q) {
              q->k = p->match_sa.k; q->l = p->match_sa.l; q->qlen = p->qlen; q->tlen = p->tlen; q->G = p->G;
              q->beg = beg; q->end = end; q->G2 = q->k == q->l? 0 : q->G;
              q->flag = q->n_seeds = 0;
          }
      }
  }
}

/* 
 "narrow hits" are node-to-node hits that have a high score and
 are not so repetitive (|SA interval|<=IS). 
 */
static void 
tmap_map2_save_narrow_hits(const tmap_bwtl_t *bwtl, tmap_map2_entry_t *u, tmap_map2_aln_t *b1, int32_t t, int32_t IS)
{
  int32_t i;
  for(i = 0; i < u->n; ++i) {
      tmap_map2_hit_t *q;
      tmap_map2_cell_t *p = u->array + i;
      if(p->G >= t && p->match_sa.l - p->match_sa.k + 1 <= IS) { // good narrow hit
          if(b1->max <= b1->n) { // check if we need more memory
              tmap_map2_aln_realloc(b1, b1->n+1);
          }
          q = &b1->hits[b1->n++];
          q->k = p->match_sa.k; q->l = p->match_sa.l;
          q->qlen = p->qlen;
          q->tlen = p->tlen;
          q->G = p->G; q->G2 = 0;
          q->beg = bwtl->sa[u->tk]; q->end = q->beg + p->tlen - 1;
          q->flag = 0;
          // delete p
          p->match_sa.k = p->match_sa.l = 0; p->G = TMAP_MAP2_MINUS_INF;
          p->match_sa.hi = p->match_sa.offset = 0;
          if(p->ppos >= 0) u->array[p->ppos].cpos[p->pj] = -3;
      }
  }
}
/* --- END: processing partial hits --- */

static inline int32_t 
tmap_map2_core_fill_cell(const tmap_map_opt_t *opt, int32_t match_score, tmap_map2_cell_t *c[4])
{
  int32_t G = c[3]? c[3]->G + match_score : TMAP_MAP2_MINUS_INF;
  if(c[1]) {
      c[0]->I = c[1]->I > c[1]->G - opt->pen_gapo? c[1]->I - opt->pen_gape : c[1]->G - opt->pen_gapo - opt->pen_gape;
      if(c[0]->I > G) G = c[0]->I;
  } else c[0]->I = TMAP_MAP2_MINUS_INF;
  if(c[2]) {
      c[0]->D = c[2]->D > c[2]->G - opt->pen_gapo? c[2]->D - opt->pen_gape : c[2]->G - opt->pen_gapo - opt->pen_gape;
      if(c[0]->D > G) G = c[0]->D;
  } else c[0]->D = TMAP_MAP2_MINUS_INF;
  return(c[0]->G = G);
}

static void 
tmap_map2_core_init(const tmap_bwtl_t *target, const tmap_bwt_t *query_bwt, tmap_map2_stack_t *s)
{
  tmap_map2_entry_t *u;
  tmap_map2_cell_t *x;

  u = tmap_map2_mempool_pop(s->pool);
  u->tk = 0; u->tl = target->seq_len;
  x = tmap_map2_core_push_array_p(u);
  *x = tmap_map2_core_default_cell;
  x->G = 0; // set to zero, no TMAP_MAP2_MINUS_INF
  x->match_sa.k = 0; x->match_sa.l = query_bwt->seq_len;
  x->match_sa.hi = 0; x->match_sa.offset = 0;
  u->n++;
  tmap_map2_stack_push0(s, u);
}

/*
 On return, ret[1] keeps not-so-repetitive hits (narrow SA hits); ret[0]
 keeps all hits (right?) 
Note: tlen may be over-estimated!
 */
tmap_map2_aln_t **
tmap_map2_core_aln(const tmap_map_opt_t *opt, const tmap_bwtl_t *target, 
               const tmap_refseq_t *refseq,
               const tmap_bwt_t *query_bwt, 
               const tmap_sa_t *query_sa, 
               tmap_bwt_match_hash_t *query_hash, 
               tmap_map2_global_mempool_t *pool)
{
  tmap_map2_stack_t *stack = (tmap_map2_stack_t*)pool->stack;
  tmap_map2_aln_t *b, *b1, **b_ret;
  int32_t i, j, score_mat[16], *heap, heap_size, n_tot = 0;
  tmap_hash_t(tmap_map2_qintv) *rhash;
  tmap_hash_t(64) *chash;

  // initialize connectivity hash (chash)
  chash = tmap_map2_core_connectivity(target);
  // calculate score matrix
  for(i = 0; i != 4; ++i)
    for(j = 0; j != 4; ++j)
      score_mat[i<<2|j] = (i == j)? opt->score_match : -opt->pen_mm;
  // initialize other variables
  rhash = tmap_hash_init(tmap_map2_qintv);
  tmap_map2_core_init(target, query_bwt, stack);
  heap_size = opt->z_best;
  heap = tmap_calloc(heap_size, sizeof(int32_t), "heap");
  // initialize the return struct
  b = tmap_map2_aln_init();
  b1 = tmap_map2_aln_init();
  tmap_map2_aln_realloc(b, target->seq_len * 2);
  b->n = target->seq_len * 2; // why?
  for(i=0;i<b->n;i++) { // set the default score to -inf
      b->hits[i].G = TMAP_MAP2_MINUS_INF;
  }
  b_ret = tmap_calloc(2, sizeof(tmap_map2_aln_t*), "b_ret");
  b_ret[0] = b; b_ret[1] = b1;
  // the main loop: traversal of the DAG
  while(!tmap_map2_stack_isempty(stack)) {
      int32_t old_n, tj;
      tmap_map2_entry_t *v;
      tmap_bwt_int_t k, l;
      uint32_t tcntk[4], tcntl[4];

      v = tmap_map2_stack_pop(stack); old_n = v->n;
      n_tot += v->n;

      for(i = 0; i < v->n; ++i) { // test max depth and band width
          tmap_map2_cell_t *p = v->array + i;
          if(p->match_sa.l == 0) continue;
          if(p->tlen - (int)p->qlen > opt->bw || (int)p->qlen - p->tlen > opt->bw) {
              p->match_sa.k = p->match_sa.l = 0;
              p->match_sa.hi = p->match_sa.offset = 0;
              if(p->ppos >= 0) v->array[p->ppos].cpos[p->pj] = -5;
          }
      }

      // get Occ for the DAG
      tmap_bwtl_2occ4(target, v->tk - 1, v->tl, tcntk, tcntl);
      for(tj = 0; tj != 4; ++tj) { // descend to the children
          tmap_bwt_match_occ_t qnext[4];
          int32_t qj, *curr_score_mat = score_mat + tj * 4;
          tmap_hash_iter_t iter;
          tmap_map2_entry_t *u;

          k = target->L2[tj] + tcntk[tj] + 1;
          l = target->L2[tj] + tcntl[tj];
          if(k > l) continue;
          // update counter
          iter = tmap_hash_get(64, chash, (uint64_t)k<<32 | l);
          --tmap_hash_value(chash, iter);
          // initialization
          u = tmap_map2_mempool_pop(stack->pool);
          u->tk = k; u->tl = l;
          memset(heap, 0, sizeof(int) * opt->z_best);
          // loop through all the nodes in v
          for(i = 0; i < v->n; ++i) {
              tmap_map2_cell_t *p = v->array + i, *x, *c[4]; // c[0]=>current, c[1]=>I, c[2]=>D, c[3]=>G
              int32_t is_added = 0;
              if(p->match_sa.l == 0) continue; // deleted node
              c[0] = x = tmap_map2_core_push_array_p(u);
              x->G = TMAP_MAP2_MINUS_INF;
              p->upos = x->upos = -1;
              if(p->ppos >= 0) { // parent has been visited
                  c[1] = (v->array[p->ppos].upos >= 0)? u->array + v->array[p->ppos].upos : 0;
                  c[3] = v->array + p->ppos; c[2] = p;
                  if(tmap_map2_core_fill_cell(opt, curr_score_mat[p->pj], c) > TMAP_MAP2_MINUS_INF) { // then update topology at p and x
                      x->ppos = v->array[p->ppos].upos; // the parent pos in u
                      p->upos = u->n++; // the current pos in u
                      if(x->ppos >= 0) u->array[x->ppos].cpos[p->pj] = p->upos; // the child pos of its parent in u
                      is_added = 1;
                  }
              } else {
                  x->D = p->D > p->G - opt->pen_gapo? p->D - opt->pen_gape : p->G - opt->pen_gapo - opt->pen_gape;
                  if(x->D > TMAP_MAP2_MINUS_INF) {
                      x->G = x->D;
                      x->I = TMAP_MAP2_MINUS_INF; x->ppos = -1;
                      p->upos = u->n++;
                      is_added = 1;
                  }
              }
              if(is_added) { // x has been added to u->array. fill the remaining variables
                  x->cpos[0] = x->cpos[1] = x->cpos[2] = x->cpos[3] = -1;
                  x->pj = p->pj; 
                  x->match_sa = p->match_sa;
                  x->qlen = p->qlen; x->tlen = p->tlen + 1;
                  if(x->G > -heap[0]) {
                      heap[0] = -x->G;
                      tmap_sort_heapadjust(int32_t, 0, heap_size, heap);
                  }
              }
              if((x->G > opt->pen_gapo + opt->pen_gape && x->G >= -heap[0]) || i < old_n) { // good node in u, or in v
                  if(p->cpos[0] == -1 || p->cpos[1] == -1 || p->cpos[2] == -1 || p->cpos[3] == -1) {
                      tmap_bwt_match_hash_2occ4(query_bwt, &p->match_sa, qnext, query_hash);
                      for(qj = 0; qj != 4; ++qj) { // descend to the prefix trie
                          if(p->cpos[qj] != -1) continue; // this node will be visited later
                          k = qnext[qj].k;
                          l = qnext[qj].l;
                          if(k > l) { p->cpos[qj] = -2; continue; }
                          x = tmap_map2_core_push_array_p(v);
                          p = v->array + i; // p may not point to the correct position after realloc
                          x->G = x->I = x->D = TMAP_MAP2_MINUS_INF;
                          x->match_sa = qnext[qj];
                          x->pj = qj; x->qlen = p->qlen + 1; x->ppos = i; x->tlen = p->tlen;
                          x->cpos[0] = x->cpos[1] = x->cpos[2] = x->cpos[3] = -1;
                          p->cpos[qj] = v->n++;
                      } // ~for(qj)
                  } // ~if(p->cpos[])
              } // ~if
          } // ~for(i)
          if(u->n) tmap_map2_core_save_hits(target, opt->score_thr, b->hits, u);
            { // push u to the stack (or to the pending array)
              uint32_t cnt, pos;
              cnt = (uint32_t)tmap_hash_value(chash, iter);
              pos = tmap_hash_value(chash, iter)>>32;
              if(pos) { // something in the pending array, then merge
                  tmap_map2_entry_t *w = tmap_vec_A(stack->pending, pos-1);
                  if(u->n) {
                      if(w->n < u->n) { // swap
                          w = u; u = tmap_vec_A(stack->pending, pos-1); tmap_vec_A(stack->pending, pos-1) = w;
                      }
                      tmap_map2_core_merge_entry(opt, w, u);
                  }
                  if(cnt == 0) { // move from pending to stack0
                      tmap_map2_core_remove_duplicate(w, rhash);
                      tmap_map2_save_narrow_hits(target, w, b1, opt->score_thr, opt->max_seed_intv);
                      tmap_map2_core_cut_tail(w, opt->z_best, u);
                      tmap_map2_stack_push0(stack, w);
                      tmap_vec_A(stack->pending, pos-1) = 0;
                      --stack->n_pending;
                  }
                  tmap_map2_mempool_push(stack->pool, u);
              } else if(cnt) { // the first time
                  if(u->n) { // push to the pending queue
                      ++stack->n_pending;
                      tmap_vec_push(tmap_map2_entry_p, stack->pending, u);
                      tmap_hash_value(chash, iter) = (uint64_t)tmap_vec_size(stack->pending)<<32 | cnt;
                  } else tmap_map2_mempool_push(stack->pool, u);
              } else { // cnt == 0, then push to the stack
                  tmap_map2_entry_t *w = tmap_map2_mempool_pop(stack->pool);
                  tmap_map2_save_narrow_hits(target, u, b1, opt->score_thr, opt->max_seed_intv);
                  tmap_map2_core_cut_tail(u, opt->z_best, w);
                  tmap_map2_mempool_push(stack->pool, w);
                  tmap_map2_stack_push0(stack, u);
              }
            }
      } // ~for(tj)
      tmap_map2_mempool_push(stack->pool, v);
  } // while(top)
  // free
  free(heap);
  tmap_hash_destroy(tmap_map2_qintv, rhash);
  tmap_hash_destroy(64, chash);
  stack->pending.n = stack->stack0.n = 0;

  return b_ret;
}
