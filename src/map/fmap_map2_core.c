#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/resource.h>
#include "../util/fmap_alloc.h"
#include "../util/fmap_sort.h"
#include "../util/fmap_vec.h"
#include "../util/fmap_hash.h"
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_bwtl.h"
#include "../index/fmap_sa.h"
#include "fmap_map_util.h"
#include "fmap_map2.h"
#include "fmap_map2_aux.h"
#include "fmap_map2_core.h"

// for hashing
FMAP_HASH_MAP_INIT_INT64(64, uint64_t)

// for sorting generically
FMAP_SORT_INIT_GENERIC(int32_t)

// TODO: document
static const fmap_map2_cell_t fmap_map2_core_default_cell = { {0, 0, 0, 0}, FMAP_MAP2_MINUS_INF, FMAP_MAP2_MINUS_INF,
     FMAP_MAP2_MINUS_INF, 0, 0, 0, -1, -1, {-1, -1, -1, -1} };

/* --- BEGIN: utilities --- */
static fmap_hash_t(64) *
fmap_map2_core_connectivity(const fmap_bwtl_t *b)
{
  fmap_hash_t(64) *h;
  uint32_t k, l, cntk[4], cntl[4];
  uint64_t x;
  fmap_hash_iter_t iter;
  int32_t j, ret;
  fmap_vec_t(uint64_t) stack;

  fmap_vec_init(stack);
  h = fmap_hash_init(64);
  fmap_hash_resize(64, h, b->seq_len * 4);
  x = b->seq_len;
  fmap_vec_push(uint64_t, stack, x);
  while(fmap_vec_size(stack)) {
      x = fmap_vec_pop(stack);
      k = x>>32; l = (uint32_t)x;
      fmap_bwtl_2occ4(b, k-1, l, cntk, cntl);
      for(j = 0; j != 4; ++j) {
          k = b->L2[j] + cntk[j] + 1;
          l = b->L2[j] + cntl[j];
          if(k > l) continue;
          x = (uint64_t)k << 32 | l;
          iter = fmap_hash_put(64, h, x, &ret);
          if(ret) { // if not present
              fmap_hash_value(h, iter) = 1;
              fmap_vec_push(uint64_t, stack, x);
          } else ++fmap_hash_value(h, iter);
      }
  }
  fmap_vec_destroy(stack);
  return h;
}

// pick up top T matches at a node
static void 
fmap_map2_core_cut_tail(fmap_map2_entry_t *u, int32_t T, fmap_map2_entry_t *aux)
{
  int32_t i, *a, n, x;
  if(u->n <= T) return;
  if(aux->max < u->n) {
      aux->max = u->n;
      aux->array = fmap_realloc(aux->array, aux->max * sizeof(fmap_map2_cell_t), "aux->array");
  }
  a = (int*)aux->array;
  for(i = n = 0; i != u->n; ++i)
    if(u->array[i].match_sa.l && u->array[i].G > 0)
      a[n++] = -u->array[i].G;
  if(n <= T) return;
  x = -fmap_sort_small(int32_t, n, a, T);
  n = 0;
  for(i = 0; i < u->n; ++i) {
      fmap_map2_cell_t *p = u->array + i;
      if(p->G == x) ++n;
      if(p->G < x || (p->G == x && n >= T)) {
          p->match_sa.k = p->match_sa.l = 0; p->G = 0;
          p->match_sa.hi = p->match_sa.offset = 0;
          if(p->ppos >= 0) u->array[p->ppos].cpos[p->pj] = -1;
      }
  }
}
// remove duplicated cells
static inline void 
fmap_map2_core_remove_duplicate(fmap_map2_entry_t *u, fmap_hash_t(64) *hash)
{
  int32_t i, ret, j;
  fmap_hash_iter_t k;
  uint64_t key;
  fmap_hash_clear(64, hash);
  for(i = 0; i != u->n; ++i) {
      fmap_map2_cell_t *p = u->array + i;
      if(p->match_sa.l== 0) continue;
      key = (uint64_t)p->match_sa.k << 32 | p->match_sa.l;
      k = fmap_hash_put(64, hash, key, &ret);
      j = -1;
      if(ret == 0) {
          if((uint32_t)fmap_hash_value(hash, k) >= p->G) j = i;
          else {
              j = fmap_hash_value(hash, k)>>32;
              fmap_hash_value(hash, k) = (uint64_t)i<<32 | p->G;
          }
      } else fmap_hash_value(hash, k) = (uint64_t)i<<32 | p->G;
      if(j >= 0) {
          p = u->array + j;
          p->match_sa.k = p->match_sa.l = 0; p->G = 0;
          p->match_sa.hi = p->match_sa.offset = 0;
          if(p->ppos >= 0) u->array[p->ppos].cpos[p->pj] = -3;
      }
  }
}

// merge two entries
static void 
fmap_map2_core_merge_entry(const fmap_map2_opt_t * __restrict opt, 
                           fmap_map2_entry_t *u, fmap_map2_entry_t *v)
{
  int32_t i;
  if(u->n + v->n >= u->max) {
      u->max = u->n + v->n;
      u->array = fmap_realloc(u->array, u->max * sizeof(fmap_map2_cell_t), "u->array");
  }
  for(i = 0; i != v->n; ++i) {
      fmap_map2_cell_t *p = v->array + i;
      if(p->ppos >= 0) p->ppos += u->n;
      if(p->cpos[0] >= 0) p->cpos[0] += u->n;
      if(p->cpos[1] >= 0) p->cpos[1] += u->n;
      if(p->cpos[2] >= 0) p->cpos[2] += u->n;
      if(p->cpos[3] >= 0) p->cpos[3] += u->n;
  }
  memcpy(u->array + u->n, v->array, v->n * sizeof(fmap_map2_cell_t));
  u->n += v->n;
}

static inline fmap_map2_cell_t *
fmap_map2_core_push_array_p(fmap_map2_entry_t *e)
{
  if(e->n == e->max) {
      e->max = e->max? e->max<<1 : 256;
      e->array = fmap_realloc(e->array, e->max * sizeof(fmap_map2_cell_t), "e->array");
  }
  return e->array + e->n;
}
/* --- END: utilities --- */

/* --- BEGIN: processing partial hits --- */
static void 
fmap_map2_core_save_hits(const fmap_bwtl_t *bwtl, int32_t thres, fmap_map2_hit_t *hits, fmap_map2_entry_t *u)
{
  int32_t i;
  uint32_t k;
  for(i = 0; i < u->n; ++i) {
      fmap_map2_cell_t *p = u->array + i;
      if(p->G < thres) continue;
      for(k = u->tk; k <= u->tl; ++k) {
          int32_t beg, end;
          fmap_map2_hit_t *q = NULL;
          beg = bwtl->sa[k]; end = beg + p->tlen;
          if(p->G > hits[beg*2].G) {
              hits[beg*2+1] = hits[beg*2];
              q = hits + beg * 2;
          } else if(p->G > hits[beg*2+1].G) q = hits + beg * 2 + 1;
          if(q) {
              q->k = p->match_sa.k; q->l = p->match_sa.l; q->len = p->qlen; q->G = p->G;
              q->beg = beg; q->end = end; q->G2 = q->k == q->l? 0 : q->G;
              q->flag = q->n_seeds = 0;
          }
      }
  }
}

static void 
fmap_map2_save_narrow_hits(const fmap_bwtl_t *bwtl, fmap_map2_entry_t *u, fmap_map2_aln_t *b1, int32_t t, int32_t IS)
{
  int32_t i;
  for(i = 0; i < u->n; ++i) {
      fmap_map2_hit_t *q;
      fmap_map2_cell_t *p = u->array + i;
      if(p->G >= t && p->match_sa.l - p->match_sa.k + 1 <= IS) { // good narrow hit
          if(b1->max == b1->n) {
              b1->max = b1->max? b1->max<<1 : 4;
              b1->hits = fmap_realloc(b1->hits, b1->max * sizeof(fmap_map2_hit_t), "b1->hits");
          }
          q = &b1->hits[b1->n++];
          q->k = p->match_sa.k; q->l = p->match_sa.l;
          q->len = p->qlen;
          q->G = p->G; q->G2 = 0;
          q->beg = bwtl->sa[u->tk]; q->end = q->beg + p->tlen;
          q->flag = 0;
          // delete p
          p->match_sa.k = p->match_sa.l = 0; p->G = 0;
          p->match_sa.hi = p->match_sa.offset = 0;
          if(p->ppos >= 0) u->array[p->ppos].cpos[p->pj] = -3;
      }
  }
}
/* --- END: processing partial hits --- */

static inline int32_t 
fmap_map2_core_fill_cell(const fmap_map2_opt_t *opt, int32_t match_score, fmap_map2_cell_t *c[4])
{
  int32_t G = c[3]? c[3]->G + match_score : FMAP_MAP2_MINUS_INF;
  if(c[1]) {
      c[0]->I = c[1]->I > c[1]->G - opt->pen_gapo? c[1]->I - opt->pen_gape : c[1]->G - opt->pen_gapo - opt->pen_gape;
      if(c[0]->I > G) G = c[0]->I;
  } else c[0]->I = FMAP_MAP2_MINUS_INF;
  if(c[2]) {
      c[0]->D = c[2]->D > c[2]->G - opt->pen_gapo? c[2]->D - opt->pen_gape : c[2]->G - opt->pen_gapo - opt->pen_gape;
      if(c[0]->D > G) G = c[0]->D;
  } else c[0]->D = FMAP_MAP2_MINUS_INF;
  return(c[0]->G = G);
}

static void 
fmap_map2_core_init(const fmap_bwtl_t *target, const fmap_bwt_t *query_bwt, fmap_map2_stack_t *s)
{
  fmap_map2_entry_t *u;
  fmap_map2_cell_t *x;

  u = fmap_map2_mempool_pop(s->pool);
  u->tk = 0; u->tl = target->seq_len;
  x = fmap_map2_core_push_array_p(u);
  *x = fmap_map2_core_default_cell;
  x->G = 0; 
  x->match_sa.k = 0; x->match_sa.l = query_bwt->seq_len;
  x->match_sa.hi = 0; x->match_sa.offset = 0;
  u->n++;
  fmap_map2_stack_push0(s, u);
}

fmap_map2_aln_t **
fmap_map2_core_aln(const fmap_map2_opt_t *opt, const fmap_bwtl_t *target, 
               const fmap_bwt_t *query_bwt, const fmap_sa_t *query_sa, 
               fmap_map2_global_mempool_t *pool)
{
  fmap_map2_stack_t *stack = (fmap_map2_stack_t*)pool->stack;
  fmap_map2_aln_t *b, *b1, **b_ret;
  int32_t i, j, score_mat[16], *heap, heap_size, n_tot = 0;
  fmap_hash_t(64) *rhash, *chash;

  // initialize connectivity hash (chash)
  chash = fmap_map2_core_connectivity(target);
  // calculate score matrix
  for(i = 0; i != 4; ++i)
    for(j = 0; j != 4; ++j)
      score_mat[i<<2|j] = (i == j)? opt->score_match : -opt->pen_mm;
  // initialize other variables
  rhash = fmap_hash_init(64);
  fmap_map2_core_init(target, query_bwt, stack);
  heap_size = opt->z_best;
  heap = fmap_calloc(heap_size, sizeof(int32_t), "heap");
  // initialize the return struct
  b = fmap_calloc(1, sizeof(fmap_map2_aln_t), "b");
  b->n = b->max = target->seq_len * 2;
  b->hits = fmap_calloc(b->max, sizeof(fmap_map2_hit_t), "b->hits");
  b1 = fmap_calloc(1, sizeof(fmap_map2_aln_t), "b1");
  b_ret = fmap_calloc(2, sizeof(fmap_map2_aln_t*), "b_ret");
  b_ret[0] = b; b_ret[1] = b1;
  // the main loop: traversal of the DAG
  while(!fmap_map2_stack_isempty(stack)) {
      int32_t old_n, tj;
      fmap_map2_entry_t *v;
      uint32_t k, l, tcntk[4], tcntl[4];

      v = fmap_map2_stack_pop(stack); old_n = v->n;
      n_tot += v->n;

      for(i = 0; i < v->n; ++i) { // test max depth and band width
          fmap_map2_cell_t *p = v->array + i;
          if(p->match_sa.l == 0) continue;
          if(p->tlen - (int)p->qlen > opt->bw || (int)p->qlen - p->tlen > opt->bw) {
              p->match_sa.k = p->match_sa.l = 0;
              p->match_sa.hi = p->match_sa.offset = 0;
              if(p->ppos >= 0) v->array[p->ppos].cpos[p->pj] = -5;
          }
      }

      // get Occ for the DAG
      fmap_bwtl_2occ4(target, v->tk - 1, v->tl, tcntk, tcntl);
      for(tj = 0; tj != 4; ++tj) { // descend to the children
          fmap_bwt_match_occ_t qnext[4];
          int32_t qj, *curr_score_mat = score_mat + tj * 4;
          fmap_hash_iter_t iter;
          fmap_map2_entry_t *u;

          k = target->L2[tj] + tcntk[tj] + 1;
          l = target->L2[tj] + tcntl[tj];
          if(k > l) continue;
          // update counter
          iter = fmap_hash_get(64, chash, (uint64_t)k<<32 | l);
          --fmap_hash_value(chash, iter);
          // initialization
          u = fmap_map2_mempool_pop(stack->pool);
          u->tk = k; u->tl = l;
          memset(heap, 0, sizeof(int) * opt->z_best);
          // loop through all the nodes in v
          for(i = 0; i < v->n; ++i) {
              fmap_map2_cell_t *p = v->array + i, *x, *c[4]; // c[0]=>current, c[1]=>I, c[2]=>D, c[3]=>G
              int32_t is_added = 0;
              if(p->match_sa.l == 0) continue; // deleted node
              c[0] = x = fmap_map2_core_push_array_p(u);
              x->G = FMAP_MAP2_MINUS_INF;
              p->upos = x->upos = -1;
              if(p->ppos >= 0) { // parent has been visited
                  c[1] = (v->array[p->ppos].upos >= 0)? u->array + v->array[p->ppos].upos : 0;
                  c[3] = v->array + p->ppos; c[2] = p;
                  if(fmap_map2_core_fill_cell(opt, curr_score_mat[p->pj], c) > 0) { // then update topology at p and x
                      x->ppos = v->array[p->ppos].upos; // the parent pos in u
                      p->upos = u->n++; // the current pos in u
                      if(x->ppos >= 0) u->array[x->ppos].cpos[p->pj] = p->upos; // the child pos of its parent in u
                      is_added = 1;
                  }
              } else {
                  x->D = p->D > p->G - opt->pen_gapo? p->D - opt->pen_gape : p->G - opt->pen_gapo - opt->pen_gape;
                  if(x->D > 0) {
                      x->G = x->D;
                      x->I = FMAP_MAP2_MINUS_INF; x->ppos = -1;
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
                      fmap_sort_heapadjust(int32_t, 0, heap_size, heap);
                  }
              }
              if((x->G > opt->pen_gapo + opt->pen_gape && x->G >= -heap[0]) || i < old_n) { // good node in u, or in v
                  if(p->cpos[0] == -1 || p->cpos[1] == -1 || p->cpos[2] == -1 || p->cpos[3] == -1) {
                      fmap_bwt_match_2occ4(query_bwt, &p->match_sa, qnext);
                      for(qj = 0; qj != 4; ++qj) { // descend to the prefix trie
                          if(p->cpos[qj] != -1) continue; // this node will be visited later
                          k = qnext[qj].k;
                          l = qnext[qj].l;
                          if(k > l) { p->cpos[qj] = -2; continue; }
                          x = fmap_map2_core_push_array_p(v);
                          p = v->array + i; // p may not point to the correct position after realloc
                          x->G = x->I = x->D = FMAP_MAP2_MINUS_INF;
                          x->match_sa = qnext[qj];
                          x->pj = qj; x->qlen = p->qlen + 1; x->ppos = i; x->tlen = p->tlen;
                          x->cpos[0] = x->cpos[1] = x->cpos[2] = x->cpos[3] = -1;
                          p->cpos[qj] = v->n++;
                      } // ~for(qj)
                  } // ~if(p->cpos[])
              } // ~if
          } // ~for(i)
          if(u->n) fmap_map2_core_save_hits(target, opt->score_thr, b->hits, u);
            { // push u to the stack (or to the pending array)
              uint32_t cnt, pos;
              cnt = (uint32_t)fmap_hash_value(chash, iter);
              pos = fmap_hash_value(chash, iter)>>32;
              if(pos) { // something in the pending array, then merge
                  fmap_map2_entry_t *w = fmap_vec_A(stack->pending, pos-1);
                  if(u->n) {
                      if(w->n < u->n) { // swap
                          w = u; u = fmap_vec_A(stack->pending, pos-1); fmap_vec_A(stack->pending, pos-1) = w;
                      }
                      fmap_map2_core_merge_entry(opt, w, u);
                  }
                  if(cnt == 0) { // move from pending to stack0
                      fmap_map2_core_remove_duplicate(w, rhash);
                      fmap_map2_save_narrow_hits(target, w, b1, opt->score_thr, opt->max_seed_intv);
                      fmap_map2_core_cut_tail(w, opt->z_best, u);
                      fmap_map2_stack_push0(stack, w);
                      fmap_vec_A(stack->pending, pos-1) = 0;
                      --stack->n_pending;
                  }
                  fmap_map2_mempool_push(stack->pool, u);
              } else if(cnt) { // the first time
                  if(u->n) { // push to the pending queue
                      ++stack->n_pending;
                      fmap_vec_push(fmap_map2_entry_p, stack->pending, u);
                      fmap_hash_value(chash, iter) = (uint64_t)fmap_vec_size(stack->pending)<<32 | cnt;
                  } else fmap_map2_mempool_push(stack->pool, u);
              } else { // cnt == 0, then push to the stack
                  fmap_map2_entry_t *w = fmap_map2_mempool_pop(stack->pool);
                  fmap_map2_save_narrow_hits(target, u, b1, opt->score_thr, opt->max_seed_intv);
                  fmap_map2_core_cut_tail(u, opt->z_best, w);
                  fmap_map2_mempool_push(stack->pool, w);
                  fmap_map2_stack_push0(stack, u);
              }
            }
      } // ~for(tj)
      fmap_map2_mempool_push(stack->pool, v);
  } // while(top)
  fmap_map2_aux_resolve_duphits(query_bwt, query_sa, b, opt->max_seed_intv, 0); 
  fmap_map2_aux_resolve_duphits(query_bwt, query_sa, b1, opt->max_seed_intv, 0); 
  // free
  free(heap);
  fmap_hash_destroy(64, rhash);
  fmap_hash_destroy(64, chash);
  stack->pending.n = stack->stack0.n = 0;
  return b_ret;
}
