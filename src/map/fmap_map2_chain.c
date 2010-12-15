#include <stdio.h>
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_sa.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_sort.h"
#include "../seq/fmap_seq.h"
#include "fmap_map_util.h"
#include "fmap_map2.h"
#include "fmap_map2_aux.h"
#include "fmap_map2_chain.h"

#define _fmap_map2_chain_lt(a, b) ((a).qbeg < (b).qbeg)
FMAP_SORT_INIT(fmap_map2_chain, fmap_map2_chain_t, _fmap_map2_chain_lt)

static int 
fmap_map2_chain_chaining(const fmap_map_opt_t *opt, int shift, int n, fmap_map2_chain_t *z, fmap_map2_chain_t *chain)
{
  int32_t j, k, m = 0;
  fmap_sort_introsort(fmap_map2_chain, n, z);
  for (j = 0; j < n; ++j) {
      fmap_map2_chain_t *p = z + j;
      for (k = m - 1; k >= 0; --k) {
          fmap_map2_chain_t *q = chain + k;
          int x = p->qbeg - q->qbeg; // always positive
          int y = p->tbeg - q->tbeg;
          if (y > 0 && x - y <= opt->bw && y - x <= opt->bw) {
              if (p->qend > q->qend) q->qend = p->qend;
              if (p->tend > q->tend) q->tend = p->tend;
              ++q->chain;
              p->chain = shift + k;
              break;
          }
      }
      if (k < 0) {
          chain[m] = *p;
          chain[m].chain = 1;
          chain[m].idx = p->chain = shift + m;
          ++m;
      }
  }
  return m;
}

void 
fmap_map2_chain_filter(const fmap_map_opt_t *opt, int len, fmap_map2_aln_t *b[2])
{
  fmap_map2_chain_t *z[2], *chain[2];
  int32_t i, j, k, n[2], m[2];
  char *flag;
  // initialization
  n[0] = b[0]->n; n[1] = b[1]->n;
  z[0] = fmap_calloc(n[0] + n[1], sizeof(fmap_map2_chain_t), "z[0]");
  z[1] = z[0] + n[0];
  chain[0] = fmap_calloc(n[0] + n[1], sizeof(fmap_map2_chain_t), "chain[0]");
  for (k = j = 0; k < 2; ++k) {
      for (i = 0; i < b[k]->n; ++i) {
          fmap_map2_hit_t *p = b[k]->hits + i;
          fmap_map2_chain_t *q = z[k] + i;
          q->flag = k; q->idx = i;
          q->tbeg = p->k; q->tend = p->k + p->len;
          q->chain = -1;
          q->qbeg = p->beg; q->qend = p->end;
      }
  }
  // chaining
  m[0] = fmap_map2_chain_chaining(opt, 0,    n[0], z[0], chain[0]);
  chain[1] = chain[0] + m[0];
  m[1] = fmap_map2_chain_chaining(opt, m[0], n[1], z[1], chain[1]);	
  // change query coordinate on the reverse strand
  for (k = 0; k < m[1]; ++k) {
      fmap_map2_chain_t *p = chain[1] + k;
      int32_t tmp = p->qbeg;
      p->qbeg = len - p->qend; p->qend = len - tmp;
  }
  // filtering
  flag = fmap_calloc(m[0] + m[1], sizeof(char), "flag");
  fmap_sort_introsort(fmap_map2_chain, m[0] + m[1], chain[0]);
  for (k = 1; k < m[0] + m[1]; ++k) {
      fmap_map2_chain_t *p = chain[0] + k;
      for (j = 0; j < k; ++j) {
          fmap_map2_chain_t *q = chain[0] + j;
          if (flag[q->idx]) continue;
          if (q->qend >= p->qend && q->chain > p->chain * opt->seeds_rev * 2) {
              flag[p->idx] = 1;
              break;
          }
      }
  }
  for (k = 0; k < n[0] + n[1]; ++k) {
      fmap_map2_chain_t *p = z[0] + k;
      if (flag[p->chain])
        b[p->flag]->hits[p->idx].G = 0;
  }
  free(flag);
  // squeeze out filtered elements in b[2]
  for (k = 0; k < 2; ++k) {
      for (j = i = 0; j < n[k]; ++j) {
          fmap_map2_hit_t *p = b[k]->hits + j;
          if (p->G) {
              if (i != j) b[k]->hits[i++] = *p;
              else ++i;
          }
      }
      b[k]->n = i;
  }
  // free
  free(z[0]); free(chain[0]);
}
