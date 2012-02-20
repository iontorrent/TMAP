/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <assert.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_vec.h"
#include "../io/tmap_file.h"
#include "tmap_bwt_gen.h"
#include "tmap_bwt.h"
#include "tmap_bwt_aux.h"
#include "tmap_bwt_check.h"
#include "tmap_bwt_match.h"
#include "tmap_bwt_smem.h"

// TODO: the primary and secondary hash (etc.)
static void
tmap_bwt_smem_extend(const tmap_bwt_t *bwt, const tmap_bwt_smem_intv_t *ik, tmap_bwt_smem_intv_t ok[4], int32_t is_back)
{
  tmap_bwt_int_t tk[4], tl[4];
  int32_t i;
  // TODO: the primary and secondary hash (etc.)
  tmap_bwt_2occ4(bwt, ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->size, tk, tl);
  for (i = 0; i != 4; ++i) {
      ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
      if(tl[i] < tk[i]) ok[i].size = 0;
      else ok[i].size = tl[i] - tk[i];
  }
  ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->size - 1 >= bwt->primary);
  ok[2].x[is_back] = ok[3].x[is_back] + ok[3].size;
  ok[1].x[is_back] = ok[2].x[is_back] + ok[2].size;
  ok[0].x[is_back] = ok[1].x[is_back] + ok[1].size;
}

static void 
tmap_bwt_smem_reverse_intvs(tmap_bwt_smem_intv_vec_t *p)
{
  if (p->n > 1) {
      int32_t j;
      for (j = 0; j < p->n>>1; ++j) {
          tmap_bwt_smem_intv_t tmp = p->a[p->n - 1 - j];
          p->a[p->n - 1 - j] = p->a[j];
          p->a[j] = tmp;
      }
  }
}

// TODO: document
int32_t
tmap_bwt_smem1(const tmap_bwt_t *bwt, int32_t len, const uint8_t *q, int32_t x, tmap_bwt_smem_intv_vec_t *mem, tmap_bwt_smem_intv_vec_t *tmpvec[2])
{
  int32_t i, j, c, ret;
  tmap_bwt_smem_intv_t ik, ok[4];
  tmap_bwt_smem_intv_vec_t a[2], *prev, *curr, *swap;

  mem->n = 0;
  if (q[x] > 3) return x + 1;
  tmap_vec_init(a[0]); tmap_vec_init(a[1]);
  prev = tmpvec[0]? tmpvec[0] : &a[0];
  curr = tmpvec[1]? tmpvec[1] : &a[1];
  tmap_bwt_smem_set_intv(bwt, q[x], ik);
  ik.info = x + 1;

  for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
      if (q[i] < 4) {
          c = 3 - q[i];
          tmap_bwt_smem_extend(bwt, &ik, ok, 0);
          if (ok[c].size != ik.size) // change of the interval size
            tmap_vec_push(tmap_bwt_smem_intv_t, *curr, ik);
          if (ok[c].size <= 0) {
              ok[c].size = 0;
              break; // cannot be extended
          }
          ik = ok[c]; ik.info = i + 1;
      } else { // an ambiguous base
          tmap_vec_push(tmap_bwt_smem_intv_t, *curr, ik);
          break; // cannot be extended; in this case, i<len always stands
      }
  }
  if (i == len) tmap_vec_push(tmap_bwt_smem_intv_t, *curr, ik); // push the last interval if we reach the end
  tmap_bwt_smem_reverse_intvs(curr); // s.t. smaller intervals visited first
  ret = curr->a[0].info; // this will be the returned value
  swap = curr; curr = prev; prev = swap;

  for (i = x - 1; i >= -1; --i) { // backward search for MEMs
      c = (i < 0) ? 0 : q[i];
      if (c > 3) break;
      for (j = 0, curr->n = 0; j < prev->n; ++j) {
          tmap_bwt_smem_intv_t *p = &prev->a[j];
          tmap_bwt_smem_extend(bwt, p, ok, 1);
          if (ok[c].size <= 0 || i == -1) { // keep the hit if reaching the beginning or not extended further
              if (curr->n == 0) { // curr->n to make sure there is no longer matches
                  if (mem->n == 0 || i + 1 < mem->a[mem->n-1].info>>32) { // skip contained matches
                      ik = *p; ik.info |= (uint64_t)(i + 1)<<32;
                      tmap_vec_push(tmap_bwt_smem_intv_t, *mem, ik);
                  }
              } // otherwise the match is contained in another longer match
          }
          if (0 < ok[c].size && (curr->n == 0 || ok[c].size != curr->a[curr->n-1].size)) {
              ok[c].info = p->info;
              tmap_vec_push(tmap_bwt_smem_intv_t, *curr, ok[c]);
          }
      }
      if (curr->n == 0) break;
      swap = curr; curr = prev; prev = swap;
  }
  tmap_bwt_smem_reverse_intvs(mem); // s.t. sorted by the start coordinate

  if (tmpvec[0] == 0) free(a[0].a);
  if (tmpvec[1] == 0) free(a[1].a);
  return ret;
}
