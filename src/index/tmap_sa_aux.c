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

/*
 * This is to incorporate the BWT optimizations from Roel Kluin:
 * https://github.com/RoelKluin/bwa
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_definitions.h"
#include "../io/tmap_file.h"
#include "tmap_bwt.h"
#include "tmap_bwt_aux.h"
#include "tmap_bwt_gen.h"
#include "tmap_sa.h"
#include "tmap_sa_aux.h"

static const uint64_t n_mask[5] = { 0xfffffffffffffffful, 0xaaaaaaaaaaaaaaaaul,
        0x5555555555555555ul, 0x0ul, 0xfffffffffffffffful };

static inline tmap_bwt_int_t 
tmap_sa_aux_cal_isa(const tmap_bwt_t *bwt, tmap_bwt_int_t isa)
{
  if (TMAP_LIKELY(isa != bwt->primary)) {
      tmap_bwt_int_t c, _i;
      _i = (isa < bwt->primary) ? isa : isa - 1;
      c = tmap_bwt_B0(bwt, _i);
      if (TMAP_LIKELY(isa < bwt->seq_len)) {
          uint64_t w;
          const uint32_t *p;
          w = n_mask[c];
          isa = bwt->L2[c] + ((const tmap_bwt_int_t *)(p = tmap_bwt_occ_intv(bwt, _i)))[c];
          p += sizeof(tmap_bwt_int_t) + ((_i&0x60)>>4);
          w = tmap_bwt_aux_occ_p(_i, w, (const tmap_bwt_int_t *)p);
          isa += w * 0x101010101010101ul >> 56;
      } else {
          isa = (isa == bwt->seq_len ? bwt->L2[c+1] : bwt->L2[c]);
      }
  } else {
      isa = 0;
  }

  return isa;
}

static inline tmap_bwt_int_t 
tmap_sa_aux_cal_isa_PleSl(const tmap_bwt_t *bwt, tmap_bwt_int_t isa)
{
  uint64_t w;
  const uint32_t *p;
  tmap_bwt_int_t c;
  if (TMAP_LIKELY(isa != bwt->primary)) {
      if (isa > bwt->primary) {
          if (TMAP_UNLIKELY(isa > bwt->seq_len))
            return 0;
          --isa;

      }
      c = tmap_bwt_B0(bwt, isa);
      w = n_mask[c];
      c = bwt->L2[c] + ((const tmap_bwt_int_t *)(p = tmap_bwt_occ_intv(bwt, isa)))[c];
      p += sizeof(tmap_bwt_int_t) + ((isa&0x60)>>4);
      w = tmap_bwt_aux_occ_p(isa, w, (const tmap_bwt_int_t *)p);
      isa = c + (w * 0x101010101010101ul >> 56);
  } else {
      c = tmap_bwt_B0(bwt, isa);
      if (isa == bwt->seq_len)
        ++c;
      isa = bwt->L2[c];
  }
  return isa;
}

tmap_bwt_int_t 
tmap_sa_pac_pos_aux(const tmap_sa_t *sa, const tmap_bwt_t *bwt, tmap_bwt_int_t k)
{
  tmap_bwt_int_t mask, isa = 0;
  mask = sa->sa_intv - 1;
  if (TMAP_LIKELY(bwt->primary <= bwt->seq_len)) {
      // not power of 2 before decrement
      while (k & mask) {
          ++isa;
          k = tmap_sa_aux_cal_isa_PleSl(bwt, k);
      }
  } else {
      while(k & mask) {
          ++isa;
          k = tmap_sa_aux_cal_isa(bwt, k);
      }
  }
  /* without setting sa->sa[0] = -1, the following line should be
     changed to (isa + sa->sa[k/sa->sa_intv]) % (bwt->seq_len + 1) */
  return isa + sa->sa[k/sa->sa_intv];
}
