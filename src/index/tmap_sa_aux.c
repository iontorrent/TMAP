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


static inline uint64_t 
tmap_bwt_aux_occ(tmap_bwt_int_t k, __m64 x, const uint32_t *const p)
{
  __m128i t, t2;
  t2 = t = _mm_set1_epi64(x);
  switch (k&0x60) {
    case 0x60:
      x = _mm_set_pi32(p[-6], p[-5]);
    case 0x40:
      t = _mm_set_epi64(x, _mm_set_pi32(p[-4], p[-3]));
    case 0x20: x = _mm_set_pi32(p[-2], p[-1]);
  }
  t = _mm_nc_combmask_epi128(t, t2, n_mask_128[0]);

  t2 = _mm_xor_si128(_mm_set_epi64(x, _mm_set_pi32(p[0], p[1])), t2);
  t2 = _mm_and_si128(t2, _mm_srli_epi64(t2, 1));

  x = _mm_srli_si64(n_mask_64[7], ((k&31)<<1));
  x = _mm_sub_si64(n_mask_64[2], x);
  t2 = _mm_and_si128(t2, _mm_set_epi64(n_mask_64[2], x));

  t = _mm_add_epi64(t, t2);

  t2 = _mm_srli_epi64(t, 2);
  t2 = _mm_and_si128(t2, n_mask_128[1]);
  t = _mm_and_si128(t, n_mask_128[1]);
  t = _mm_add_epi64(t, t2);

  t2 = _mm_srli_epi64(t, 4);
  t2 = _mm_and_si128(t2, n_mask_128[2]);
  t = _mm_and_si128(t, n_mask_128[2]);
  t = _mm_add_epi64(t, t2);

  return _mm_extract_epi64(t, 1) + _mm_extract_epi64(t, 0);
}

static inline tmap_bwt_int_t 
tmap_bwt_aux_invPsi(const tmap_bwt_t *bwt, tmap_bwt_int_t isa)
{
  if (TMAP_LIKELY(isa != bwt->primary)) {
      tmap_bwt_int_t c, _i;
      _i = (isa < bwt->primary) ? isa : isa - 1;
      c = bwt_B0(bwt, _i);
      if (likely(isa < bwt->seq_len)) {
          const uint32_t *p;
          isa = bwt->L2[c] + ((const tmap_bwt_int_t *)(p = tmap_bwt_occ_intv(bwt, _i)))[c];
          p += sizeof(tmap_bwt_int_t) + ((_i&0x60)>>4);
          c = tmap_bwt_aux_occ(_i, n_mask_64[c], p);
          isa += c * 0x101010101010101ul >> 56;
      } else {
          isa = (isa == bwt->seq_len ? bwt->L2[c+1] : bwt->L2[c]);
      }
  } else {
      isa = 0;
  }

  return isa;
}

tmap_bwt_int_t 
tmap_sa_pac_pos_aux(const tmap_sa_t *sa, const tmap_bwt_t *bwt, tmap_bwt_int_t k)
{
  tmap_bwt_int_t mask, sa = 0;
  mask = sa->sa_intv - 1;
  while(k & mask) {
      ++sa;
      k = tmap_bwt_aux_invPsi(bwt, k);
  }
  /* without setting bwt->sa[0] = -1, the following line should be
     changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
  return sa + sa->sa[k/sa->sa_intv];
}
