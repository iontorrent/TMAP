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
#include "tmap_bwt_gen.h"
#include "tmap_bwt.h"
#include "tmap_bwt_match.h"
#include "tmap_bwt_aux.h"

const uint32_t tmap_bwt_aux_occ_mask[16] = {
    0xc0000000u, 0xf0000000u, 0xfc000000u, 0xff000000u,
    0xffc00000u, 0xfff00000u, 0xfffc0000u, 0xffff0000u, 
    0xffffc000u, 0xfffff000u, 0xfffffc00u, 0xffffff00u,
    0xffffffc0u, 0xfffffff0u, 0xfffffffcu, 0xffffffffu
};


const uint64_t tmap_bwt_aux_n_mask[5] = { 0xfffffffffffffffful, 0xaaaaaaaaaaaaaaaaul,
    0x5555555555555555ul, 0x0ul, 0xfffffffffffffffful };

__m128i tmap_bwt_aux_n_mask_128[3];
__m64 tmap_bwt_aux_n_mask_64[9];

void
tmap_bwt_aux_set_mask()
{
  tmap_bwt_aux_n_mask_64[0] = _mm_cvtsi64x_si64(0xfffffffffffffffful);
  tmap_bwt_aux_n_mask_64[1] = _mm_cvtsi64x_si64(0xaaaaaaaaaaaaaaaaul);
  tmap_bwt_aux_n_mask_64[2] = _mm_cvtsi64x_si64(0x5555555555555555ul);
  tmap_bwt_aux_n_mask_64[3] = _mm_setzero_si64();
  tmap_bwt_aux_n_mask_64[4] = _mm_cvtsi64x_si64(0xfffffffffffffffful);
  tmap_bwt_aux_n_mask_64[5] = _mm_cvtsi64x_si64(0x3333333333333333ul);
  tmap_bwt_aux_n_mask_64[6] = _mm_cvtsi64x_si64(0x0f0f0f0f0f0f0f0ful);
  tmap_bwt_aux_n_mask_64[7] = _mm_cvtsi64x_si64(0x1555555555555555ul);
  tmap_bwt_aux_n_mask_64[8] = _mm_cvtsi64x_si64(0x1111111111111111ul);
  tmap_bwt_aux_n_mask_128[0] = _mm_set1_epi64(tmap_bwt_aux_n_mask_64[2]);
  tmap_bwt_aux_n_mask_128[1] = _mm_set1_epi64(tmap_bwt_aux_n_mask_64[5]);
  tmap_bwt_aux_n_mask_128[2] = _mm_set1_epi64(tmap_bwt_aux_n_mask_64[6]);
}

#define occ_mask2(n) (0x5555555555555555ul - (0x1555555555555555ul >> ((n&31)<<1)))

//reduce nucleotide to bit counting - after w xor n_mask[c].
#define nucleo_5mask(w)         (w & (w >> 1) & 0x5555555555555555ul)

//transforms in five stages from bits to number of bits.
//these do not allow addition of several of the previous bit-count stage
#define nucleo_3mask(w)         ((w + (w >> 2)) & 0x3333333333333333ul)
#define nucleo_f0mask(w)        ((w + (w >> 4)) & 0x0f0f0f0f0f0f0f0ful)

// The combined three of these do the same as w * 0x101010101010101ull >> 56
#define nucleo_ffmask(w)        ((w + (w >> 8)) & 0x00ff00ff00ff00fful)
#define nucleo_4fmask(w)        ((w + (w >> 16)) & 0x0000ffff0000fffful)
#define nucleo_8fmask(w)        ((w + (w >> 32)) & 0x00000000fffffffful)

//for convenience, these macros do several stages in one call.
#define nucleo_upto5mask(p, x, w) ({            \
                                   w = (p) ^ (x);                          \
                                   nucleo_5mask(w);                        \
                                   })

#define nucleo_upto3mask(p, x, w) ({            \
                                   w = nucleo_upto5mask(p, x, w);          \
                                   nucleo_3mask(w);                        \
                                   })

#define nucleo_uptof0mask(p, x, w) ({           \
                                    w = nucleo_upto3mask(p, x, w);          \
                                    nucleo_f0mask(w);                       \
                                    })

// These combine macros are alternative stages from bits to number of bits,
// these once allow the addition of several of the previous bit-count stage
#define nucleo_combine_3mask(v, w) ({           \
                                    v = w & 0x3333333333333333ul;           \
                                    v + ((w ^ v) >> 2);                     \
                                    })

#define nucleo_combine_f0mask(v, w) ({          \
                                     v = w & 0xf0f0f0f0f0f0f0ful;            \
                                     v + ((w ^ v) >> 4);                     \
                                     })

#define nucleo_combine_ffmask(v, w) ({          \
                                     v = w & 0x00ff00ff00ff00fful;           \
                                     v + ((w ^ v) >> 8);                     \
                                     })

#define nucleo_combine_4fmask(v, w) ({          \
                                     v = w & 0x0000ffff0000fffful;           \
                                     v + ((w ^ v) >> 16);                    \
                                     })

#define nucleo_combine_8fmask(v, w) ({          \
                                     v = w & 0x00000000fffffffful;           \
                                     v + ((w ^ v) >> 32);                    \
                                     })

//nucleo_upto5mask
#define _mm_nc_combmask_xxx(x, xor, m, t1, t2) ({               \
                                                x = _mm_xor_##t1(x, xor);                               \
                                                x = _mm_and_##t1(x, _mm_srli_##t2(x, 1));               \
                                                _mm_and_##t1(x, m);                                     \
                                                })

#define _mm_nc_combmask_si64(q, x, m) _mm_nc_combmask_xxx(q, x, m, si64, si64)
#define _mm_nc_combmask_epi128(t, x, m)                                 \
  _mm_nc_combmask_xxx(t, x, m, si128, epi64)

// e.g. nucleo_3mask()
#define _mm_nmask_xxx(q, q2, shft, m, t1, t2) ({\
                                               q2 = _mm_srli_##t2(q, shft);            \
                                               q = _mm_add_##t2(q, q2);                \
                                               _mm_and_##t1(q, m);                     \
                                               })

#define _mm_nmask_si64(q, q2, shft, m)          \
  _mm_nmask_xxx(q, q2, shft, m, si64, si64)

#define _mm_nmask_epi128(q, q2, shft, m)        \
  _mm_nmask_xxx(q, q2, shft, m, si128, epi64)

// e.g. als nucleo_combine_3mask()
#define _mm_ctmap_bwt_aux_n_mask_xxx(x, y, m, shft, t1, t2) ({       \
                                                             x = _mm_and_##t1(y, m);                         \
                                                             y = _mm_xor_##t1(y, x);                         \
                                                             y = _mm_srli_##t2(y, shft);                     \
                                                             _mm_add_##t2(x, y);                             \
                                                             })

#define _mm_ctmap_bwt_aux_n_mask_si64(x, q, m, shft) _mm_ctmap_bwt_aux_n_mask_xxx(x, q, m, shft, si64, si64)
#define _mm_ctmap_bwt_aux_n_mask_epi128(t2, t, m, shft) _mm_ctmap_bwt_aux_n_mask_xxx(t2, t, m, shft, si128, epi64)

#define _mm_sum2si128_si64(t)                   \
  _mm_add_si64(_mm_movepi64_pi64(t), _mm_cvtsi64x_si64(_mm_extract_epi64(t, 1)))


// an analogy to bwt_occ() but more efficient, requiring k <= l
inline tmap_bwt_int_t 
tmap_bwt_aux_2occ(const tmap_bwt_t *bwt, tmap_bwt_int_t k, tmap_bwt_int_t *l, uint8_t c)
{
  if(TMAP_BWT_INT_MAX == k || TMAP_BWT_INT_MAX == *l) {
      *l = bwt->L2[c];
      return bwt->L2[c] + 1;
  }

  if (*l >= bwt->primary) {
      if (k > bwt->primary) {
          --k;
      } else if (k == 0) {
          *l = bwt->L2[c+1];
          k = bwt->L2[c] + 1;
          goto tmap_bwt_aux_2occ_out;
      }
      --*l;
  }
  tmap_bwt_int_t w, y, z;
  const uint32_t *p, *p2;
  tmap_bwt_int_t n = *l;

  --k;
  *l = ((const tmap_bwt_int_t *)(p = tmap_bwt_occ_intv(bwt, n)))[c] + bwt->L2[c];
  const uint64_t x = (c == 1 ? 0xaaaaaaaaaaaaaaaaul :
                      (c == 2 ? 0x5555555555555555ul :
                       (c == 3 ? 0x0ul : 0xfffffffffffffffful)));

  p += sizeof(tmap_bwt_int_t) + ((n&0x60)>>4);
  w = ((uint64_t)p[0]<<32 | p[1]) ^ x;
  y = w & (w >> 1) & occ_mask2(n);
  n = ((n^k)&~31) | ((k&0x60) >> 4);
  w = 0ul;
  z = occ_mask2(k);
  switch (n) {
    case 0x6: w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, k);
    case 0x4: w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, k);
    case 0x2: w += nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, k);
              w = nucleo_combine_3mask(k, w);
    case 0x0: z &= y;
              if (y == z) {
                  k = (tmap_bwt_int_t)(-1);
                  goto tmap_bwt_aux_2occ_out;
              }
              y = nucleo_3mask(y);
              z = nucleo_3mask(z) + w;
              k = *l;
              break;
    case 0x24:w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, k);
    case 0x22:w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, k);
              w = nucleo_combine_3mask(k, w);
    case 0x20: k = nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, k);
               y += k;
               k &= z;
               z = nucleo_3mask(k);
               y = nucleo_combine_3mask(k, y);
               if (y == z) {
                   k = (tmap_bwt_int_t)(-1);
                   goto tmap_bwt_aux_2occ_out;
               }
               z += w;
               k = *l;
               break;
    default:
               k = ((const tmap_bwt_int_t *)(p2 = tmap_bwt_occ_intv(bwt, k)))[c] + bwt->L2[c];
               n &= 0x66;
               p2 += sizeof(tmap_bwt_int_t) + (n & 0x6); //is really k
               w = ((uint64_t)p2[0]<<32 | p2[1]) ^ x;
               z &= w & (w >> 1);
               w = 0ul;
               switch (n) {
                 case 0x06: w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, n);
                 case 0x26: w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, n);
                 case 0x46: w += nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, n);
                            w = nucleo_combine_3mask(n, w);
                 case 0x66: z += nucleo_upto5mask(((uint64_t)p2[-6]<<32 | p2[-5]), x, n);
                            z += nucleo_upto5mask(((uint64_t)p2[-4]<<32 | p2[-3]), x, n);
                            z = nucleo_combine_3mask(n, z);
                            n = nucleo_upto5mask(((uint64_t)p2[-2]<<32 | p2[-1]), x, n);
                            z += nucleo_3mask(n);
                            break;
                 case 0x60: w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, n);
                 case 0x40: w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, n);
                 case 0x20: w += nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, n);
                            w = nucleo_combine_3mask(n, w);
                 case 0x00: z = nucleo_3mask(z);
                            break;
                 default:
                            switch (n) {
                              case 0x24: w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, n);
                              case 0x04: w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, n);
                              case 0x64: w += nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, n);
                                         w = nucleo_combine_3mask(n, w);
                              case 0x44: z += nucleo_upto5mask(((uint64_t)p2[-4]<<32 | p2[-3]), x, n);
                                         break;
                              case 0x42: w = nucleo_upto5mask(((uint64_t)p[-6]<<32 | p[-5]), x, n);
                              case 0x62: w += nucleo_upto5mask(((uint64_t)p[-4]<<32 | p[-3]), x, n);
                              case 0x02: w += nucleo_upto5mask(((uint64_t)p[-2]<<32 | p[-1]), x, n);
                                         w = nucleo_combine_3mask(n, w);
                            }
                            z += nucleo_upto5mask(((uint64_t)p2[-2]<<32 | p2[-1]), x, n);
                            z = nucleo_combine_3mask(n, z);
               }
               y = nucleo_3mask(y);
  }
  y += w;
  z = nucleo_combine_f0mask(w, z);
  y = nucleo_combine_f0mask(w, y);
  /*y = nucleo_ffmask(y);
    z = nucleo_ffmask(z);
    z = nucleo_4fmask(z);
    y = nucleo_4fmask(y);
   *l += nucleo_8fmask(y);
   k += nucleo_8fmask(z) + 1;*/
  k += (z * 0x101010101010101ul >> 56) + 1;
  *l += y * 0x101010101010101ul >> 56;
tmap_bwt_aux_2occ_out:
  return k;
}
