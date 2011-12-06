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

static const uint32_t occ_mask[16] = {
    0xc0000000u, 0xf0000000u, 0xfc000000u, 0xff000000u,
    0xffc00000u, 0xfff00000u, 0xfffc0000u, 0xffff0000u, 
    0xffffc000u, 0xfffff000u, 0xfffffc00u, 0xffffff00u,
    0xffffffc0u, 0xfffffff0u, 0xfffffffcu, 0xffffffffu
};

static const uint64_t occ_mask2[32] = {
    0x40000000ul, 	0x50000000ul, 	0x54000000ul,
    0x55000000ul, 	0x55400000ul, 	0x55500000ul,
    0x55540000ul, 	0x55550000ul, 	0x55554000ul,
    0x55555000ul, 	0x55555400ul, 	0x55555500ul,
    0x55555540ul, 	0x55555550ul, 	0x55555554ul,
    0x55555555ul, 	0x4000000055555555ul, 	0x5000000055555555ul,
    0x5400000055555555ul, 	0x5500000055555555ul, 	0x5540000055555555ul,
    0x5550000055555555ul, 	0x5554000055555555ul, 	0x5555000055555555ul,
    0x5555400055555555ul, 	0x5555500055555555ul, 	0x5555540055555555ul,
    0x5555550055555555ul, 	0x5555554055555555ul, 	0x5555555055555555ul,
    0x5555555455555555ul, 0x5555555555555555ul
};

#define nucleo_5mask(w) 	(w & (w >> 1) & 0x5555555555555555ul)
#define nucleo_3mask(w)		((w + (w >> 2)) & 0x3333333333333333ul)
#define nucleo_f0mask(w)	((w + (w >> 4)) & 0x0f0f0f0f0f0f0f0ful)

#define nucleo_upto5mask(p, x, w) ({		\
                                   w = *(p) ^ (x);				\
                                   nucleo_5mask(w);			\
                                   })

#define nucleo_upto3mask(p, x, w) ({		\
                                   w = nucleo_upto5mask(p, x, w);		\
                                   nucleo_3mask(w);			\
                                   })

#define nucleo_uptof0mask(p, x, w) ({		\
                                    w = nucleo_upto3mask(p, x, w);		\
                                    nucleo_f0mask(w);			\
                                    })


#define nucleo_pre_combine_3mask(v, w) ({	\
                                        v = w;					\
                                        w &= 0x3333333333333333ul;		\
                                        })

// doesn't work?
#define nucleo_pre_combine_f0mask(v, w) ({	\
                                         v = w;					\
                                         w &= 0x0f0f0f0f0f0f0f0ful;		\
                                         w ^ v;					\
                                         })

#define nucleo_combine_3mask(v, w) ({		\
                                    v = w & 0x3333333333333333ul;		\
                                    v + ((w ^ v) >> 2);			\
                                    })

#define nucleo_combine_f0mask(v, w) ({		\
                                     v = w & 0xf0f0f0f0f0f0f0ful;		\
                                     v + ((w ^ v) >> 4);			\
                                     })

inline uint64_t 
tmap_bwt_aux_occ_p(const tmap_bwt_int_t k, uint64_t x, const tmap_bwt_int_t *const p)
{
  uint64_t w = 0ul;
  uint64_t y = *p ^ x;
  y &= (y >> 1) & occ_mask2[k&31];
  switch (k&0x60) {
    case 0x60: y += nucleo_upto5mask(p - 3, x, w);
               y += nucleo_upto5mask(p - 2, x, w);
               y = nucleo_combine_3mask(w, y);
               y += nucleo_upto3mask(p - 1, x, w);
               return nucleo_combine_f0mask(w, y);
    case 0x00: y = nucleo_3mask(y);
               break;
    case 0x40: y += nucleo_upto5mask(p - 2, x, w);
    case 0x20: y += nucleo_upto5mask(p - 1, x, w);
               y = nucleo_combine_3mask(w, y);
  }
  return nucleo_f0mask(y);
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
inline tmap_bwt_int_t 
tmap_bwt_aux_2occ(const tmap_bwt_t *bwt, tmap_bwt_int_t k, tmap_bwt_int_t *l, uint8_t c)
{
  if(TMAP_BWT_INT_MAX == k || TMAP_BWT_INT_MAX == *l) {
      if(*l == TMAP_BWT_INT_MAX) k = TMAP_BWT_INT_MAX;
      *l = tmap_bwt_occ(bwt, *l, c) + bwt->L2[c];
      k = tmap_bwt_occ(bwt, k, c) + bwt->L2[c] + 1;
      goto tmap_bwt_aux_2occ_out;
  }

  if (*l >= bwt->primary) {
      if (k > bwt->primary) {
          --k;
      } 
      else if (k == 0) {
          *l = tmap_bwt_occ(bwt, *l, c) + bwt->L2[c];
          //*l = bwt->L2[c+1];
          k = bwt->L2[c] + 1;
          goto tmap_bwt_aux_2occ_out;
      }
      --*l;
  }

  tmap_bwt_int_t w, y, z;
  const tmap_bwt_int_t *p, *p2;
  tmap_bwt_int_t n = *l;

  --k;
  p = (const tmap_bwt_int_t*)tmap_bwt_occ_intv(bwt, n);
  *l = p[c] + bwt->L2[c];
  const uint64_t x = (c == 1 ? 0xaaaaaaaaaaaaaaaaul :
                      (c == 2 ? 0x5555555555555555ul :
                       (c == 3 ? 0x0ul : 0xfffffffffffffffful)));

  p += (sizeof(tmap_bwt_int_t)/2) + ((n&0x60)>>5);
  w = *p ^ x;
  y = w & (w >> 1) & occ_mask2[n&31];
  n = ((n^k)&~31) | ((k&0x60) >> 5);
  w = 0ul;
  z = occ_mask2[k&31];
  switch (n) {
    case 0x3: w = nucleo_upto5mask(p - 3, x, k);
    case 0x2: w += nucleo_upto5mask(p - 2, x, k);
    case 0x1: w += nucleo_upto5mask(p - 1, x, k);
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
    case 0x22:w = nucleo_upto5mask(p - 3, x, k);
    case 0x21:w += nucleo_upto5mask(p - 2, x, k);
              w = nucleo_combine_3mask(k, w);
    case 0x20: k = nucleo_upto5mask(p - 1, x, k);
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
               p2 = (const tmap_bwt_int_t *)tmap_bwt_occ_intv(bwt, k);
               k = p2[c] + bwt->L2[c];
               n &= 0x63;
               p2 += (sizeof(tmap_bwt_int_t)/2) + (n & 0x3); //is really k
               w = *p2 ^ x;
               z &= w & (w >> 1);
               w = 0ul;
               switch (n) {
                 case 0x03: w = nucleo_upto5mask(p - 3, x, n);
                 case 0x23: w += nucleo_upto5mask(p - 2, x, n);
                 case 0x43: w += nucleo_upto5mask(p - 1, x, n);
                            w = nucleo_combine_3mask(n, w);
                 case 0x63: z += nucleo_upto5mask(p2 - 3, x, n);
                            z += nucleo_upto5mask(p2 - 2, x, n);
                            z = nucleo_combine_3mask(n, z);
                            n = nucleo_upto5mask(p2 - 1, x, n);
                            z += nucleo_3mask(n);
                            break;
                 case 0x60: w = nucleo_upto5mask(p - 3, x, n);
                 case 0x40: w += nucleo_upto5mask(p - 2, x, n);
                 case 0x20: w += nucleo_upto5mask(p - 1, x, n);
                            w = nucleo_combine_3mask(n, w);
                 case 0x00: z = nucleo_3mask(z);
                            break;
                 default:
                            switch (n) {
                              case 0x22: w = nucleo_upto5mask(p - 3, x, n);
                              case 0x02: w += nucleo_upto5mask(p - 2, x, n);
                              case 0x62: w += nucleo_upto5mask(p - 1, x, n);
                                         w = nucleo_combine_3mask(n, w);
                              case 0x42: z += nucleo_upto5mask(p2 - 2, x, n);
                                         break;
                              case 0x41: w = nucleo_upto5mask(p - 3, x, n);
                              case 0x61: w += nucleo_upto5mask(p - 2, x, n);
                              case 0x01: w += nucleo_upto5mask(p - 1, x, n);
                                         w = nucleo_combine_3mask(n, w);
                            }
                            z += nucleo_upto5mask(p2 - 1, x, n);
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
