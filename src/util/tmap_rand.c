#include <stdlib.h>
#include <stdint.h>
#include "tmap_error.h"
#include "tmap_alloc.h"
#include "tmap_rand.h"

/**
 * 64-bit Mersenne Twister pseudorandom number generator. Adapted from:
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
 * which was written by Takuji Nishimura and Makoto Matsumoto and released
 * under the 3-clause BSD license.
*/

static void 
tmap_rand_srand0(uint64_t seed, tmap_rand_t *r)
{
  r->mt[0] = seed;
  for (r->mti = 1; r->mti < TMAP_RAND_NN; ++r->mti)
    r->mt[r->mti] = 6364136223846793005ULL * (r->mt[r->mti - 1] ^ (r->mt[r->mti - 1] >> 62)) + r->mti;
}

tmap_rand_t *
tmap_rand_init(uint64_t seed)
{
  tmap_rand_t *r;
  r = tmap_calloc(1, sizeof(tmap_rand_t), "r"); 
  tmap_rand_srand0(seed, r);
  return r;
}

uint64_t 
tmap_rand_int(tmap_rand_t *r)
{
  uint64_t x;
  static const uint64_t mag01[2] = { 0, 0xB5026F5AA96619E9ULL };
  if (r->mti >= TMAP_RAND_NN) {
      int i;
      if (r->mti == TMAP_RAND_NN + 1) tmap_rand_srand0(5489ULL, r);
      for (i = 0; i < TMAP_RAND_NN - TMAP_RAND_MM; ++i) {
          x = (r->mt[i] & TMAP_RAND_UM) | (r->mt[i+1] & TMAP_RAND_LM);
          r->mt[i] = r->mt[i + TMAP_RAND_MM] ^ (x>>1) ^ mag01[(int)(x&1)];
      }
      for (; i < TMAP_RAND_NN - 1; ++i) {
          x = (r->mt[i] & TMAP_RAND_UM) | (r->mt[i+1] & TMAP_RAND_LM);
          r->mt[i] = r->mt[i + (TMAP_RAND_MM - TMAP_RAND_NN)] ^ (x>>1) ^ mag01[(int)(x&1)];
      }
      x = (r->mt[TMAP_RAND_NN - 1] & TMAP_RAND_UM) | (r->mt[0] & TMAP_RAND_LM);
      r->mt[TMAP_RAND_NN - 1] = r->mt[TMAP_RAND_MM - 1] ^ (x>>1) ^ mag01[(int)(x&1)];
      r->mti = 0;
  }
  x = r->mt[r->mti++];
  x ^= (x >> 29) & 0x5555555555555555ULL;
  x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
  x ^= (x << 37) & 0xFFF7EEE000000000ULL;
  x ^= (x >> 43);
  return x;
}

double
tmap_rand_get(tmap_rand_t *r)
{
  return tmap_rand_int(r) / (double)UINT64_MAX;
}

void
tmap_rand_destroy(tmap_rand_t *r)
{
  free(r);
}
