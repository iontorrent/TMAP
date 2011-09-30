#include <stdlib.h>
#include <stdint.h>
#include "tmap_error.h"
#include "tmap_alloc.h"
#include "tmap_rand.h"

tmap_rand_t*
tmap_rand_init(uint32_t seed)
{
  int32_t i;
  tmap_rand_t *r = NULL;

  r = tmap_calloc(1, sizeof(tmap_rand_t), "r"); 
  r->Q[0] = seed;
  r->Q[1] = seed + TMAP_RAND_PHI;
  r->Q[2] = seed + TMAP_RAND_PHI + TMAP_RAND_PHI;
  for(i=3;i<TMAP_RAND_R_LAG;i++) {
    r->Q[i] = r->Q[i-3] ^ r->Q[i-2] ^ TMAP_RAND_PHI ^ i;
  }
  r->i = TMAP_RAND_PHI - 1;
  r->c = TMAP_RAND_C;

  return r;
}

double
tmap_rand_get(tmap_rand_t *r)
{
  uint64_t t;
  uint32_t x;
  r->i = (r->i + 1) & (TMAP_RAND_R_LAG-1);
  t = TMAP_RAND_A * r->Q[r->i] + r->c;
  r->c = (t >> 32);
  x = t + r->c;
  if (x < r->c) {
      x++;
      r->c++;
  }
  return ((double)(r->Q[r->i] = TMAP_RAND_R - x) / UINT32_MAX);
}

void
tmap_rand_destroy(tmap_rand_t *r)
{
  free(r);
}
