#include <stdlib.h>
#include "fmap_bwt.h"
#include "fmap_bwt_match.h"

int
fmap_bwt_match_exact(const fmap_bwt_t *bwt, int len, const uint8_t *str, uint32_t *sa_begin, uint32_t *sa_end)
{
  uint32_t k, l, ok, ol;
  int i;
  k = 0; l = bwt->seq_len;
  for(i=0;i<len;i++) {
      uint8_t c = str[i];
      if(c > 3) return 0; // no match
      fmap_bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
      k = bwt->L2[c] + ok + 1;
      l = bwt->L2[c] + ol;
      if(k > l) break; // no match
  }
  if(k > l) return 0; // no match
  if(sa_begin) *sa_begin = k;
  if(sa_end)   *sa_end = l;
  return l - k + 1;
}

int
fmap_bwt_match_exact_alt(const fmap_bwt_t *bwt, int len, const uint8_t *str, uint32_t *k0, uint32_t *l0)
{
  int i;
  uint32_t k, l, ok, ol;
  k = *k0; l = *l0;
  for(i=0;i<len;i++) {
      uint8_t c = str[i];
      if(c > 3) return 0; // there is an N here. no match
      fmap_bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
      k = bwt->L2[c] + ok + 1;
      l = bwt->L2[c] + ol;
      if(k > l) return 0; // no match
  }
  *k0 = k; *l0 = l;
  return l - k + 1;
}
