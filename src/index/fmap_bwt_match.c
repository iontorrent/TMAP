#include <stdlib.h>
#include "fmap_bwt.h"
#include "fmap_bwt_match.h"

inline void
fmap_bwt_match_occ(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, uint8_t c, fmap_bwt_match_occ_t *next)
{
  if(bwt->hash_width <= prev->offset) { // do not use the hash
      next->offset = prev->offset + 1;
      next->hi = UINT32_MAX;
      next->k = fmap_bwt_occ(bwt, prev->k, c);
      next->l = UINT32_MAX; 
  }
  else { // use the hash
      next->offset = prev->offset + 1;
      next->hi = (prev->hi << 2) + c;
      next->k = bwt->hash_k[next->offset-1][next->hi];
      next->l = UINT32_MAX; 
  }
}

inline void
fmap_bwt_match_2occ(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, uint8_t c, fmap_bwt_match_occ_t *next)
{
  if(bwt->hash_width <= prev->offset) { // do not use the hash
      next->offset = prev->offset + 1;
      next->hi = UINT32_MAX;
      fmap_bwt_2occ(bwt, prev->k, prev->l, c, &next->k, &next->l);
  }
  else { // use the hash
      next->offset = prev->offset + 1;
      next->hi = (prev->hi << 2) + c;
      next->k = bwt->hash_k[next->offset-1][next->hi];
      next->l = bwt->hash_l[next->offset-1][next->hi];
  }
}

inline void
fmap_bwt_match_occ4(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, fmap_bwt_match_occ_t *next[4])
{
  uint32_t i;
  if(bwt->hash_width <= prev->offset) { // do not use the hash
      uint32_t cntk[4];
      fmap_bwt_occ4(bwt, prev->k, cntk);
      for(i=0;i<4;i++) {
          next[i]->offset = prev->offset + 1;
          next[i]->hi = UINT32_MAX;
          next[i]->k = cntk[i];
          next[i]->l = UINT32_MAX;
      }
  }
  else { // use the hash
      for(i=0;i<4;i++) {
          next[i]->offset = prev->offset + 1;
          next[i]->hi = (prev->hi << 2) + i; // the next hash index
          next[i]->k = bwt->hash_k[next[i]->offset-1][next[i]->hi];
          next[i]->l = UINT32_MAX;
      }
  }
}

inline void
fmap_bwt_match_2occ4(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, fmap_bwt_match_occ_t *next[4])
{
  uint32_t i;
  if(bwt->hash_width <= prev->offset) { // do not use the hash
      uint32_t cntk[4], cntl[4];
      fmap_bwt_2occ4(bwt, prev->k, prev->l, cntk, cntl);
      for(i=0;i<4;i++) {
          next[i]->offset = prev->offset + 1;
          next[i]->hi = UINT32_MAX;
          next[i]->k = cntk[i];
          next[i]->l = cntl[i];
      }
  }
  else { // use the hash
      for(i=0;i<4;i++) {
          next[i]->offset = prev->offset + 1;
          next[i]->hi = (prev->hi << 2) + i; // the next hash index
          next[i]->k = bwt->hash_k[next[i]->offset-1][next[i]->hi];
          next[i]->l = bwt->hash_l[next[i]->offset-1][next[i]->hi];
      }
  }
}

void
fmap_bwt_match_cal_width(const fmap_bwt_t *bwt, int len, const char *str, fmap_bwt_match_width_t *width)
{
  // 'width[i]' is the lower bound of the number of differences in str[i,len]
  uint32_t k, l, ok, ol;
  int i, bid;

  bid = 0;
  k = 0; l = bwt->seq_len;
  for(i=len-1;0<=i;i--) {
      uint8_t c = (int)str[i];
      if(c < 4) {
          fmap_bwt_2occ(bwt, k-1, l, c, &ok, &ol);
          k = bwt->L2[c] + ok + 1;
          l = bwt->L2[c] + ol;
      }
      if(l < k || 3 < c) { // new widith
          k = 0;
          l = bwt->seq_len;
          bid++;
      }
      width[i].w = l - k + 1;
      width[i].bid = bid;
  }
}

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
