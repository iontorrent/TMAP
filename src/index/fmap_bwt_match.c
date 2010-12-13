#include <stdlib.h>
#include <unistd.h>
#include "fmap_bwt.h"
#include "fmap_bwt_match.h"

inline void
fmap_bwt_match_occ(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, uint8_t c, fmap_bwt_match_occ_t *next)
{
  uint32_t offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the hash
      uint32_t prev_k;
      prev_k = (NULL == prev) ? 0 : prev->k;
      next->k = fmap_bwt_occ(bwt, prev_k-1, c) + bwt->L2[c] + 1;
      next->offset = offset + 1;
      next->hi = UINT32_MAX;
      next->l = UINT32_MAX; 
  }
  else { // use the hash
      uint32_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      next->offset = offset + 1;
      next->hi = (prev_hi << 2) + c;
      next->k = bwt->hash_k[next->offset-1][next->hi];
      next->l = UINT32_MAX; 
  }
}

inline void
fmap_bwt_match_2occ(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, uint8_t c, fmap_bwt_match_occ_t *next)
{
  uint32_t offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the hash
      uint32_t prev_k, prev_l;
      prev_k = (NULL == prev) ? 0 : prev->k;
      prev_l = (NULL == prev) ? bwt->seq_len : prev->l;
      fmap_bwt_2occ(bwt, prev_k-1, prev_l, c, &next->k, &next->l);
      next->offset = offset + 1;
      next->hi = UINT32_MAX;
      next->k += bwt->L2[c] + 1;
      next->l += bwt->L2[c];
  }
  else { // use the hash
      uint32_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      next->offset = offset + 1;
      next->hi = (prev_hi << 2) + c;
      next->k = bwt->hash_k[next->offset-1][next->hi];
      next->l = bwt->hash_l[next->offset-1][next->hi];
  }
}

inline void
fmap_bwt_match_occ4(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, fmap_bwt_match_occ_t next[4])
{
  uint32_t i, offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the hash
      uint32_t cntk[4];
      uint32_t prev_k;
      prev_k = (NULL == prev) ? 0 : prev->k;
      fmap_bwt_occ4(bwt, prev_k-1, cntk);
      for(i=0;i<4;i++) {
          next[i].offset = offset + 1;
          next[i].hi = UINT32_MAX;
          next[i].k = cntk[i] + bwt->L2[i] + 1;
          next[i].l = UINT32_MAX;
      }
  }
  else { // use the hash
      uint32_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      for(i=0;i<4;i++) {
          next[i].offset = offset + 1;
          next[i].hi = (prev_hi << 2) + i;
          next[i].k = bwt->hash_k[next[i].offset-1][next[i].hi];
          next[i].l = UINT32_MAX;
      }
  }
}

inline void
fmap_bwt_match_2occ4(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, fmap_bwt_match_occ_t next[4])
{
  uint32_t i, offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the hash
      uint32_t cntk[4], cntl[4];
      uint32_t prev_k, prev_l;
      prev_k = (NULL == prev) ? 0 : prev->k;
      prev_l = (NULL == prev) ? bwt->seq_len : prev->l;
      fmap_bwt_2occ4(bwt, prev_k-1, prev_l, cntk, cntl);
      for(i=0;i<4;i++) {
          next[i].offset = offset + 1;
          next[i].hi = UINT32_MAX;
          next[i].k = cntk[i] + bwt->L2[i] + 1;
          next[i].l = cntl[i] + bwt->L2[i];
      }
  }
  else { // use the hash
      uint32_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      for(i=0;i<4;i++) {
          next[i].offset = offset + 1;
          next[i].hi = (prev_hi << 2) + i;
          next[i].k = bwt->hash_k[next[i].offset-1][next[i].hi];
          next[i].l = bwt->hash_l[next[i].offset-1][next[i].hi];
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
      if(l < k || 3 < c) { // new width
          k = 0;
          l = bwt->seq_len;
          bid++;
      }
      width[i].w = l - k + 1;
      width[i].bid = bid;
  }
}

uint32_t
fmap_bwt_match_exact(const fmap_bwt_t *bwt, int len, const uint8_t *str, fmap_bwt_match_occ_t *match_sa)
{
  int32_t i;
  fmap_bwt_match_occ_t prev, next;

  prev.k = 0; prev.l = bwt->seq_len;
  prev.offset = 0;
  prev.hi = 0;

  for(i=0;i<len;i++) {
      uint8_t c = str[i];
      if(c > 3) return 0; // no match
      fmap_bwt_match_2occ(bwt, &prev, c, &next);
      prev = next;
      if(next.k > next.l) break; // no match
  }
  if(prev.k > prev.l) return 0; // no match
  if(NULL != match_sa) {
      (*match_sa) = prev;
  }
  return prev.l - prev.k + 1;
}

uint32_t
fmap_bwt_match_exact_alt(const fmap_bwt_t *bwt, int len, const uint8_t *str, fmap_bwt_match_occ_t *match_sa)
{
  int i;
  fmap_bwt_match_occ_t prev, next;

  prev = (*match_sa);

  for(i=0;i<len;i++) {
      uint8_t c = str[i];
      if(c > 3) return 0; // there is an N here. no match
      fmap_bwt_match_2occ(bwt, &prev, c, &next);
      prev = next;
      if(next.k > next.l) return 0; // no match
  }
  (*match_sa) = prev;
  return prev.l - prev.k + 1;
}
