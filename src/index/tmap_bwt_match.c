/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <unistd.h>
#include "tmap_bwt.h"
#include "tmap_bwt_match.h"

inline void
tmap_bwt_match_occ(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, uint8_t c, tmap_bwt_match_occ_t *next)
{
  uint32_t offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the hash
      tmap_bwt_int_t prev_k;
      prev_k = (NULL == prev) ? 0 : prev->k;
      next->k = tmap_bwt_occ(bwt, prev_k-1, c) + bwt->L2[c] + 1;
      next->offset = offset + 1;
      next->hi = TMAP_BWT_INT_MAX;
      next->l = TMAP_BWT_INT_MAX; 
  }
  else { // use the hash
      uint64_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      next->offset = offset + 1;
      next->hi = (prev_hi << 2) + c;
      next->k = bwt->hash_k[next->offset-1][next->hi];
      next->l = TMAP_BWT_INT_MAX; 
  }
}

inline void
tmap_bwt_match_2occ(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, uint8_t c, tmap_bwt_match_occ_t *next)
{
  uint32_t offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the hash
      tmap_bwt_int_t prev_k, prev_l;
      prev_k = (NULL == prev) ? 0 : prev->k;
      prev_l = (NULL == prev) ? bwt->seq_len : prev->l;
      tmap_bwt_2occ(bwt, prev_k-1, prev_l, c, &next->k, &next->l);
      next->offset = offset + 1;
      next->hi = TMAP_BWT_INT_MAX;
      next->k += bwt->L2[c] + 1;
      next->l += bwt->L2[c];
  }
  else { // use the hash
      uint64_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      next->offset = offset + 1;
      next->hi = (prev_hi << 2) + c;
      next->k = bwt->hash_k[next->offset-1][next->hi];
      next->l = bwt->hash_l[next->offset-1][next->hi];
  }
}

inline void
tmap_bwt_match_occ4(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, tmap_bwt_match_occ_t next[4])
{
  uint32_t i, offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the hash
      tmap_bwt_int_t cntk[4];
      tmap_bwt_int_t prev_k;
      prev_k = (NULL == prev) ? 0 : prev->k;
      tmap_bwt_occ4(bwt, prev_k-1, cntk);
      for(i=0;i<4;i++) {
          next[i].offset = offset + 1;
          next[i].hi = TMAP_BWT_INT_MAX;
          next[i].k = cntk[i] + bwt->L2[i] + 1;
          next[i].l = TMAP_BWT_INT_MAX;
      }
  }
  else { // use the hash
      uint64_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      for(i=0;i<4;i++) {
          next[i].offset = offset + 1;
          next[i].hi = (prev_hi << 2) + i;
          next[i].k = bwt->hash_k[next[i].offset-1][next[i].hi];
          next[i].l = TMAP_BWT_INT_MAX;
      }
  }
}

inline void
tmap_bwt_match_2occ4(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, tmap_bwt_match_occ_t next[4])
{
  uint32_t i, offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the hash
      tmap_bwt_int_t cntk[4], cntl[4];
      tmap_bwt_int_t prev_k, prev_l;
      prev_k = (NULL == prev) ? 0 : prev->k;
      prev_l = (NULL == prev) ? bwt->seq_len : prev->l;
      tmap_bwt_2occ4(bwt, prev_k-1, prev_l, cntk, cntl);
      for(i=0;i<4;i++) {
          next[i].offset = offset + 1;
          next[i].hi = TMAP_BWT_INT_MAX;
          next[i].k = cntk[i] + bwt->L2[i] + 1;
          next[i].l = cntl[i] + bwt->L2[i];
      }
  }
  else { // use the hash
      uint64_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      for(i=0;i<4;i++) {
          next[i].offset = offset + 1;
          next[i].hi = (prev_hi << 2) + i;
          next[i].k = bwt->hash_k[next[i].offset-1][next[i].hi];
          next[i].l = bwt->hash_l[next[i].offset-1][next[i].hi];
      }
  }
}

void
tmap_bwt_match_cal_width_forward(const tmap_bwt_t *bwt, int len, const char *str, tmap_bwt_match_width_t *width)
{
  // 'width[i]' is the lower bound of the number of differences in str[i,len]
  tmap_bwt_int_t k, l, ok, ol;
  int32_t i, bid;

  bid = 0;
  k = 0; l = bwt->seq_len;
  for(i=len-1;0<=i;i--) {
      uint8_t c = (int)str[i];
      if(c < 4) {
          tmap_bwt_2occ(bwt, k-1, l, c, &ok, &ol);
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

void
tmap_bwt_match_cal_width_reverse(const tmap_bwt_t *bwt, int len, const char *str, tmap_bwt_match_width_t *width)
{
  // 'width[i]' is the lower bound of the number of differences in str[0,i]
  tmap_bwt_int_t k, l, ok, ol;
  int32_t i, bid;

  bid = 0;
  k = 0; l = bwt->seq_len;
  for(i=0;i<len;i++) {
      uint8_t c = (int)str[i];
      if(c < 4) {
          tmap_bwt_2occ(bwt, k-1, l, c, &ok, &ol);
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
  width[len].w = 0;
  width[len].bid = ++bid;
}

tmap_bwt_int_t
tmap_bwt_match_exact(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa)
{
  int32_t i;
  uint8_t c = 0;
  tmap_bwt_match_occ_t prev, next;

  prev.k = 0; prev.l = bwt->seq_len;
  prev.offset = 0;
  prev.hi = 0;

  for(i=0;i<len;i++) {
      c = str[i];
      if(TMAP_UNLIKELY(3 < c)) { 
          prev.offset++;
          prev.k = prev.l + 1;
          break;
      }
      tmap_bwt_match_2occ(bwt, &prev, c, &next);
      prev = next;
      if(TMAP_BWT_INT_MAX == next.k || next.k > next.l) break; // no match
  }
  if(NULL != match_sa) {
      (*match_sa) = prev;
  }
  if(TMAP_BWT_INT_MAX == prev.k || TMAP_UNLIKELY(3 < c) || prev.k > prev.l) return 0; // no match
  return prev.l - prev.k + 1;
}

tmap_bwt_int_t
tmap_bwt_match_exact_reverse(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa)
{
  int32_t i;
  uint8_t c = 0;
  tmap_bwt_match_occ_t prev, next;

  prev.k = 0; prev.l = bwt->seq_len;
  prev.offset = 0;
  prev.hi = 0;

  for(i=len-1;0<=i;i--) {
      c = str[i];
      if(TMAP_UNLIKELY(3 < c)) { 
          prev.offset++; 
          prev.k = prev.l + 1;
          break;
      }
      tmap_bwt_match_2occ(bwt, &prev, c, &next);
      prev = next;
      if(TMAP_BWT_INT_MAX == next.k || next.k > next.l) break; // no match
  }
  if(NULL != match_sa) {
      (*match_sa) = prev;
  }
  if(TMAP_BWT_INT_MAX == prev.k || TMAP_UNLIKELY(3 < c) || prev.k > prev.l) return 0; // no match
  return prev.l - prev.k + 1;
}

tmap_bwt_int_t
tmap_bwt_match_exact_alt(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa)
{
  int i;
  tmap_bwt_match_occ_t next;

  for(i=0;i<len;i++) {
      uint8_t c = str[i];
      if(TMAP_UNLIKELY(c > 3)) return 0; // there is an N here. no match
      tmap_bwt_match_2occ(bwt, match_sa, c, &next);
      (*match_sa) = next;
      if(TMAP_BWT_INT_MAX == next.k || next.k > next.l) return 0; // no match
  }
  return match_sa->l - match_sa->k + 1;
}

tmap_bwt_int_t
tmap_bwt_match_exact_alt_reverse(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa)
{
  int i;
  tmap_bwt_match_occ_t next;

  for(i=len-1;0<=i;i--) {
      uint8_t c = str[i];
      if(TMAP_UNLIKELY(c > 3)) return 0; // there is an N here. no match
      tmap_bwt_match_2occ(bwt, match_sa, c, &next);
      (*match_sa) = next;
      if(TMAP_BWT_INT_MAX == next.k || next.k > next.l) return 0; // no match
  }
  return match_sa->l - match_sa->k + 1;
}
