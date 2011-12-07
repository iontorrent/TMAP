/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "tmap_bwt.h"
#include "../util/tmap_hash.h"
#include "../util/tmap_definitions.h"
#include "tmap_bwt_match.h"
#include "tmap_bwt_match_hash.h"

// TODO: use the register keyword on pointers

// TODO: change this if the size of the pacpos
TMAP_HASH_MAP_INIT_INT(tmap_bwt_match_hash, uint32_t) 

tmap_bwt_match_hash_t*
tmap_bwt_match_hash_init()
{
  int32_t i;
  tmap_bwt_match_hash_t *h = NULL;
  h = tmap_calloc(sizeof(tmap_bwt_match_hash_t), 1, "hash");
  for(i=0;i<4;i++) {
      h->hash[i] = (void*)tmap_hash_init(tmap_bwt_match_hash);
  }
  return h;
}

void
tmap_bwt_match_hash_destroy(tmap_bwt_match_hash_t *h)
{
  int32_t i;
  for(i=0;i<4;i++) {
      tmap_hash_t(tmap_bwt_match_hash) *hash = (tmap_hash_t(tmap_bwt_match_hash)*)(h->hash[i]);
      tmap_hash_destroy(tmap_bwt_match_hash, hash); 
  }
  free(h);
}

void
tmap_bwt_match_hash_clear(tmap_bwt_match_hash_t *h)
{
  int32_t i;
  for(i=0;i<4;i++) {
      tmap_hash_t(tmap_bwt_match_hash) *hash = (tmap_hash_t(tmap_bwt_match_hash)*)(h->hash[i]);
      tmap_hash_clear(tmap_bwt_match_hash, hash); 
  }
}

int32_t
tmap_bwt_match_hash_put(tmap_bwt_match_hash_t *h, uint32_t key, uint8_t c, uint32_t val)
{
  int32_t ret;
  tmap_hash_int_t iter;
  tmap_hash_t(tmap_bwt_match_hash) *hash = (tmap_hash_t(tmap_bwt_match_hash)*)(h->hash[c]);
  // get the iter
  iter = tmap_hash_put(tmap_bwt_match_hash, hash, key, &ret); 
  // set the value
  tmap_hash_val(hash, iter) = val;
  return ret;
}

uint32_t
tmap_bwt_match_hash_get(tmap_bwt_match_hash_t *h, uint32_t key, uint8_t c, uint32_t *found)
{
  tmap_hash_int_t iter;
  tmap_hash_t(tmap_bwt_match_hash) *hash = (tmap_hash_t(tmap_bwt_match_hash)*)(h->hash[c]);
  // get the iter
  iter = tmap_hash_get(tmap_bwt_match_hash, hash, key); 
  // check it was found
  if(tmap_hash_end(hash) == iter) { // not found
      *found = 0;
      return UINT32_MAX;
  }
  else {
      *found = 1;
      // return the value
      return (uint32_t)(tmap_hash_val(hash, iter));
  }
}

// NB: assumes we pass in prev_k, not prev_k-1, and that val had 'L2[c] + 1' added
static inline int32_t
tmap_bwt_match_hash_put_k(tmap_bwt_match_hash_t *h, uint32_t k, uint8_t c, uint32_t val)
{
  return tmap_bwt_match_hash_put(h, k-1, c, val-1);
}

// NB: assumes we pass in prev_l, and that val had 'L2[c]' added
static inline int32_t
tmap_bwt_match_hash_put_l(tmap_bwt_match_hash_t *h, uint32_t l, uint8_t c, uint32_t val)
{
  return tmap_bwt_match_hash_put(h, l, c, val);
}

// NB: assumes we pass in prev_k, not prev_k-1, and that the stored val had 'L2[c]' added
static inline uint32_t
tmap_bwt_match_hash_get_k(tmap_bwt_match_hash_t *h, uint32_t k, uint8_t c, uint32_t *found)
{
  return tmap_bwt_match_hash_get(h, k-1, c, found) + 1;
}

// NB: assumes we pass in prev_l, and that the stored val had 'L2[c]' added
static inline uint32_t
tmap_bwt_match_hash_get_l(tmap_bwt_match_hash_t *h, uint32_t l, uint8_t c, uint32_t *found)
{
  return tmap_bwt_match_hash_get(h, l, c, found);
}

inline void
tmap_bwt_match_hash_occ(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, uint8_t c, 
                          tmap_bwt_match_occ_t *next, tmap_bwt_match_hash_t *hash)
{
  uint32_t offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the bwt hash
      uint32_t prev_k, found_k;
      prev_k = (NULL == prev) ? 0 : prev->k;
      if(NULL == hash) {
          next->k = tmap_bwt_occ(bwt, prev_k-1, c) + bwt->L2[c] + 1;
      }
      else { // test the user "hash"
          next->k = tmap_bwt_match_hash_get_k(hash, prev_k, c, &found_k);
          // compute the k, value, if necessary
          if(0 == found_k) {
              next->k = tmap_bwt_occ(bwt, prev_k-1, c) + bwt->L2[c] + 1;
              // put it in the hash
              tmap_bwt_match_hash_put_k(hash, prev_k, c, next->k);
          }
      }
      next->offset = offset + 1;
      next->hi = UINT32_MAX;
      next->l = UINT32_MAX; 
  }
  else { // use the bwt hash
      uint64_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      next->offset = offset + 1;
      next->hi = (prev_hi << 2) + c;
      next->k = bwt->hash_k[next->offset-1][next->hi];
      next->l = UINT32_MAX; 
  }
}

inline void
tmap_bwt_match_hash_2occ(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, uint8_t c, 
                    tmap_bwt_match_occ_t *next, tmap_bwt_match_hash_t *hash)
{
  uint32_t offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the bwt hash
      uint32_t prev_k, prev_l, found_k, found_l;
      prev_k = (NULL == prev) ? 0 : prev->k;
      prev_l = (NULL == prev) ? bwt->seq_len : prev->l;
      if(NULL == hash) {
          tmap_bwt_2occ(bwt, prev_k-1, prev_l, c, &next->k, &next->l); 
          next->k += bwt->L2[c] + 1;
          next->l += bwt->L2[c];
      }
      else { // test the user "hash"
          next->k = tmap_bwt_match_hash_get_k(hash, prev_k, c, &found_k); // compute k
          // compute the k value, if necessary
          if(0 == found_k) { 
              tmap_bwt_2occ(bwt, prev_k-1, prev_l, c, &next->k, &next->l); 
              next->k += bwt->L2[c] + 1;
              next->l += bwt->L2[c];
              // put it in the hash
              tmap_bwt_match_hash_put_k(hash, prev_k, c, next->k);
              tmap_bwt_match_hash_put_l(hash, prev_l, c, next->l);
          }
          else {
              next->l = tmap_bwt_match_hash_get_l(hash, prev_l, c, &found_l); // compute -l
              // compute the l value, if necessary
              if(0 == found_l) { // for l
                  next->l = tmap_bwt_occ(bwt, prev_l, c) + bwt->L2[c];
                  // put it in the hash
                  tmap_bwt_match_hash_put_l(hash, prev_l, c, next->l);
              }
          }
      }
      next->offset = offset + 1;
      next->hi = UINT32_MAX;
  }
  else { // use the bwt hash
      uint64_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      next->offset = offset + 1;
      next->hi = (prev_hi << 2) + c;
      next->k = bwt->hash_k[next->offset-1][next->hi];
      next->l = bwt->hash_l[next->offset-1][next->hi];
  }
}

inline void
tmap_bwt_match_hash_occ4(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, 
                    tmap_bwt_match_occ_t next[4], tmap_bwt_match_hash_t *hash)
{
  uint32_t i, offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the bwt hash
      uint32_t cntk[4];
      uint32_t prev_k, found_k;
      prev_k = (NULL == prev) ? 0 : prev->k;
      if(NULL == hash) {
          tmap_bwt_occ4(bwt, prev_k-1, cntk);
          for(i=0;i<4;i++) {
              next[i].offset = offset + 1;
              next[i].hi = UINT32_MAX;
              next[i].k = cntk[i] + bwt->L2[i] + 1;
              next[i].l = UINT32_MAX;
          }
      }
      else { // test the user "hash" for k
          for(i=0;i<4;i++) {
              next[i].k = tmap_bwt_match_hash_get_k(hash, prev_k, i, &found_k); // compute k
              if(0 == found_k) { // for k
                  break;
              }
          } // TODO: should we use tmap_bwt_occ if we break out of this loop?
          // compute the k value, if necessary
          if(0 == found_k) {
              tmap_bwt_occ4(bwt, prev_k-1, cntk);
              for(i=0;i<4;i++) {
                  next[i].k = cntk[i] + bwt->L2[i] + 1;
                  // put it in the hash
                  tmap_bwt_match_hash_put_k(hash, prev_k, i, next[i].k);
              }
          }
          for(i=0;i<4;i++) {
              next[i].offset = offset + 1;
              next[i].hi = UINT32_MAX;
              // ignore k;
              next[i].l = UINT32_MAX;
          }
      }
  }
  else { // use the bwt hash
      uint64_t prev_hi = (NULL == prev) ? 0 : prev->hi;
      for(i=0;i<4;i++) {
          next[i].offset = offset + 1;
          next[i].hi = (prev_hi << 2) + i;
          next[i].k = bwt->hash_k[next[i].offset-1][next[i].hi];
          next[i].l = UINT32_MAX;
      }
  }
}

inline void
tmap_bwt_match_hash_2occ4(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, 
                            tmap_bwt_match_occ_t next[4], tmap_bwt_match_hash_t *hash)
{
  uint32_t i, offset;
  offset = (NULL == prev) ? 0 : prev->offset;
  if(bwt->hash_width <= offset) { // do not use the bwt hash
      uint32_t cntk[4], cntl[4];
      uint32_t prev_k, prev_l, found_k, found_l;
      prev_k = (NULL == prev) ? 0 : prev->k;
      prev_l = (NULL == prev) ? bwt->seq_len : prev->l;
      if(NULL == hash) {
          tmap_bwt_2occ4(bwt, prev_k-1, prev_l, cntk, cntl);
          for(i=0;i<4;i++) {
              next[i].offset = offset + 1;
              next[i].hi = UINT32_MAX;
              next[i].k = cntk[i] + bwt->L2[i] + 1;
              next[i].l = cntl[i] + bwt->L2[i];
          }
      }
      else { // test the user "hash" for k
          for(i=0;i<4;i++) {
              next[i].k = tmap_bwt_match_hash_get_k(hash, prev_k, i, &found_k); // compute k
              if(0 == found_k) { // for k
                  break;
              }
              next[i].l = tmap_bwt_match_hash_get_l(hash, prev_l, i, &found_l); // compute k
              if(0 == found_l) { // for l
                  break;
              }
          } // TODO: should we use tmap_bwt_occ if we break out of this loop?
          // compute the k value, if necessary
          if(0 == found_k || 0 == found_l) {
              tmap_bwt_2occ4(bwt, prev_k-1, prev_l, cntk, cntl);
              for(i=0;i<4;i++) {
                  next[i].k = cntk[i] + bwt->L2[i] + 1;
                  next[i].l = cntl[i] + bwt->L2[i];
                  // put it in the hash
                  tmap_bwt_match_hash_put_k(hash, prev_k, i, next[i].k);
                  tmap_bwt_match_hash_put_l(hash, prev_l, i, next[i].l);
              }
          }
          for(i=0;i<4;i++) {
              next[i].offset = offset + 1;
              next[i].hi = UINT32_MAX;
              // ignore k and l
          }
      }
  }
  else { // use the bwt hash
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
tmap_bwt_match_hash_cal_width_forward(const tmap_bwt_t *bwt, int len, const char *str, 
                                      tmap_bwt_match_width_t *width, tmap_bwt_match_hash_t *hash)
{
  // 'width[i]' is the lower bound of the number of differences in str[i,len]
  tmap_bwt_match_occ_t prev, next;
  int32_t i, bid;
  
  prev.k = 0; prev.l = bwt->seq_len;
  prev.offset = 0;
  prev.hi = 0;

  bid = 0;
  for(i=len-1;0<=i;i--) {
      uint8_t c = (int)str[i];
      if(c < 4) {
          tmap_bwt_match_hash_2occ(bwt, &prev, c, &next, hash);
      }
      if(next.l < next.k || UINT32_MAX == next.k || 3 < c) { // new width
          next.k = 0;
          next.l = bwt->seq_len;
          next.offset = 0;
          next.hi = 0;
          bid++;
      }
      width[i].w = next.l - next.k + 1;
      width[i].bid = bid;
      prev = next;
  }
}

void
tmap_bwt_match_hash_cal_width_reverse(const tmap_bwt_t *bwt, int len, const char *str, 
                                 tmap_bwt_match_width_t *width, tmap_bwt_match_hash_t *hash)
{
  // 'width[i]' is the lower bound of the number of differences in str[0,i]
  tmap_bwt_match_occ_t prev, next;
  int32_t i, bid;
  
  prev.k = 0; prev.l = bwt->seq_len;
  prev.offset = 0;
  prev.hi = 0;

  bid = 0;
  for(i=0;i<len;i++) {
      uint8_t c = (int)str[i];
      if(c < 4) {
          tmap_bwt_match_hash_2occ(bwt, &prev, c, &next, hash);
      }
      if(next.l < next.k || UINT32_MAX == next.k || 3 < c) { // new width
          next.k = 0;
          next.l = bwt->seq_len;
          next.offset = 0;
          next.hi = 0;
          bid++;
      }
      width[i].w = next.l - next.k + 1;
      width[i].bid = bid;
      prev = next;
  }
  width[len].w = 0;
  width[len].bid = ++bid;
}

static inline uint32_t
tmap_bwt_match_hash_forward_init(const tmap_bwt_t *bwt, int len, const uint8_t *str,
                                    tmap_bwt_match_occ_t *match_sa)
{
  int32_t i;
  uint8_t c = 0;
  len = (len < bwt->hash_width) ? len : bwt->hash_width;
  match_sa->hi = 0;
  for(i=0;i<len;i++) {
      c = str[i];
      if(TMAP_UNLIKELY(3 < c)) {
          match_sa->offset = i+1;
          match_sa->k = bwt->hash_k[match_sa->offset-1][match_sa->hi];
          match_sa->l = bwt->hash_l[match_sa->offset-1][match_sa->hi];
          return 0; 
      }
      match_sa->hi = (match_sa->hi << 2) + c;
  }
  match_sa->offset = len;
  match_sa->k = bwt->hash_k[match_sa->offset-1][match_sa->hi];
  match_sa->l = bwt->hash_l[match_sa->offset-1][match_sa->hi];
  if(match_sa->l < match_sa->k || UINT32_MAX == match_sa->k) {
      for(i=len-1;0<=i&&1<match_sa->offset;i--) {
          if(bwt->hash_k[match_sa->offset-2][match_sa->hi>>2] <= bwt->hash_l[match_sa->offset-2][match_sa->hi>>2] 
             && UINT32_MAX != bwt->hash_k[match_sa->offset-2][match_sa->hi>>2]) {
              break;
          }
          match_sa->hi >>= 2;
          match_sa->offset--;
          match_sa->k = bwt->hash_k[match_sa->offset-1][match_sa->hi];
          match_sa->l = bwt->hash_l[match_sa->offset-1][match_sa->hi];
      }
      return 0;
  }
  return (match_sa->l - match_sa->k + 1);
}

static inline uint32_t
tmap_bwt_match_hash_reverse_init(const tmap_bwt_t *bwt, int len, const uint8_t *str,
                                    tmap_bwt_match_occ_t *match_sa)
{
  int32_t i, low;
  uint8_t c = 0;
  low = (len < bwt->hash_width) ? 0 : (len - bwt->hash_width);
  match_sa->hi = 0;
  for(i=len-1;low<=i;i--) {
      c = str[i];
      if(TMAP_UNLIKELY(3 < c)) {
          match_sa->offset = len-i;
          match_sa->k = bwt->hash_k[match_sa->offset-1][match_sa->hi];
          match_sa->l = bwt->hash_l[match_sa->offset-1][match_sa->hi];
          return 0; 
      }
      match_sa->hi = (match_sa->hi << 2) + c;
  }
  match_sa->offset = len-low;
  match_sa->k = bwt->hash_k[match_sa->offset-1][match_sa->hi];
  match_sa->l = bwt->hash_l[match_sa->offset-1][match_sa->hi];
  if(match_sa->l < match_sa->k || UINT32_MAX == match_sa->k) {
      for(i=low;i<len&&1<match_sa->offset;i++) {
          if(bwt->hash_k[match_sa->offset-2][match_sa->hi>>2] <= bwt->hash_l[match_sa->offset-2][match_sa->hi>>2]
             && UINT32_MAX != bwt->hash_k[match_sa->offset-2][match_sa->hi>>2]) {
              break;
          }
          match_sa->hi >>= 2;
          match_sa->offset--;
          match_sa->k = bwt->hash_k[match_sa->offset-1][match_sa->hi];
          match_sa->l = bwt->hash_l[match_sa->offset-1][match_sa->hi];
      }
      return 0;
  }
  return (match_sa->l - match_sa->k + 1);
}

uint32_t
tmap_bwt_match_hash_exact(const tmap_bwt_t *bwt, int len, const uint8_t *str, 
                          tmap_bwt_match_occ_t *match_sa, tmap_bwt_match_hash_t *hash)
{
  int32_t i;
  uint8_t c = 0;
  tmap_bwt_match_occ_t prev, next;

  if(0 < bwt->hash_width && 0 < len) { 
      if(0 == tmap_bwt_match_hash_forward_init(bwt, len, str, &prev)) {
          if(NULL != match_sa) {
              (*match_sa) = prev;
          }
          return 0;
      }
  }
  else {
      prev.k = 0; prev.l = bwt->seq_len;
      prev.offset = 0;
      prev.hi = 0;
  }

  for(i=prev.offset;i<len;i++) {
      c = str[i];
      if(TMAP_UNLIKELY(3 < c)) {
          prev.offset++;
          prev.k = prev.l + 1;
          break;
      }
      tmap_bwt_match_hash_2occ(bwt, &prev, c, &next, hash);
      prev = next;
      if(next.k > next.l || UINT32_MAX == next.k) break; // no match
  }
  if(NULL != match_sa) {
      (*match_sa) = prev;
  }
  if(3 < c || prev.k > prev.l || UINT32_MAX == prev.k) return 0; // no match
  return prev.l - prev.k + 1;
}

uint32_t
tmap_bwt_match_hash_exact_reverse(const tmap_bwt_t *bwt, int len, const uint8_t *str, 
                             tmap_bwt_match_occ_t *match_sa, tmap_bwt_match_hash_t *hash)
{
  int32_t i;
  uint8_t c = 0;
  tmap_bwt_match_occ_t prev, next;
  
  if(0 < bwt->hash_width && 0 < len) { 
      if(0 == tmap_bwt_match_hash_reverse_init(bwt, len, str, &prev)) {
          if(NULL != match_sa) {
              (*match_sa) = prev;
          }
          return 0;
      }
  }
  else {
      prev.k = 0; prev.l = bwt->seq_len;
      prev.offset = 0;
      prev.hi = 0;
  }

  for(i=len-prev.offset-1;0<=i;i--) {
      c = str[i];
      if(TMAP_UNLIKELY(3 < c)) {
          prev.offset++; 
          prev.k = prev.l + 1;
          break;
      }
      tmap_bwt_match_hash_2occ(bwt, &prev, c, &next, hash);
      prev = next;
      if(next.k > next.l || UINT32_MAX == next.k) break; // no match
  }
  if(NULL != match_sa) {
      (*match_sa) = prev;
  }
  if(3 < c || prev.k > prev.l || UINT32_MAX == prev.k) return 0; // no match
  return prev.l - prev.k + 1;
}

uint32_t
tmap_bwt_match_hash_exact_alt(const tmap_bwt_t *bwt, int len, const uint8_t *str, 
                         tmap_bwt_match_occ_t *match_sa, tmap_bwt_match_hash_t *hash)
{
  int i;
  tmap_bwt_match_occ_t next;

  for(i=0;i<len;i++) {
      uint8_t c = str[i];
      if(TMAP_UNLIKELY(c > 3)) return 0; // there is an N here. no match
      tmap_bwt_match_hash_2occ(bwt, match_sa, c, &next, hash);
      (*match_sa) = next;
      if(next.k > next.l || UINT32_MAX == next.k) return 0; // no match
  }
  return match_sa->l - match_sa->k + 1;
}

uint32_t
tmap_bwt_match_hash_exact_alt_reverse(const tmap_bwt_t *bwt, int len, const uint8_t *str, 
                                 tmap_bwt_match_occ_t *match_sa, tmap_bwt_match_hash_t *hash)
{
  int i;
  tmap_bwt_match_occ_t next;

  for(i=len-1;0<=i;i--) {
      uint8_t c = str[i];
      if(TMAP_UNLIKELY(c > 3)) return 0; // there is an N here. no match
      tmap_bwt_match_hash_2occ(bwt, match_sa, c, &next, hash);
      (*match_sa) = next;
      if(next.k > next.l || UINT32_MAX == next.k) return 0; // no match
  }
  return match_sa->l - match_sa->k + 1;
}

uint32_t
tmap_bwt_match_hash_invPsi(const tmap_bwt_t *bwt, uint32_t sa_intv, uint32_t k, uint32_t *s, tmap_bwt_match_hash_t *hash)
{
  uint8_t b0;
  tmap_bwt_match_occ_t prev, next;

  *s = 0;
  prev.k = k;
  prev.l = UINT32_MAX;
  prev.offset = UINT32_MAX; // do not use the bwt hash
  prev.hi = UINT32_MAX;
  while(0 != (prev.k % sa_intv)) {
      (*s)++;
      if(TMAP_LIKELY(prev.k != bwt->primary)) { // likely
          if(prev.k < bwt->primary) {
              b0 = tmap_bwt_B0(bwt, prev.k);
          }
          else {
              b0 = tmap_bwt_B0(bwt, prev.k-1);
          }

          prev.k++; // since k-1 will be used in tmap_bwt_match_hash_occ
          tmap_bwt_match_hash_occ(bwt, &prev, b0, &next, hash); 
          prev.k = next.k - 1; // since there was an extra one added in tmap_bwt_match_hash_occ
      }
      else {
          prev.k = 0;
      }
  }

  return prev.k;
}
