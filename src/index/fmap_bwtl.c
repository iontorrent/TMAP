#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

#include "../util/tmap_alloc.h"
#include "tmap_sa.h"
#include "tmap_bwtl.h"

tmap_bwtl_t *
tmap_bwtl_seq2bwtl(int32_t len, const uint8_t *seq)
{
  tmap_bwtl_t *bwtl = NULL;
  int32_t i;

  bwtl = tmap_calloc(1, sizeof(tmap_bwtl_t), "bwtl");
  bwtl->seq_len = len;

  // calculate bwtl->bwt
  uint8_t *s;
  bwtl->sa = tmap_calloc(len + 1, sizeof(uint32_t), "bwtl->sa");
  tmap_sa_gen_short(seq, (int32_t*)bwtl->sa, len);
  s = tmap_calloc(len + 1, sizeof(uint8_t), "s");
  for(i = 0; i <= len; ++i) {
      if(bwtl->sa[i] == 0) bwtl->primary = i;
      else s[i] = seq[bwtl->sa[i] - 1];
  }
  for(i = bwtl->primary; i < len; ++i) s[i] = s[i + 1];
  bwtl->bwt_size = (len + 15) / 16;
  bwtl->bwt = tmap_calloc(bwtl->bwt_size, sizeof(uint32_t), "bwtl->bwt");
  for(i = 0; i < len; ++i)
    bwtl->bwt[i>>4] |= s[i] << ((15 - (i&15)) << 1);
  free(s);

  // calculate bwtl->occ
  uint32_t c[4];
  bwtl->n_occ = (len + 15) / 16 * 4;
  bwtl->occ = tmap_calloc(bwtl->n_occ, sizeof(uint32_t), "bwtl->occ");
  memset(c, 0, 16);
  for(i = 0; i < len; ++i) {
      if(i % 16 == 0)
        memcpy(bwtl->occ + (i/16) * 4, c, 16);
      ++c[tmap_bwtl_B0(bwtl, i)];
  }
  memcpy(bwtl->L2+1, c, 16);
  for(i = 2; i < 5; ++i) bwtl->L2[i] += bwtl->L2[i-1];

  // generate cnt_table
  for(i = 0; i != 256; ++i) {
      u_int32_t j, x = 0;
      for(j = 0; j != 4; ++j)
        x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
      bwtl->cnt_table[i] = x;
  }
  return bwtl;
}

inline uint32_t 
tmap_bwtl_occ(const tmap_bwtl_t *bwtl, uint32_t k, uint8_t c)
{
  uint32_t n, b;
  if(k == bwtl->seq_len) return bwtl->L2[c+1] - bwtl->L2[c];
  if(k == (uint32_t)(-1)) return 0;
  if(k >= bwtl->primary) --k; // because $ is not in bwt
  n = bwtl->occ[k/16<<2|c];
  b = bwtl->bwt[k/16] & ~((1U<<((15-(k&15))<<1)) - 1);
  n += (bwtl->cnt_table[b&0xff] + bwtl->cnt_table[b>>8&0xff]
        + bwtl->cnt_table[b>>16&0xff] + bwtl->cnt_table[b>>24]) >> (c<<3) & 0xff;
  if(c == 0) n -= 15 - (k&15); // corrected for the masked bits
  return n;
}

inline void 
tmap_bwtl_occ4(const tmap_bwtl_t *bwtl, uint32_t k, uint32_t cnt[4])
{
  uint32_t x, b;
  if(k == (uint32_t)(-1)) {
      memset(cnt, 0, 16);
      return;
  }
  if(k >= bwtl->primary) --k; // because $ is not in bwt
  memcpy(cnt, bwtl->occ + (k>>4<<2), 16);
  b = bwtl->bwt[k>>4] & ~((1U<<((~k&15)<<1)) - 1);
  x = bwtl->cnt_table[b&0xff] + bwtl->cnt_table[b>>8&0xff]
    + bwtl->cnt_table[b>>16&0xff] + bwtl->cnt_table[b>>24];
  x -= 15 - (k&15);
  cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

inline void 
tmap_bwtl_2occ4(const tmap_bwtl_t *bwtl, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4])
{
  tmap_bwtl_occ4(bwtl, k, cntk);
  tmap_bwtl_occ4(bwtl, l, cntl);
}

void 
tmap_bwtl_destroy(tmap_bwtl_t *bwtl)
{
  if(NULL != bwtl) {
      free(bwtl->occ); free(bwtl->bwt); free(bwtl->sa);
      free(bwtl);
  }
}
