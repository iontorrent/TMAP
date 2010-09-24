#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "../util/fmap_alloc.h"
#include "fmap_sa.h"
#include "fmap_bwtl.h"

fmap_bwtl_t *
fmap_bwtl_seq2bwtl(int32_t len, const uint8_t *seq)
{
  fmap_bwtl_t *bwtl = NULL;
  int32_t i;

  bwtl = fmap_calloc(1, sizeof(fmap_bwtl_t), "bwtl");
  bwtl->seq_len = len;

  // calculate bwtl->bwt
  uint8_t *s;
  bwtl->sa = fmap_calloc(len + 1, sizeof(uint32_t), "bwtl->sa");
  fmap_sa_gen_short(seq, (int32_t*)bwtl->sa, len);
  s = fmap_calloc(len + 1, sizeof(uint8_t), "s");
  for(i = 0; i <= len; ++i) {
      if(bwtl->sa[i] == 0) bwtl->primary = i;
      else s[i] = seq[bwtl->sa[i] - 1];
  }
  for(i = bwtl->primary; i < len; ++i) s[i] = s[i + 1];
  bwtl->bwt_size = (len + 15) / 16;
  bwtl->bwt = fmap_calloc(bwtl->bwt_size, sizeof(uint32_t), "bwtl->bwt");
  for(i = 0; i < len; ++i)
    bwtl->bwt[i>>4] |= s[i] << ((15 - (i&15)) << 1);
  free(s);

  // calculate bwtl->occ
  uint32_t c[4];
  bwtl->n_occ = (len + 15) / 16 * 4;
  bwtl->occ = fmap_calloc(bwtl->n_occ, sizeof(uint32_t), "bwtl->occ");
  memset(c, 0, 16);
  for(i = 0; i < len; ++i) {
      if(i % 16 == 0)
        memcpy(bwtl->occ + (i/16) * 4, c, 16);
      ++c[fmap_bwtl_B0(bwtl, i)];
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
fmap_bwtl_occ(const fmap_bwtl_t *bwtl, uint32_t k, uint8_t c)
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
fmap_bwtl_occ4(const fmap_bwtl_t *bwtl, uint32_t k, uint32_t cnt[4])
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
fmap_bwtl_2occ4(const fmap_bwtl_t *bwtl, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4])
{
  fmap_bwtl_occ4(bwtl, k, cntk);
  fmap_bwtl_occ4(bwtl, l, cntl);
}

void 
fmap_bwtl_destroy(fmap_bwtl_t *bwtl)
{
  if(NULL != bwtl) {
      free(bwtl->occ); free(bwtl->bwt); free(bwtl->sa);
      free(bwtl);
  }
}
