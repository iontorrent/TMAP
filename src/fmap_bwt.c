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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "fmap_io.h"
#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_refseq.h"
#include "fmap_bwt.h"

fmap_bwt_t *
fmap_bwt_read(const char *fn_fasta)
{
  fmap_bwt_t *bwt;
  char *fn_bwt = NULL;
  fmap_file_t *fp_bwt = NULL;

  fn_bwt = fmap_refseq_get_file_name(fn_fasta, FMAP_REFSEQ_BWT_FILE);
  fp_bwt = fmap_file_fopen(fn_bwt, "rb", FMAP_REFSEQ_BWT_COMPRESSION);

  bwt = fmap_calloc(1, sizeof(fmap_bwt_t), "bwt");

  if(1 != fmap_file_fread(&bwt->bwt_size, sizeof(uint32_t), 1, fp_bwt)) {
      fmap_error(NULL, Exit, ReadFileError);
  }

  bwt->bwt = fmap_calloc(bwt->bwt_size, sizeof(uint32_t), "bwt->bwt");
  if(1 != fmap_file_fread(&bwt->primary, sizeof(uint32_t), 1, fp_bwt)
     || 4 != fmap_file_fread(bwt->L2+1, sizeof(uint32_t), 4, fp_bwt)
     || 1 != fmap_file_fread(&bwt->occ_interval, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fread(&bwt->seq_len, sizeof(uint32_t), 1, fp_bwt)
     || bwt->bwt_size != fmap_file_fread(bwt->bwt, sizeof(uint32_t), bwt->bwt_size, fp_bwt)) {
      fmap_error(NULL, Exit, ReadFileError);
  }
  fmap_bwt_gen_cnt_table(bwt);

  fmap_file_fclose(fp_bwt);
  free(fn_bwt);

  return bwt;
}

void 
fmap_bwt_write(const char *fn_fasta, fmap_bwt_t *bwt)
{
  char *fn_bwt = NULL;
  fmap_file_t *fp_bwt = NULL;

  fn_bwt = fmap_refseq_get_file_name(fn_fasta, FMAP_REFSEQ_BWT_FILE);
  fp_bwt = fmap_file_fopen(fn_bwt, "wb", FMAP_REFSEQ_BWT_COMPRESSION);

  if(1 != fmap_file_fwrite(&bwt->bwt_size, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fwrite(&bwt->primary, sizeof(uint32_t), 1, fp_bwt) 
     || 4 != fmap_file_fwrite(bwt->L2+1, sizeof(uint32_t), 4, fp_bwt)
     || 1 != fmap_file_fwrite(&bwt->occ_interval, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fwrite(&bwt->seq_len, sizeof(uint32_t), 1, fp_bwt)
     || bwt->bwt_size != fmap_file_fwrite(bwt->bwt, sizeof(uint32_t), bwt->bwt_size, fp_bwt)) {
      fmap_error(NULL, Exit, WriteFileError);
  }

  fmap_file_fclose(fp_bwt);
  free(fn_bwt);
}

void 
fmap_bwt_destroy(fmap_bwt_t *bwt)
{
  if(bwt == 0) return;
  free(bwt->bwt);
  free(bwt);
}

void
fmap_bwt_update(fmap_bwt_t *bwt, uint32_t occ_interval)
{
  uint32_t i, k, c[4], n_occ;
  uint32_t *buf = NULL;

#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

  bwt->occ_interval = occ_interval;
  n_occ = (bwt->seq_len + occ_interval - 1) / occ_interval + 1;
  bwt->bwt_size += n_occ * 4; // the new size
  buf = fmap_calloc(bwt->bwt_size, sizeof(uint32_t), "buf"); // will be the new bwt
  c[0] = c[1] = c[2] = c[3] = 0;
  for (i = k = 0; i < bwt->seq_len; ++i) {
      if (i % occ_interval == 0) {
          memcpy(buf + k, c, sizeof(uint32_t) * 4);
          k += 4;
      }
      if (i % 16 == 0) buf[k++] = bwt->bwt[i>>4];
      ++c[bwt_B00(bwt, i)];
  }
  // the last element
  memcpy(buf + k, c, sizeof(uint32_t) * 4);
  if(k + 4 != bwt->bwt_size) {
      fmap_error(NULL, Exit, OutOfRange);
  }
  // update bwt
  free(bwt->bwt); bwt->bwt = buf;
}

void 
fmap_bwt_gen_cnt_table(fmap_bwt_t *bwt)
{
  int i, j;
  for(i = 0; i != 256; ++i) {
      uint32_t x = 0;
      for(j = 0; j != 4; ++j)
        x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
      bwt->cnt_table[i] = x;
  }
}

static inline int 
__occ_aux(uint64_t y, int c)
{
  // reduce nucleotide counting to bits counting
  y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
  // count the number of 1s in y
  y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
  return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

inline uint32_t 
bwt_occ(const fmap_bwt_t *bwt, uint32_t k, uint8_t c)
{
  uint32_t n, l, j;
  uint32_t *p;

  if(k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
  if(k == (uint32_t)(-1)) return 0;
  if(k >= bwt->primary) --k; // because $ is not in bwt

  // retrieve Occ at k/bwt->occ_interval
  n = (p = fmap_bwt_occ_intv(bwt, k))[c];

  p += 4; // jump to the start of the first BWT cell

  // calculate Occ up to the last k/32
  j = k >> 5 << 5;
  for(l = k/bwt->occ_interval*bwt->occ_interval; l < j; l += 32, p += 2)
    n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);

  // calculate Occ
  n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
  if(c == 0) n -= ~k&31; // corrected for the masked bits

  return n;
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
inline void 
bwt_2occ(const fmap_bwt_t *bwt, uint32_t k, uint32_t l, uint8_t c, uint32_t *ok, uint32_t *ol)
{
  uint32_t _k, _l;
  if(k == l) {
      *ok = *ol = bwt_occ(bwt, k, c);
      return;
  }
  _k = (k >= bwt->primary)? k-1 : k;
  _l = (l >= bwt->primary)? l-1 : l;
  if(_l/bwt->occ_interval != _k/bwt->occ_interval || k == (uint32_t)(-1) || l == (uint32_t)(-1)) {
      *ok = bwt_occ(bwt, k, c);
      *ol = bwt_occ(bwt, l, c);
  } else {
      uint32_t m, n, i, j;
      uint32_t *p;
      if(k >= bwt->primary) --k;
      if(l >= bwt->primary) --l;
      n = (p = fmap_bwt_occ_intv(bwt, k))[c];
      p += 4;
      // calculate *ok
      j = k >> 5 << 5;
      for(i = k/bwt->occ_interval*bwt->occ_interval; i < j; i += 32, p += 2)
        n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
      m = n;
      n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
      if(c == 0) n -= ~k&31; // corrected for the masked bits
      *ok = n;
      // calculate *ol
      j = l >> 5 << 5;
      for(; i < j; i += 32, p += 2)
        m += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
      m += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~l&31)<<1)) - 1), c);
      if(c == 0) m -= ~l&31; // corrected for the masked bits
      *ol = m;
  }
}

#define __occ_aux4(bwt, b)											\
  ((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
   + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

inline void 
bwt_occ4(const fmap_bwt_t *bwt, uint32_t k, uint32_t cnt[4])
{
  uint32_t l, j, x;
  uint32_t *p;
  if(k == (uint32_t)(-1)) {
      memset(cnt, 0, 4 * sizeof(uint32_t));
      return;
  }
  if(k >= bwt->primary) --k; // because $ is not in bwt
  p = fmap_bwt_occ_intv(bwt, k);
  memcpy(cnt, p, 16);
  p += 4;
  j = k >> 4 << 4;
  for(l = k / bwt->occ_interval * bwt->occ_interval, x = 0; l < j; l += 16, ++p)
    x += __occ_aux4(bwt, *p);
  x += __occ_aux4(bwt, *p & ~((1U<<((~k&15)<<1)) - 1)) - (~k&15);
  cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

// an analogy to bwt_occ4() but more efficient, requiring k <= l
inline void 
bwt_2occ4(const fmap_bwt_t *bwt, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4])
{
  uint32_t _k, _l;
  if(k == l) {
      bwt_occ4(bwt, k, cntk);
      memcpy(cntl, cntk, 4 * sizeof(uint32_t));
      return;
  }
  _k = (k >= bwt->primary)? k-1 : k;
  _l = (l >= bwt->primary)? l-1 : l;
  if(_l/bwt->occ_interval != _k/bwt->occ_interval || k == (uint32_t)(-1) || l == (uint32_t)(-1)) {
      bwt_occ4(bwt, k, cntk);
      bwt_occ4(bwt, l, cntl);
  } else {
      uint32_t i, j, x, y;
      uint32_t *p;
      int cl[4];
      if(k >= bwt->primary) --k; // because $ is not in bwt
      if(l >= bwt->primary) --l;
      cl[0] = cl[1] = cl[2] = cl[3] = 0;
      p = fmap_bwt_occ_intv(bwt, k);
      memcpy(cntk, p, 4 * sizeof(uint32_t));
      p += 4;
      // prepare cntk[]
      j = k >> 4 << 4;
      for(i = k / bwt->occ_interval * bwt->occ_interval, x = 0; i < j; i += 16, ++p)
        x += __occ_aux4(bwt, *p);
      y = x;
      x += __occ_aux4(bwt, *p & ~((1U<<((~k&15)<<1)) - 1)) - (~k&15);
      // calculate cntl[] and finalize cntk[]
      j = l >> 4 << 4;
      for(; i < j; i += 16, ++p) y += __occ_aux4(bwt, *p);
      y += __occ_aux4(bwt, *p & ~((1U<<((~l&15)<<1)) - 1)) - (~l&15);
      memcpy(cntl, cntk, 16);
      cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
      cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
  }
}

int 
bwt_match_exact(const fmap_bwt_t *bwt, int len, const uint8_t *str, uint32_t *sa_begin, uint32_t *sa_end)
{
  uint32_t k, l, ok, ol;
  int i;
  k = 0; l = bwt->seq_len;
  for(i = len - 1; i >= 0; --i) {
      uint8_t c = str[i];
      if(c > 3) return 0; // no match
      bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
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
bwt_match_exact_alt(const fmap_bwt_t *bwt, int len, const uint8_t *str, uint32_t *k0, uint32_t *l0)
{
  int i;
  uint32_t k, l, ok, ol;
  k = *k0; l = *l0;
  for(i = len - 1; i >= 0; --i) {
      uint8_t c = str[i];
      if(c > 3) return 0; // there is an N here. no match
      bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
      k = bwt->L2[c] + ok + 1;
      l = bwt->L2[c] + ol;
      if(k > l) return 0; // no match
  }
  *k0 = k; *l0 = l;
  return l - k + 1;
}
