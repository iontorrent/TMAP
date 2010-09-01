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
#include "fmap_bwt_gen.h"
#include "fmap_progress.h"
#include "fmap_definitions.h"
#include "fmap_bwt.h"

fmap_bwt_t *
fmap_bwt_read(const char *fn_fasta, uint32_t is_rev)
{
  fmap_bwt_t *bwt;
  char *fn_bwt = NULL;
  fmap_file_t *fp_bwt = NULL;

  fn_bwt = fmap_get_file_name(fn_fasta, (0 == is_rev) ? FMAP_BWT_FILE : FMAP_REV_BWT_FILE);
  fp_bwt = fmap_file_fopen(fn_bwt, "rb", (0 == is_rev) ? FMAP_BWT_COMPRESSION : FMAP_REV_BWT_COMPRESSION);

  bwt = fmap_calloc(1, sizeof(fmap_bwt_t), "bwt");

  if(1 != fmap_file_fread(&bwt->version_id, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fread(&bwt->bwt_size, sizeof(uint32_t), 1, fp_bwt)) {
      fmap_error(NULL, Exit, ReadFileError);
  }

  if(bwt->version_id != FMAP_VERSION_ID) {
      fmap_error("version id did not match", Exit, ReadFileError);
  }

  bwt->bwt = fmap_calloc(bwt->bwt_size, sizeof(uint32_t), "bwt->bwt");

  if(1 != fmap_file_fread(&bwt->hash_width, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fread(&bwt->primary, sizeof(uint32_t), 1, fp_bwt)
     || 4 != fmap_file_fread(bwt->L2+1, sizeof(uint32_t), 4, fp_bwt)
     || 1 != fmap_file_fread(&bwt->occ_interval, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fread(&bwt->seq_len, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fread(&bwt->is_rev, sizeof(uint32_t), 1, fp_bwt)
     || bwt->bwt_size != fmap_file_fread(bwt->bwt, sizeof(uint32_t), bwt->bwt_size, fp_bwt)) {
      fmap_error(NULL, Exit, ReadFileError);
  }
  if(0 < bwt->hash_width) {
      bwt->hash_length = 1 << (bwt->hash_width * 2); // 4^{hash_width} entries
      bwt->hash_k = fmap_malloc(bwt->hash_length*sizeof(uint32_t), "bwt->hash_k");
      bwt->hash_l = fmap_malloc(bwt->hash_length*sizeof(uint32_t), "bwt->hash_l");
      if(bwt->hash_length != fmap_file_fread(bwt->hash_k, sizeof(uint32_t), bwt->hash_length, fp_bwt)
         || bwt->hash_length != fmap_file_fread(bwt->hash_l, sizeof(uint32_t), bwt->hash_length, fp_bwt)) {
          fmap_error(NULL, Exit, ReadFileError);
      }
  }
  else {
      bwt->hash_length = 0;
      bwt->hash_k = bwt->hash_l = NULL;
  }

  fmap_bwt_gen_cnt_table(bwt);

  if(is_rev != bwt->is_rev) {
      fmap_error("is_rev != bwt->is_rev", Exit, OutOfRange);
  }

  fmap_file_fclose(fp_bwt);
  free(fn_bwt);

  return bwt;
}

void 
fmap_bwt_write(const char *fn_fasta, fmap_bwt_t *bwt, uint32_t is_rev)
{
  char *fn_bwt = NULL;
  fmap_file_t *fp_bwt = NULL;

  if(is_rev != bwt->is_rev) {
      fmap_error("is_rev != bwt->is_rev", Exit, OutOfRange);
  }

  fn_bwt = fmap_get_file_name(fn_fasta, (0 == is_rev) ? FMAP_BWT_FILE : FMAP_REV_BWT_FILE);
  fp_bwt = fmap_file_fopen(fn_bwt, "wb", (0 == is_rev) ? FMAP_BWT_COMPRESSION : FMAP_REV_BWT_COMPRESSION);

  if(1 != fmap_file_fwrite(&bwt->version_id, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fwrite(&bwt->bwt_size, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fwrite(&bwt->hash_width, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fwrite(&bwt->primary, sizeof(uint32_t), 1, fp_bwt) 
     || 4 != fmap_file_fwrite(bwt->L2+1, sizeof(uint32_t), 4, fp_bwt)
     || 1 != fmap_file_fwrite(&bwt->occ_interval, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fwrite(&bwt->seq_len, sizeof(uint32_t), 1, fp_bwt)
     || 1 != fmap_file_fwrite(&bwt->is_rev, sizeof(uint32_t), 1, fp_bwt)
     || bwt->bwt_size != fmap_file_fwrite(bwt->bwt, sizeof(uint32_t), bwt->bwt_size, fp_bwt)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
  if(0 < bwt->hash_width) {
      if(bwt->hash_length != fmap_file_fwrite(bwt->hash_k, sizeof(uint32_t), bwt->hash_length, fp_bwt)
         || bwt->hash_length != fmap_file_fwrite(bwt->hash_l, sizeof(uint32_t), bwt->hash_length, fp_bwt)) {
          fmap_error(NULL, Exit, WriteFileError);
      }
  }

  fmap_file_fclose(fp_bwt);
  free(fn_bwt);
}

void 
fmap_bwt_destroy(fmap_bwt_t *bwt)
{
  if(bwt == 0) return;
  free(bwt->hash_k);
  free(bwt->hash_l);
  free(bwt->bwt);
  free(bwt);
}

void
fmap_bwt_update_occ_interval(fmap_bwt_t *bwt, uint32_t occ_interval)
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

static inline void
fmap_bwt_gen_helper(fmap_bwt_t *bwt, uint32_t hash_value, uint32_t level, 
                    uint32_t **k, uint32_t **l)
{
  uint32_t i, j;
  // TODO: could make this tail-recursive for speed

  if(bwt->hash_width < level) {
      // do nothing
      return;
  }
  else if(bwt->hash_width == level) { // at the leaves
      // store values
      for(i=0;i<4;i++) {
          bwt->hash_k[(hash_value << 2) + i] = k[level-1][i] + bwt->L2[i] + 1;
          bwt->hash_l[(hash_value << 2) + i] = l[level-1][i] + bwt->L2[i]; 
      }
      return;
  }
  else { // at a non-leaf
      for(i=0;i<4;i++) {
          if(l[level-1][i] < k[level-1][i]) {
              uint32_t lower, upper;

              lower = (hash_value << 2) + i;
              upper = ((hash_value + 1) << 2);
              if(upper <= lower) {
                  fmap_error("Control reached unexpected point", Exit, OutOfRange);
              }

              // all at this level or below are NULL
              for(j=lower;j<upper;j++) {
                  bwt->hash_k[j] = UINT32_MAX;
                  bwt->hash_l[j] = UINT32_MAX;
              }

              level--;
              return;
          }
          else { // descend deeper
              bwt_2occ4(bwt, 
                        k[level-1][i]+bwt->L2[i], l[level-1][i]+bwt->L2[i], 
                        k[level], l[level]); // get occ values
              fmap_bwt_gen_helper(bwt, (hash_value << 2) + i, level+1, k, l); 
          }
      }
      return;
  }
}

void
fmap_bwt_gen_hash(fmap_bwt_t *bwt, uint32_t hash_width)
{
  uint32_t **k=NULL, **l=NULL;
  uint32_t i;

  fmap_progress_set_start_time(clock());
  fmap_progress_print("constructing the occurence hash for the BWT string");

  bwt->hash_length = 1 << (hash_width * 2); // 4^{hash_width} entries

  // TODO check if hash_width falls within acceptable limits
  if(UINT32_MAX <= bwt->hash_length-1) {
      fmap_error("UINT32_MAX <= bwt->hash_length-1", Exit, OutOfRange);
  }

  bwt->hash_width = hash_width;
  bwt->hash_k = fmap_malloc(sizeof(uint32_t)*bwt->hash_length, "bwt->hash_k");
  bwt->hash_l = fmap_malloc(sizeof(uint32_t)*bwt->hash_length, "bwt->hash_l");

  // working space
  k = fmap_malloc(sizeof(uint32_t*)*hash_width, "k");
  l = fmap_malloc(sizeof(uint32_t*)*hash_width, "l");
  for(i=0;i<hash_width;i++) {
      k[i] = fmap_malloc(sizeof(uint32_t)*4, "k[i]");
      l[i] = fmap_malloc(sizeof(uint32_t)*4, "l[i]");
  }

  // Get first level
  bwt_2occ4(bwt, -1, bwt->seq_len, k[0], l[0]);
  fmap_bwt_gen_helper(bwt, 0, 1, k, l); 

  uint32_t sum = 0;
  for(i=0;i<bwt->hash_length;i++) {
      if(UINT32_MAX != bwt->hash_k[i]) {
          sum += bwt->hash_l[i] - bwt->hash_k[i] + 1;
          //fprintf(stderr, "[%d] k[%d]=%u l[%u]=%u\n", i, i, bwt->hash_k[i], i, bwt->hash_l[i]);
      }
  }
  //fprintf(stderr, "bwt->seq_len=%u sum=%u\n", bwt->seq_len, sum);
  if(sum != bwt->seq_len - bwt->hash_width + 1) {
      fmap_error("sum != bwt->seq_len - bwt->hash_width + 1", Exit, OutOfRange);
  }

  for(i=0;i<hash_width;i++) {
      free(k[i]);
      free(l[i]);
  }
  free(k);
  free(l);

  fmap_progress_print2("constructed the occurence hash for the BWT string");
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
  for(i=0;i<len;i++) {
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
  for(i=0;i<len;i++) {
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

int
fmap_bwt_pac2bwt_main(int argc, char *argv[])
{
  int c, is_large = 0, occ_interval = FMAP_BWT_OCC_INTERVAL, hash_width = FMAP_BWT_HASH_WIDTH;
  while((c = getopt(argc, argv, "o:l")) >= 0) {
      switch(c) {
        case 'l': is_large= 1; break;
        case 'o': occ_interval = atoi(optarg); break;
        case 'w': hash_width = atoi(optarg); break;
        default: return 1;
      }
  }
  if (argc < 2) {
      fprintf(stderr, "Usage: %s %s [-l -o INT -w INT] <in.fasta>\n", PACKAGE, argv[0]);
      return 1;
  }
  if(0 < occ_interval && 0 != (occ_interval % 16)) {
      fmap_error("option -o out of range", Exit, CommandLineArgument);
  }

  fmap_bwt_pac2bwt(argv[2], is_large, occ_interval, hash_width);

  return 0;
}
