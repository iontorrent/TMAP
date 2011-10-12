/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
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
#include <stdint.h>
#include <unistd.h>
#include <assert.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_definitions.h"
#include "../io/tmap_file.h"
#include "tmap_bwt_gen.h"
#include "tmap_bwt.h"
#include "tmap_bwt_match.h" // HERE

static uint64_t
tmap_bwt_get_hash_length(uint64_t i)
{
  return ((uint64_t)1) << (i << 1); // 4^{hash_width} entries
}

tmap_bwt_t *
tmap_bwt_read(const char *fn_fasta, uint32_t is_rev)
{
  tmap_bwt_t *bwt = NULL;
  char *fn_bwt = NULL;
  tmap_file_t *fp_bwt = NULL;

  fn_bwt = tmap_get_file_name(fn_fasta, (0 == is_rev) ? TMAP_BWT_FILE : TMAP_REV_BWT_FILE);
  fp_bwt = tmap_file_fopen(fn_bwt, "rb", (0 == is_rev) ? TMAP_BWT_COMPRESSION : TMAP_REV_BWT_COMPRESSION);

  bwt = tmap_calloc(1, sizeof(tmap_bwt_t), "bwt");

  if(1 != tmap_file_fread(&bwt->version_id, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fread(&bwt->bwt_size, sizeof(uint32_t), 1, fp_bwt)) {
      tmap_error(NULL, Exit, ReadFileError);
  }

  if(bwt->version_id != TMAP_VERSION_ID) {
      tmap_error("version id did not match", Exit, ReadFileError);
  }

  bwt->bwt = tmap_calloc(bwt->bwt_size, sizeof(uint32_t), "bwt->bwt");

  if(1 != tmap_file_fread(&bwt->hash_width, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fread(&bwt->primary, sizeof(uint32_t), 1, fp_bwt)
     || 4 != tmap_file_fread(bwt->L2+1, sizeof(uint32_t), 4, fp_bwt)
     || 1 != tmap_file_fread(&bwt->occ_interval, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fread(&bwt->seq_len, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fread(&bwt->is_rev, sizeof(uint32_t), 1, fp_bwt)
     || bwt->bwt_size != tmap_file_fread(bwt->bwt, sizeof(uint32_t), bwt->bwt_size, fp_bwt)) {
      tmap_error(NULL, Exit, ReadFileError);
  }

  if(0 < bwt->hash_width) {
      int32_t i;
      bwt->hash_k = tmap_malloc(bwt->hash_width*sizeof(uint32_t*), "bwt->hash_k");
      bwt->hash_l = tmap_malloc(bwt->hash_width*sizeof(uint32_t*), "bwt->hash_l");
      for(i=1;i<=bwt->hash_width;i++) {
          uint64_t hash_length = tmap_bwt_get_hash_length(i);
          bwt->hash_k[i-1] = tmap_malloc(hash_length*sizeof(uint32_t), "bwt->hash_k[i-1]");
          bwt->hash_l[i-1] = tmap_malloc(hash_length*sizeof(uint32_t), "bwt->hash_l[i-1]");
          if(hash_length != tmap_file_fread(bwt->hash_k[i-1], sizeof(uint32_t), hash_length, fp_bwt)
             || hash_length != tmap_file_fread(bwt->hash_l[i-1], sizeof(uint32_t), hash_length, fp_bwt)) {
              tmap_error(NULL, Exit, ReadFileError);
          }
      }
  }
  else {
      bwt->hash_k = bwt->hash_l = NULL;
  }

  tmap_bwt_gen_cnt_table(bwt);

  if(is_rev != bwt->is_rev) {
      tmap_error("is_rev != bwt->is_rev", Exit, OutOfRange);
  }

  tmap_file_fclose(fp_bwt);
  free(fn_bwt);

  bwt->is_shm = 0;

  return bwt;
}

void 
tmap_bwt_write(const char *fn_fasta, tmap_bwt_t *bwt, uint32_t is_rev)
{
  char *fn_bwt = NULL;
  tmap_file_t *fp_bwt = NULL;

  if(is_rev != bwt->is_rev) {
      tmap_error("is_rev != bwt->is_rev", Exit, OutOfRange);
  }

  fn_bwt = tmap_get_file_name(fn_fasta, (0 == is_rev) ? TMAP_BWT_FILE : TMAP_REV_BWT_FILE);
  fp_bwt = tmap_file_fopen(fn_bwt, "wb", (0 == is_rev) ? TMAP_BWT_COMPRESSION : TMAP_REV_BWT_COMPRESSION);

  if(1 != tmap_file_fwrite(&bwt->version_id, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fwrite(&bwt->bwt_size, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fwrite(&bwt->hash_width, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fwrite(&bwt->primary, sizeof(uint32_t), 1, fp_bwt) 
     || 4 != tmap_file_fwrite(bwt->L2+1, sizeof(uint32_t), 4, fp_bwt)
     || 1 != tmap_file_fwrite(&bwt->occ_interval, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fwrite(&bwt->seq_len, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fwrite(&bwt->is_rev, sizeof(uint32_t), 1, fp_bwt)
     || bwt->bwt_size != tmap_file_fwrite(bwt->bwt, sizeof(uint32_t), bwt->bwt_size, fp_bwt)) {
      tmap_error(NULL, Exit, WriteFileError);
  }
  if(0 < bwt->hash_width) {
      int32_t i;
      for(i=1;i<=bwt->hash_width;i++) {
          uint64_t hash_length = tmap_bwt_get_hash_length(i);
          if(hash_length != tmap_file_fwrite(bwt->hash_k[i-1], sizeof(uint32_t), hash_length, fp_bwt)
             || hash_length != tmap_file_fwrite(bwt->hash_l[i-1], sizeof(uint32_t), hash_length, fp_bwt)) {
              tmap_error(NULL, Exit, WriteFileError);
          }
      }
  }

  tmap_file_fclose(fp_bwt);
  free(fn_bwt);
}

size_t
tmap_bwt_shm_num_bytes(tmap_bwt_t *bwt)
{
  // returns the number of bytes to allocate for shared memory
  int32_t i;
  size_t n = 0;

  // fixed length data
  n += sizeof(uint32_t); // version id
  n += sizeof(uint32_t); // primary
  n += 5*sizeof(uint32_t); // L2[5]
  n += sizeof(uint32_t); // seq_len
  n += sizeof(uint32_t); // bwt_size;
  n += sizeof(uint32_t); // occ_interval
  n += 256*sizeof(uint32_t); // cnt_table[256]
  n += sizeof(uint32_t); // is_rev
  n += sizeof(uint32_t); // hash_width

  //variable length data
  n += sizeof(uint32_t)*bwt->bwt_size; // bwt
  for(i=1;i<=bwt->hash_width;i++) {
      uint64_t hash_length = tmap_bwt_get_hash_length(i);
      n += sizeof(uint32_t)*hash_length; // hash_k[i-1]
      n += sizeof(uint32_t)*hash_length; // hash_l[i-1]
  }

  return n;
}

size_t
tmap_bwt_shm_read_num_bytes(const char *fn_fasta, uint32_t is_rev)
{
  size_t n = 0;
  tmap_bwt_t *bwt = NULL;
  char *fn_bwt = NULL;
  tmap_file_t *fp_bwt = NULL;

  fn_bwt = tmap_get_file_name(fn_fasta, (0 == is_rev) ? TMAP_BWT_FILE : TMAP_REV_BWT_FILE);
  fp_bwt = tmap_file_fopen(fn_bwt, "rb", (0 == is_rev) ? TMAP_BWT_COMPRESSION : TMAP_REV_BWT_COMPRESSION);

  bwt = tmap_calloc(1, sizeof(tmap_bwt_t), "bwt");

  if(1 != tmap_file_fread(&bwt->version_id, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fread(&bwt->bwt_size, sizeof(uint32_t), 1, fp_bwt)) {
      tmap_error(NULL, Exit, ReadFileError);
  }

  if(bwt->version_id != TMAP_VERSION_ID) {
      tmap_error("version id did not match", Exit, ReadFileError);
  }

  if(1 != tmap_file_fread(&bwt->hash_width, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fread(&bwt->primary, sizeof(uint32_t), 1, fp_bwt)
     || 4 != tmap_file_fread(bwt->L2+1, sizeof(uint32_t), 4, fp_bwt)
     || 1 != tmap_file_fread(&bwt->occ_interval, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fread(&bwt->seq_len, sizeof(uint32_t), 1, fp_bwt)
     || 1 != tmap_file_fread(&bwt->is_rev, sizeof(uint32_t), 1, fp_bwt)) {
      tmap_error(NULL, Exit, ReadFileError);
  }

  // No need to read in bwt->bwt, bwt->hash_k, bwt->hash_l
  bwt->bwt = NULL;
  bwt->hash_k = bwt->hash_l = NULL;

  if(is_rev != bwt->is_rev) {
      tmap_error("is_rev != bwt->is_rev", Exit, OutOfRange);
  }

  tmap_file_fclose(fp_bwt);
  free(fn_bwt);

  bwt->is_shm = 0;

  // get the number of bytes
  n = tmap_bwt_shm_num_bytes(bwt);

  tmap_bwt_destroy(bwt);

  return n;
}

uint8_t *
tmap_bwt_shm_pack(tmap_bwt_t *bwt, uint8_t *buf)
{
  int32_t i;
  // fixed length data
  memcpy(buf, &bwt->version_id, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(buf, &bwt->primary, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(buf, bwt->L2, 5*sizeof(uint32_t)); buf += 5*sizeof(uint32_t);
  memcpy(buf, &bwt->seq_len, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(buf, &bwt->bwt_size, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(buf, &bwt->occ_interval, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(buf, bwt->cnt_table, 256*sizeof(uint32_t)); buf += 256*sizeof(uint32_t);
  memcpy(buf, &bwt->is_rev, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(buf, &bwt->hash_width, sizeof(uint32_t)); buf += sizeof(uint32_t);
  // variable length data
  memcpy(buf, bwt->bwt, bwt->bwt_size*sizeof(uint32_t)); buf += bwt->bwt_size*sizeof(uint32_t);
  for(i=1;i<=bwt->hash_width;i++) {
      uint64_t hash_length = tmap_bwt_get_hash_length(i);
      memcpy(buf, bwt->hash_k[i-1], hash_length*sizeof(uint32_t)); buf += hash_length*sizeof(uint32_t);
      memcpy(buf, bwt->hash_l[i-1], hash_length*sizeof(uint32_t)); buf += hash_length*sizeof(uint32_t);
  }
  return buf;
}

tmap_bwt_t *
tmap_bwt_shm_unpack(uint8_t *buf)
{
  tmap_bwt_t *bwt = NULL;
  int32_t i;

  if(NULL == buf) return NULL;

  bwt = tmap_calloc(1, sizeof(tmap_bwt_t), "bwt");

  // fixed length data
  memcpy(&bwt->version_id, buf, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(&bwt->primary, buf, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(bwt->L2, buf, 5*sizeof(uint32_t)); buf += 5*sizeof(uint32_t);
  memcpy(&bwt->seq_len, buf, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(&bwt->bwt_size, buf, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(&bwt->occ_interval, buf, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(bwt->cnt_table, buf, 256*sizeof(uint32_t)); buf += 256*sizeof(uint32_t);
  memcpy(&bwt->is_rev, buf, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(&bwt->hash_width, buf, sizeof(uint32_t)); buf += sizeof(uint32_t);

  // allocate memory 
  bwt->hash_k = tmap_calloc(bwt->hash_width, sizeof(uint32_t*), "bwt->hash_k");
  bwt->hash_l = tmap_calloc(bwt->hash_width, sizeof(uint32_t*), "bwt->hash_l");

  // variable length data
  bwt->bwt = (uint32_t*)buf; buf += bwt->bwt_size*sizeof(uint32_t);
  for(i=1;i<=bwt->hash_width;i++) {
      uint64_t hash_length = tmap_bwt_get_hash_length(i);
      bwt->hash_k[i-1] = (uint32_t*)buf; buf += hash_length*sizeof(uint32_t);
      bwt->hash_l[i-1] = (uint32_t*)buf; buf += hash_length*sizeof(uint32_t);
  }

  bwt->is_shm = 1;

  return bwt;
}

void 
tmap_bwt_destroy(tmap_bwt_t *bwt)
{
  int32_t i;
  if(bwt == 0) return;
  if(1 == bwt->is_shm) {
      free(bwt->hash_k);
      free(bwt->hash_l);
      free(bwt);
  }
  else {
      if(NULL != bwt->hash_k) {
          for(i=0;i<bwt->hash_width;i++) {
              free(bwt->hash_k[i]);
          }
      }
      if(NULL != bwt->hash_l) {
          for(i=0;i<bwt->hash_width;i++) {
              free(bwt->hash_l[i]);
          }
      }
      free(bwt->hash_k);
      free(bwt->hash_l);
      free(bwt->bwt);
      free(bwt);
  }
}

void
tmap_bwt_update_occ_interval(tmap_bwt_t *bwt, uint32_t occ_interval)
{
  uint32_t i, k, c[4], n_occ;
  uint32_t *buf = NULL;

  // 2-bit DNA packed
#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

  if(occ_interval < 16 || 0 != occ_interval % 16) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }

  n_occ = ((bwt->seq_len + occ_interval - 1) / occ_interval) + 1; // the number of occurences to store, on top of the bwt string
  bwt->occ_interval = occ_interval;
  bwt->bwt_size += n_occ * 4; // the new size
  buf = tmap_calloc(bwt->bwt_size, sizeof(uint32_t), "buf"); // will be the new bwt
  c[0] = c[1] = c[2] = c[3] = 0;
  for (i = k = 0; i < bwt->seq_len; ++i) {
      // store the occurrences
      if (i % occ_interval == 0) {
          memcpy(buf + k, c, sizeof(uint32_t) * 4);
          k += 4;
      }
      // store the bwt string
      if (i % 16 == 0) buf[k++] = bwt->bwt[i>>4];
      ++c[bwt_B00(bwt, i)];
  }
  // the last element
  memcpy(buf + k, c, sizeof(uint32_t) * 4);
  if(k + 4 != bwt->bwt_size) {
      tmap_error(NULL, Exit, OutOfRange);
  }
  // update bwt
  free(bwt->bwt); bwt->bwt = buf;
}

void 
tmap_bwt_gen_cnt_table(tmap_bwt_t *bwt)
{
  int i, j;
  for(i = 0; i != 256; ++i) {
      uint32_t x = 0;
      for(j = 0; j != 4; ++j)
        x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
      bwt->cnt_table[i] = x;
  }
}

static void
tmap_bwt_gen_hash_helper(tmap_bwt_t *bwt, int32_t len)
{
  int32_t i, n;
  uint8_t *seq = NULL;
  tmap_bwt_match_occ_t match_sa;
  uint32_t sum = 0;
  uint64_t hash_i = 0, low, high;
  //int32_t j;

  seq = tmap_calloc(len, sizeof(uint8_t), "seq");

  i = 0; // length of the current k-mer
  hash_i = 0;
  while(1) {
      if(len < i) {
          tmap_error("ABORT", Exit, OutOfRange);
      }
      if(i == len) {
          // test
          n = tmap_bwt_match_exact(bwt, len, seq, &match_sa);
          bwt->hash_k[len-1][hash_i] = match_sa.k;
          bwt->hash_l[len-1][hash_i] = match_sa.l;
          if(match_sa.k <= match_sa.l) {
              sum += n;
          }

          // find the next base
          i--;
          while(0 <= i && 3 == seq[i]) {
              seq[i] = 0;
              hash_i >>= 2;
              i--;
          }
          if(i < 0) break;
          seq[i]++;
          hash_i++;
          i++;
      }
      else {
          // check the previous hash
          if(0 < i && (UINT32_MAX == bwt->hash_k[i-1][hash_i] || bwt->hash_l[i-1][hash_i] < bwt->hash_k[i-1][hash_i])) {
              // no need to search further
              low = hash_i << (2 * (len - i));
              high = (hash_i+1) << (2 * (len - i));
              while(low < high) {
                  bwt->hash_k[len-1][low] = UINT32_MAX;
                  bwt->hash_l[len-1][low] = UINT32_MAX;
                  low++;
              }
              // find the next base
              i--;
              while(0 <= i && 3 == seq[i]) {
                  seq[i] = 0;
                  hash_i >>= 2;
                  i--;
              }
              if(i < 0) break;
              seq[i]++;
              hash_i++;
              i++;
          }
          else {
              hash_i <<= 2;
              i++;
          }
      }
  }
      
  //fprintf(stderr, "sum=%d (bwt->seq_len - len + 1)=%d\n", sum, bwt->seq_len - len + 1);
  if(sum != bwt->seq_len - len + 1) {
      tmap_error("sum != bwt->seq_len - len + 1", Exit, OutOfRange);
  }
  
  free(seq);
}

void
tmap_bwt_gen_hash(tmap_bwt_t *bwt, uint32_t hash_width)
{
  uint32_t i;
  /*
  uint32_t sum k[4], l[4];
  uint64_t b, j;
  */

  tmap_progress_print("constructing the occurence hash for the BWT string");

  if(bwt->seq_len < hash_width) {
      tmap_error("Hash width was greater than the sequence length, defaulting to the sequence length", Warn, OutOfRange);
      hash_width = bwt->seq_len;
  }

  if(UINT32_MAX <= (1 << (2 * bwt->hash_width)) - 1) {
      tmap_error("Hash width is too great to fit in memory", Exit, OutOfRange);
  }


  // memory for each level
  bwt->hash_k = tmap_malloc(sizeof(uint32_t*)*hash_width, "bwt->hash_k");
  bwt->hash_l = tmap_malloc(sizeof(uint32_t*)*hash_width, "bwt->hash_l");

  bwt->hash_width = 0;
  for(i=1;i<=hash_width;i++) {
      uint64_t hash_length = tmap_bwt_get_hash_length(i);
      
      // allocate memory for this level
      bwt->hash_k[i-1] = tmap_malloc(sizeof(uint32_t)*hash_length, "bwt->hash_k[i-1]");
      bwt->hash_l[i-1] = tmap_malloc(sizeof(uint32_t)*hash_length, "bwt->hash_l[i-1]");

      tmap_bwt_gen_hash_helper(bwt, i);
      bwt->hash_width = i; // updated the hash width
  }
  
  tmap_progress_print2("constructed the occurence hash for the BWT string");
}

static inline int 
__occ_aux32(uint64_t y, int c)
{
  // reduce nucleotide counting to bits counting
  y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
  // count the number of 1s in y
  y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
  return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

static inline int 
__occ_aux16(uint32_t y, int c)
{
  // reduce nucleotide counting to bits counting
  y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x55555555;
  // count the number of 1s in y
  // from http://graphics.stanford.edu/~seander/bithacks.html
  y = y - ((y >> 1) & 0x55555555);                    // reuse input as temporary
  y = (y & 0x33333333) + ((y >> 2) & 0x33333333);     // temp
  return (((y + (y >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24; // count
}

inline uint32_t 
tmap_bwt_occ(const tmap_bwt_t *bwt, uint32_t k, uint8_t c)
{
  uint32_t n, j, l;
  uint32_t *p;

  if(k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
  if(UINT32_MAX == k) return 0;
  if(k >= bwt->primary) --k; // because $ is not in bwt

  // retrieve Occ at k/bwt->occ_interval
  n = (p = tmap_bwt_occ_intv(bwt, k))[c];
  p += 4; // jump to the start of the first BWT cell

  // calculate Occ
  j = k >> 4 << 4; // divide by 16, then multiply by 16, to subtract k % 16.
  for(l = (k/bwt->occ_interval)*bwt->occ_interval; l < j; l += 16, p++) {
      n += __occ_aux16(p[0], c);
  }
  n += __occ_aux16(p[0] & ~((1ul<<((~k&15)<<1)) - 1), c);
  if(c == 0) n -= ~k&15; // corrected for the masked bits

  /*
  // calculate Occ up to the last k/32
  n = (p = tmap_bwt_occ_intv(bwt, k))[c];
  p += 4; // jump to the start of the first BWT cell
  j = k >> 5 << 5;
  for(l = k/bwt->occ_interval*bwt->occ_interval; l < j; l += 32, p += 2) {
      n += __occ_aux32((uint64_t)p[0]<<32 | p[1], c);
  }
  // calculate Occ
  n += __occ_aux32(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
  if(c == 0) n -= ~k&31; // corrected for the masked bits
  */

  return n;
}

// an analogy to tmap_bwt_occ() but more efficient, requiring k <= l
inline void 
tmap_bwt_2occ(const tmap_bwt_t *bwt, uint32_t k, uint32_t l, uint8_t c, uint32_t *ok, uint32_t *ol)
{
  uint32_t _k, _l;
  if(k == l) {
      *ok = *ol = tmap_bwt_occ(bwt, k, c);
      return;
  }
  _k = (k >= bwt->primary)? k-1 : k;
  _l = (l >= bwt->primary)? l-1 : l;
  if(_l/bwt->occ_interval != _k/bwt->occ_interval || UINT32_MAX == k || UINT32_MAX == l) {
      if(l == UINT32_MAX) k = UINT32_MAX; 
      *ok = tmap_bwt_occ(bwt, k, c);
      *ol = tmap_bwt_occ(bwt, l, c);
  } else {
      uint32_t m, n, i, j;
      uint32_t *p;
      if(k >= bwt->primary) --k;
      if(l >= bwt->primary) --l;
      n = (p = tmap_bwt_occ_intv(bwt, k))[c];
      p += 4;
      // calculate *ok
      j = k >> 4 << 4; // divide by 16, then multiply by 16, to subtract k % 16.
      for(i = (k/bwt->occ_interval)*bwt->occ_interval; i < j; i += 16, p++) {
          n += __occ_aux16(p[0], c);
      }
      m = n; // save for ol
      n += __occ_aux16(p[0] & ~((1ul<<((~k&15)<<1)) - 1), c);
      if(c == 0) n -= ~k&15; // corrected for the masked bits
      *ok = n;
      // calculate *ol
      j = l >> 4 << 4;
      for(; i < j; i += 16, p += 1) {
          m += __occ_aux16(p[0], c);
      }
      m += __occ_aux16(p[0] & ~((1ul<<((~l&15)<<1)) - 1), c);
      if(c == 0) m -= ~l&15; // corrected for the masked bits
      *ol = m;
  }
}

#define __occ_aux4(bwt, b)											\
  ((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
   + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

inline void 
tmap_bwt_occ4(const tmap_bwt_t *bwt, uint32_t k, uint32_t cnt[4])
{
  uint32_t l, j, x;
  uint32_t *p;
  if(UINT32_MAX == k) {
      memset(cnt, 0, 4 * sizeof(uint32_t));
      return;
  }
  if(k >= bwt->primary) --k; // because $ is not in bwt
  p = tmap_bwt_occ_intv(bwt, k);
  memcpy(cnt, p, 16);
  p += 4; // move to the first bwt cell
  j = (k >> 4) << 4;
  for(l = (k / bwt->occ_interval) * bwt->occ_interval, x = 0; l < j; l += 16, ++p) {
    x += __occ_aux4(bwt, *p);
  }
  x += __occ_aux4(bwt, *p & ~((1U<<((~k&15)<<1)) - 1)) - (~k&15);
  cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

// an analogy to tmap_bwt_occ4() but more efficient, requiring k <= l
inline void 
tmap_bwt_2occ4(const tmap_bwt_t *bwt, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4])
{
  uint32_t _k, _l;
  if(k == l) {
      tmap_bwt_occ4(bwt, k, cntk);
      memcpy(cntl, cntk, 4 * sizeof(uint32_t));
      return;
  }
  _k = (k >= bwt->primary)? k-1 : k;
  _l = (l >= bwt->primary)? l-1 : l;
  if(_l/bwt->occ_interval != _k/bwt->occ_interval || UINT32_MAX == k || UINT32_MAX == l) {
      if(l == UINT32_MAX) k = UINT32_MAX; 
      tmap_bwt_occ4(bwt, k, cntk);
      tmap_bwt_occ4(bwt, l, cntl);
  } else {
      uint32_t i, j, x, y;
      uint32_t *p;
      int cl[4];
      if(k >= bwt->primary) --k; // because $ is not in bwt
      if(l >= bwt->primary) --l;
      cl[0] = cl[1] = cl[2] = cl[3] = 0;
      p = tmap_bwt_occ_intv(bwt, k);
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
tmap_bwt_pac2bwt_main(int argc, char *argv[])
{
  int c, is_large = 0, occ_interval = TMAP_BWT_OCC_INTERVAL, hash_width = TMAP_BWT_HASH_WIDTH, help = 0;

  while((c = getopt(argc, argv, "o:lw:vh")) >= 0) {
      switch(c) {
        case 'l': is_large = 1; break;
        case 'o': occ_interval = atoi(optarg); break;
        case 'w': hash_width = atoi(optarg); break;
        case 'v': tmap_progress_set_verbosity(1); break;
        case 'h': help = 1; break;
        default: return 1;
      }
  }
  if(1 != argc - optind || 1 == help) {
      tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [-l -o INT -w INT -v -h] <in.fasta>\n", PACKAGE, argv[0]);
      return 1;
  }
  if(occ_interval < 16 || 0 != (occ_interval % 16)) {
      tmap_error("option -o out of range", Exit, CommandLineArgument);
  }

  tmap_bwt_pac2bwt(argv[optind], is_large, occ_interval, hash_width);

  return 0;
}
