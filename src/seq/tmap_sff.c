/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <string.h>
#include <netinet/in.h>
#include <config.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../io/tmap_file.h"
#include "tmap_sff.h"

static inline uint32_t
tmap_sff_read_padding(tmap_file_t *fp, uint32_t n)
{
  char padding[8]="\0";
  n = (n & 7); // (n % 8)
  if(0 != n) {
      n = 8 - n; // number of bytes of padding
      if(n != tmap_file_fread(padding, sizeof(char), n, fp)) {
          tmap_error("tmap_file_fread", Exit, ReadFileError);
      }
  }
  return n;
}

#ifdef TMAP_SFF_DEBUG
static void
tmap_sff_header_print(FILE *fp, tmap_sff_header_t *h)
{
  fprintf(stderr, "** SFF HEADER - START **\n");
  fprintf(stderr, "magic=%u\n", h->magic);
  fprintf(stderr, "version=%u\n", h->version);
  fprintf(stderr, "index_offset=%llu\n", h->index_offset);
  fprintf(stderr, "index_length=%u\n", h->index_length);
  fprintf(stderr, "n_reads=%u\n", h->n_reads);
  fprintf(stderr, "gheader_length=%u\n", h->gheader_length);
  fprintf(stderr, "key_length=%u\n", h->key_length);
  fprintf(stderr, "flow_length=%u\n", h->flow_length);
  fprintf(stderr, "flowgram_format=%u\n", h->flowgram_format);
  fprintf(stderr, "flow=%s\n", h->flow->s);
  fprintf(stderr, "key=%s\n", h->key->s);
  fprintf(stderr, "** SFF HEADER - END **\n");
}
#endif

tmap_sff_header_t *
tmap_sff_header_read(tmap_file_t *fp)
{
  tmap_sff_header_t *h = NULL;
  uint32_t n = 0;

  h = tmap_calloc(1, sizeof(tmap_sff_header_t), "h");

  if(1 != tmap_file_fread(&h->magic, sizeof(uint32_t), 1, fp)
     || 1 != tmap_file_fread(&h->version, sizeof(uint32_t), 1, fp)
     || 1 != tmap_file_fread(&h->index_offset, sizeof(uint64_t), 1, fp)
     || 1 != tmap_file_fread(&h->index_length, sizeof(uint32_t), 1, fp)
     || 1 != tmap_file_fread(&h->n_reads, sizeof(uint32_t), 1, fp)
     || 1 != tmap_file_fread(&h->gheader_length, sizeof(uint16_t), 1, fp)
     || 1 != tmap_file_fread(&h->key_length, sizeof(uint16_t), 1, fp)
     || 1 != tmap_file_fread(&h->flow_length, sizeof(uint16_t), 1, fp)
     || 1 != tmap_file_fread(&h->flowgram_format, sizeof(uint8_t), 1, fp)) {
      tmap_error("tmap_file_fread", Exit, ReadFileError);
  }
  n += 4*sizeof(uint32_t) + sizeof(uint64_t) + 3*sizeof(uint16_t) + sizeof(uint8_t);

  // convert values from big-endian
  h->magic = ntohl(h->magic);
  h->version = ntohl(h->version);
  h->index_offset = ntohll(h->index_offset);
  h->index_length = ntohl(h->index_length);
  h->n_reads = ntohl(h->n_reads);
  h->gheader_length = ntohs(h->gheader_length);
  h->key_length = ntohs(h->key_length);
  h->flow_length = ntohs(h->flow_length);

  if(TMAP_SFF_MAGIC != h->magic) {
      tmap_error("SFF magic number did not match", Exit, ReadFileError);
  }
  if(h->version != TMAP_SFF_VERSION) {
      tmap_error("SFF version number did not match", Exit, ReadFileError);
  }

  h->flow = tmap_string_init(h->flow_length+1);
  h->key = tmap_string_init(h->key_length+1);

  if(h->flow_length != tmap_file_fread(h->flow->s, sizeof(char), h->flow_length, fp)
     || h->key_length != tmap_file_fread(h->key->s, sizeof(char), h->key_length, fp)) {
      tmap_error("tmap_file_fread", Exit, ReadFileError);
  }
  n += sizeof(char)*(h->flow_length + h->key_length);

  // set the length and null-terminator
  h->flow->l = h->flow_length;
  h->key->l = h->key_length;
  h->flow->s[h->flow->l]='\0';
  h->key->s[h->key->l]='\0';

  n += tmap_sff_read_padding(fp, n);

#ifdef TMAP_SFF_DEBUG
  tmap_sff_header_print(stderr, h);
#endif

  if(h->gheader_length != n) {
      tmap_error("SFF global header length did not match", Exit, ReadFileError);
  }

  return h;
}

void
tmap_sff_header_destroy(tmap_sff_header_t *h)
{
  if(NULL == h) return;
  tmap_string_destroy(h->flow);
  tmap_string_destroy(h->key);
  free(h);
}

#ifdef TMAP_SFF_DEBUG
static void
tmap_sff_read_header_print(FILE *fp, tmap_sff_read_header_t *rh)
{
  fprintf(stderr, "** SFF READ HEADER - START **\n");
  fprintf(stderr, "rheader_length=%u\n", rh->rheader_length);
  fprintf(stderr, "name_length=%u\n", rh->name_length);
  fprintf(stderr, "n_bases=%u\n", rh->n_bases);
  fprintf(stderr, "clip_qual_left=%u\n", rh->clip_qual_left);
  fprintf(stderr, "clip_qual_right=%u\n", rh->clip_qual_right);
  fprintf(stderr, "clip_adapter_left=%u\n", rh->clip_adapter_left);
  fprintf(stderr, "clip_adapter_right=%u\n", rh->clip_adapter_right);
  fprintf(stderr, "name=%s\n", rh->name->s);
  fprintf(stderr, "** SFF READ HEADER - END **\n");
}
#endif

tmap_sff_read_header_t *
tmap_sff_read_header_read(tmap_file_t *fp)
{
  tmap_sff_read_header_t *rh = NULL;
  uint32_t n = 0;

  rh = tmap_calloc(1, sizeof(tmap_sff_read_header_t), "rh");

  if(1 != tmap_file_fread(&rh->rheader_length, sizeof(uint16_t), 1, fp)
     || 1 != tmap_file_fread(&rh->name_length, sizeof(uint16_t), 1, fp)
     || 1 != tmap_file_fread(&rh->n_bases, sizeof(uint32_t), 1, fp)
     || 1 != tmap_file_fread(&rh->clip_qual_left, sizeof(uint16_t), 1, fp)
     || 1 != tmap_file_fread(&rh->clip_qual_right, sizeof(uint16_t), 1, fp)
     || 1 != tmap_file_fread(&rh->clip_adapter_left, sizeof(uint16_t), 1, fp)
     || 1 != tmap_file_fread(&rh->clip_adapter_right, sizeof(uint16_t), 1, fp)) {
      tmap_error("tmap_file_fread", Exit, ReadFileError);
  }
  n += sizeof(uint32_t) + 6*sizeof(uint16_t);

  // convert values from big-endian
  rh->rheader_length = ntohs(rh->rheader_length);
  rh->name_length = ntohs(rh->name_length);
  rh->n_bases = ntohl(rh->n_bases);
  rh->clip_qual_left = ntohs(rh->clip_qual_left);
  rh->clip_qual_right = ntohs(rh->clip_qual_right);
  rh->clip_adapter_left = ntohs(rh->clip_adapter_left);
  rh->clip_adapter_right = ntohs(rh->clip_adapter_right);

  rh->name = tmap_string_init(rh->name_length+1);

  if(rh->name_length != tmap_file_fread(rh->name->s, sizeof(char), rh->name_length, fp)) {
      tmap_error("tmap_file_fread", Exit, ReadFileError);
  }
  n += sizeof(char)*rh->name_length;

  // set read name length and null-terminator
  rh->name->l = rh->name_length;
  rh->name->s[rh->name->l]='\0';

  n += tmap_sff_read_padding(fp, n);

#ifdef TMAP_SFF_DEBUG
  tmap_sff_read_header_print(stderr, rh);
#endif

  if(rh->rheader_length != n) {
      tmap_error("SFF read header length did not match", Exit, ReadFileError);
  }

  return rh;
}

void
tmap_sff_read_header_destroy(tmap_sff_read_header_t *rh)
{
  if(NULL == rh) return;
  tmap_string_destroy(rh->name);
  free(rh);
}

static tmap_sff_read_header_t *
tmap_sff_read_header_clone(tmap_sff_read_header_t *rh)
{
  tmap_sff_read_header_t *ret = NULL;

  ret = tmap_calloc(1, sizeof(tmap_sff_read_header_t), "rh");

  ret->rheader_length = rh->rheader_length; 
  ret->name_length = rh->name_length; 
  ret->n_bases = rh->n_bases; 
  ret->clip_qual_left = rh->clip_qual_left; 
  ret->clip_qual_right = rh->clip_qual_right;
  ret->clip_adapter_left = rh->clip_adapter_left; 
  ret->clip_adapter_right = rh->clip_adapter_right;
  ret->clip_left = ret->clip_left;
  ret->clip_right = ret->clip_right;
  ret->name = tmap_string_clone(rh->name);

  return ret;
}

#ifdef TMAP_SFF_DEBUG
static void
tmap_sff_read_print(FILE *fp, tmap_sff_read_t *r, tmap_sff_header_t *gh, tmap_sff_read_header_t *rh)
{
  uint32_t i;

  fprintf(stderr, "** SFF READ - START **\n");
  fprintf(stderr, "flowgram:\n");
  for(i=0;i<gh->flow_length;i++) {
      if(0 < i) fputc(',', stderr);
      fprintf(stderr, "%u", r->flowgram[i]); 
  }
  fputc('\n', stderr);
  fprintf(stderr, "flow_index:\n");
  for(i=0;i<rh->n_bases;i++) {
      if(0 < i) fputc(',', stderr);
      fprintf(stderr, "%u", r->flow_index[i]); 
  }
  fputc('\n', stderr);
  fprintf(stderr, "bases:\n");
  fprintf(stderr, "%s\n", r->bases->s);
  fprintf(stderr, "quality:\n");
  fprintf(stderr, "%s\n", r->quality->s);
  fprintf(stderr, "** SFF READ - END **\n");
}
#endif

tmap_sff_read_t *
tmap_sff_read_read(tmap_file_t *fp, tmap_sff_header_t *gh, tmap_sff_read_header_t *rh)
{
  tmap_sff_read_t *r = NULL;
  uint32_t i, n = 0;

  r = tmap_calloc(1, sizeof(tmap_sff_read_t), "r");

  r->flowgram = tmap_malloc(sizeof(uint16_t)*gh->flow_length, "r->flowgram");
  r->flow_index = tmap_malloc(sizeof(uint8_t)*rh->n_bases, "r->flow_index");

  r->bases = tmap_string_init(rh->n_bases+1);
  r->quality = tmap_string_init(rh->n_bases+1);

  if(gh->flow_length != tmap_file_fread(r->flowgram, sizeof(uint16_t), gh->flow_length, fp)
     || rh->n_bases != tmap_file_fread(r->flow_index, sizeof(uint8_t), rh->n_bases, fp)
     || rh->n_bases != tmap_file_fread(r->bases->s, sizeof(char), rh->n_bases, fp)
     || rh->n_bases != tmap_file_fread(r->quality->s, sizeof(char), rh->n_bases, fp)) {
      tmap_error("tmap_file_fread", Exit, ReadFileError);
  }
  n += sizeof(uint16_t)*gh->flow_length + 3*sizeof(uint8_t)*rh->n_bases;

  // set length and null-terminators
  r->bases->l = rh->n_bases;
  r->quality->l = rh->n_bases;
  r->bases->s[r->bases->l]='\0';
  r->quality->s[r->quality->l]='\0';

  // convert qualities from int to char
  for(i=0;i<r->quality->l;i++) {
      r->quality->s[i] = QUAL2CHAR(r->quality->s[i]);
  }

  // convert flowgram to host order
  for(i=0;i<gh->flow_length;i++) {
      r->flowgram[i] = ntohs(r->flowgram[i]);
  }

  n += tmap_sff_read_padding(fp, n);

#ifdef TMAP_SFF_DEBUG
  tmap_sff_read_print(stderr, r, gh, rh);
#endif

  return r;
}

void
tmap_sff_read_destroy(tmap_sff_read_t *r)
{
  if(NULL == r) return;
  free(r->flowgram);
  free(r->flow_index);
  tmap_string_destroy(r->bases);
  tmap_string_destroy(r->quality);
  free(r);

}

static tmap_sff_read_t *
tmap_sff_read_clone(tmap_sff_read_t *r, tmap_sff_header_t *gh, tmap_sff_read_header_t *rh)
{
  tmap_sff_read_t *ret = NULL;
  int32_t i;

  ret = tmap_calloc(1, sizeof(tmap_sff_read_t), "r");

  ret->flowgram = tmap_malloc(sizeof(uint16_t)*gh->flow_length, "ret->flowgram");
  for(i=0;i<gh->flow_length;i++) {
      ret->flowgram[i] = r->flowgram[i];
  }

  ret->flow_index = tmap_malloc(sizeof(uint8_t)*rh->n_bases, "ret->flow_index");
  for(i=0;i<rh->n_bases;i++) {
      ret->flow_index[i] = r->flow_index[i];
  }

  ret->bases = tmap_string_clone(r->bases);
  ret->quality = tmap_string_clone(r->quality);

  return ret;
}

tmap_sff_t *
tmap_sff_init()
{
  tmap_sff_t *sff = NULL;

  sff = tmap_calloc(1, sizeof(tmap_sff_t), "sff");
  sff->gheader = NULL;
  sff->rheader = NULL;
  sff->read = NULL;

  return sff;
}

void 
tmap_sff_destroy(tmap_sff_t *sff)
{
  if(NULL == sff) return;
  tmap_sff_read_header_destroy(sff->rheader);
  tmap_sff_read_destroy(sff->read);
  free(sff);
}

tmap_sff_t *
tmap_sff_clone(tmap_sff_t *sff)
{
  tmap_sff_t *ret = NULL;

  ret = tmap_sff_init();

  ret->gheader = sff->gheader;
  ret->rheader = tmap_sff_read_header_clone(sff->rheader);
  ret->read = tmap_sff_read_clone(sff->read, sff->gheader, sff->rheader);

  return ret;

}

void
tmap_sff_reverse(tmap_sff_t *sff)
{
  int32_t i;

  // reverse flowgram
  for(i=0;i<(sff->gheader->flow_length>>1);i++) {
      uint16_t tmp = sff->read->flowgram[sff->gheader->flow_length-1-i];
      sff->read->flowgram[sff->gheader->flow_length-1-i] = sff->read->flowgram[i];
      sff->read->flowgram[i] = tmp;
  }
  // reverse flow index
  for(i=0;i<(sff->rheader->n_bases>>1);i++) {
      uint8_t tmp = sff->read->flow_index[sff->rheader->n_bases-1-i];
      sff->read->flow_index[sff->rheader->n_bases-1-i] = sff->read->flow_index[i];
      sff->read->flow_index[i] = tmp;
  }
  // reverse compliment the bases
  tmap_string_reverse(sff->read->bases);
  // reverse the qualities
  tmap_string_reverse(sff->read->quality);
}

void
tmap_sff_reverse_compliment(tmap_sff_t *sff)
{
  int32_t i;

  // reverse flowgram
  for(i=0;i<(sff->gheader->flow_length>>1);i++) {
      uint16_t tmp = sff->read->flowgram[sff->gheader->flow_length-1-i];
      sff->read->flowgram[sff->gheader->flow_length-1-i] = sff->read->flowgram[i];
      sff->read->flowgram[i] = tmp;
  }
  // reverse flow index
  for(i=0;i<(sff->rheader->n_bases>>1);i++) {
      uint8_t tmp = sff->read->flow_index[sff->rheader->n_bases-1-i];
      sff->read->flow_index[sff->rheader->n_bases-1-i] = sff->read->flow_index[i];
      sff->read->flow_index[i] = tmp;
  }
  // reverse compliment the bases
  tmap_string_reverse_compliment(sff->read->bases, sff->is_int);
  // reverse the qualities
  tmap_string_reverse(sff->read->quality);
}

void
tmap_sff_compliment(tmap_sff_t *sff)
{
  // reverse compliment the bases
  tmap_string_compliment(sff->read->bases, sff->is_int);
}

void
tmap_sff_to_int(tmap_sff_t *sff)
{
  int32_t i;
  if(1 == sff->is_int) return;
  for(i=0;i<sff->read->bases->l;i++) {
      sff->read->bases->s[i] = tmap_nt_char_to_int[(int)sff->read->bases->s[i]];
  }
  sff->is_int = 1;
}

void
tmap_sff_to_char(tmap_sff_t *sff)
{
  int32_t i;
  if(0 == sff->is_int) return;
  for(i=0;i<sff->read->bases->l;i++) {
      sff->read->bases->s[i] = "ACGTN"[(int)sff->read->bases->s[i]];
  }
  sff->read->bases->s[sff->read->bases->l] = '\0';
  sff->is_int = 0;
}

inline tmap_string_t *
tmap_sff_get_bases(tmap_sff_t *sff)
{
  return sff->read->bases;
}

inline tmap_string_t *
tmap_sff_get_qualities(tmap_sff_t *sff)
{
  return sff->read->quality;
}

inline void
tmap_sff_remove_key_sequence(tmap_sff_t *sff, int32_t remove_clipping)
{
  int32_t i;

  int32_t left, right; // zero-based

  if(1 == remove_clipping) {
      // left clipping
      if(0 < sff->rheader->clip_adapter_left || 0 < sff->rheader->clip_qual_left) {
          if(sff->rheader->clip_adapter_left < sff->rheader->clip_qual_left) { // choose the largest
              sff->rheader->clip_left = sff->rheader->clip_qual_left - 1;
          }
          else {
              sff->rheader->clip_left = sff->rheader->clip_adapter_left - 1;
          }
          // do not include the key sequence 
          if(sff->gheader->key_length < sff->rheader->clip_qual_left) {
              left = sff->rheader->clip_left;
              sff->rheader->clip_left -= sff->gheader->key_length; 
          }
          else {
              left = sff->gheader->key_length; 
              sff->rheader->clip_left = 0;
          }
      }
      else {
          left = sff->gheader->key_length;
      }

      // right clipping
      if(0 < sff->rheader->clip_adapter_right || 0 < sff->rheader->clip_qual_right) {
          if(0 == sff->rheader->clip_qual_right 
             || (0 < sff->rheader->clip_adapter_right && sff->rheader->clip_adapter_right < sff->rheader->clip_qual_right)) {
              sff->rheader->clip_right = sff->rheader->n_bases - sff->rheader->clip_adapter_right;
              right = sff->rheader->clip_adapter_right-1;
          }
          else {
              sff->rheader->clip_right = sff->rheader->n_bases - sff->rheader->clip_qual_right;
              right = sff->rheader->clip_qual_right-1;
          }
      }
      else {
          right = sff->rheader->n_bases-1;
      }
  }
  else {
      left = sff->gheader->key_length;
      right = sff->rheader->n_bases-1;
  }

  // extract the bases
  for(i=0;i<right-left+1;i++) {
      sff->read->bases->s[i] = sff->read->bases->s[i+left];
      sff->read->quality->s[i] = sff->read->quality->s[i+left];
  }
  sff->read->bases->l = (right-left+1);
  sff->read->quality->l = (right-left+1);
  if(0 == sff->is_int) {
      sff->read->bases->s[sff->read->bases->l] = '\0';
      sff->read->quality->s[sff->read->quality->l] = '\0';
  }
}

int32_t
tmap_sff_get_flow_order_int(tmap_sff_t *sff, uint8_t **flow_order)
{
  int32_t i;
  int32_t flow_order_len = sff->gheader->flow->l;
  (*flow_order) = tmap_malloc(sizeof(uint8_t) * flow_order_len, "flow_order");
  for(i=0;i<flow_order_len;i++) {
      (*flow_order)[i] = tmap_nt_char_to_int[(int)sff->gheader->flow->s[i]];
  }
  return flow_order_len;
}

int32_t
tmap_sff_get_key_seq_int(tmap_sff_t *sff, uint8_t **key_seq)
{
  int32_t i;
  int32_t key_seq_len = sff->gheader->key->l;
  (*key_seq) = tmap_malloc(sizeof(uint8_t) * key_seq_len, "key_seq");
  for(i=0;i<key_seq_len;i++) {
      (*key_seq)[i] = tmap_nt_char_to_int[(int)sff->gheader->key->s[i]];
  }
  return key_seq_len;
}

int32_t
tmap_sff_get_flowgram(tmap_sff_t *sff, uint16_t **flowgram, int32_t mem)
{
  int32_t i;
  if(mem <= sff->gheader->flow_length) {
      (*flowgram) = tmap_realloc((*flowgram), sizeof(uint16_t) * sff->gheader->flow_length, "flowgram");
  }
  for(i=0;i<sff->gheader->flow_length;i++) {
      (*flowgram)[i] = sff->read->flowgram[i];
  }
  return sff->gheader->flow_length;
}
