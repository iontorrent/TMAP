#include <stdlib.h>
#include <string.h>
#include <netinet/in.h>

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "../io/fmap_file.h"
#include "fmap_sff.h"

static inline uint32_t
fmap_sff_read_padding(fmap_file_t *fp, uint32_t n)
{
  char padding[8]="\0";
  n = (n & 7); // (n % 8)
  if(0 != n) {
      n = 8 - n; // number of bytes of padding
      if(n != fmap_file_fread(padding, sizeof(char), n, fp)) {
          fmap_error("fmap_file_fread", Exit, ReadFileError);
      }
  }
  return n;
}

#ifdef FMAP_SFF_DEBUG
static void
fmap_sff_header_print(FILE *fp, fmap_sff_header_t *h)
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

fmap_sff_header_t *
fmap_sff_header_read(fmap_file_t *fp)
{
  fmap_sff_header_t *h = NULL;
  uint32_t n = 0;

  h = fmap_calloc(1, sizeof(fmap_sff_header_t), "h");

  if(1 != fmap_file_fread(&h->magic, sizeof(uint32_t), 1, fp)
     || 1 != fmap_file_fread(&h->version, sizeof(uint32_t), 1, fp)
     || 1 != fmap_file_fread(&h->index_offset, sizeof(uint64_t), 1, fp)
     || 1 != fmap_file_fread(&h->index_length, sizeof(uint32_t), 1, fp)
     || 1 != fmap_file_fread(&h->n_reads, sizeof(uint32_t), 1, fp)
     || 1 != fmap_file_fread(&h->gheader_length, sizeof(uint16_t), 1, fp)
     || 1 != fmap_file_fread(&h->key_length, sizeof(uint16_t), 1, fp)
     || 1 != fmap_file_fread(&h->flow_length, sizeof(uint16_t), 1, fp)
     || 1 != fmap_file_fread(&h->flowgram_format, sizeof(uint8_t), 1, fp)) {
      fmap_error("fmap_file_fread", Exit, ReadFileError);
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

  if(FMAP_SFF_MAGIC != h->magic) {
      fmap_error("SFF magic number did not match", Exit, ReadFileError);
  }
  if(h->version != FMAP_SFF_VERSION) {
      fmap_error("SFF version number did not match", Exit, ReadFileError);
  }

  h->flow = fmap_string_init(h->flow_length+1);
  h->key = fmap_string_init(h->key_length+1);

  if(h->flow_length != fmap_file_fread(h->flow->s, sizeof(char), h->flow_length, fp)
     || h->key_length != fmap_file_fread(h->key->s, sizeof(char), h->key_length, fp)) {
      fmap_error("fmap_file_fread", Exit, ReadFileError);
  }
  n += sizeof(char)*(h->flow_length + h->key_length);

  // set the length and null-terminator
  h->flow->l = h->flow_length;
  h->key->l = h->key_length;
  h->flow->s[h->flow->l]='\0';
  h->key->s[h->key->l]='\0';

  n += fmap_sff_read_padding(fp, n);

#ifdef FMAP_SFF_DEBUG
  fmap_sff_header_print(stderr, h);
#endif

  if(h->gheader_length != n) {
      fmap_error("SFF global header length did not match", Exit, ReadFileError);
  }

  return h;
}

void
fmap_sff_header_destroy(fmap_sff_header_t *h)
{
  if(NULL == h) return;
  fmap_string_destroy(h->flow);
  fmap_string_destroy(h->key);
  free(h);
}

#ifdef FMAP_SFF_DEBUG
static void
fmap_sff_read_header_print(FILE *fp, fmap_sff_read_header_t *rh)
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

fmap_sff_read_header_t *
fmap_sff_read_header_read(fmap_file_t *fp)
{
  fmap_sff_read_header_t *rh = NULL;
  uint32_t n = 0;

  rh = fmap_calloc(1, sizeof(fmap_sff_read_header_t), "rh");

  if(1 != fmap_file_fread(&rh->rheader_length, sizeof(uint16_t), 1, fp)
     || 1 != fmap_file_fread(&rh->name_length, sizeof(uint16_t), 1, fp)
     || 1 != fmap_file_fread(&rh->n_bases, sizeof(uint32_t), 1, fp)
     || 1 != fmap_file_fread(&rh->clip_qual_left, sizeof(uint16_t), 1, fp)
     || 1 != fmap_file_fread(&rh->clip_qual_right, sizeof(uint16_t), 1, fp)
     || 1 != fmap_file_fread(&rh->clip_adapter_left, sizeof(uint16_t), 1, fp)
     || 1 != fmap_file_fread(&rh->clip_adapter_right, sizeof(uint16_t), 1, fp)) {
      fmap_error("fmap_file_fread", Exit, ReadFileError);
  }
  n += sizeof(uint32_t) + 6*sizeof(uint16_t);

  // convert values from big-endian
  rh->rheader_length = ntohs(rh->rheader_length);
  rh->name_length = ntohs(rh->name_length);
  rh->n_bases = ntohl(rh->n_bases);
  rh->clip_qual_left = ntohs(rh->clip_qual_left);
  rh->clip_qual_right = ntohs(rh->clip_qual_left);
  rh->clip_adapter_left = ntohs(rh->clip_adapter_left);
  rh->clip_adapter_right = ntohs(rh->clip_adapter_left);

  rh->name = fmap_string_init(rh->name_length+1);

  if(rh->name_length != fmap_file_fread(rh->name->s, sizeof(char), rh->name_length, fp)) {
      fmap_error("fmap_file_fread", Exit, ReadFileError);
  }
  n += sizeof(char)*rh->name_length;

  // set read name length and null-terminator
  rh->name->l = rh->name_length;
  rh->name->s[rh->name->l]='\0';

  n += fmap_sff_read_padding(fp, n);

#ifdef FMAP_SFF_DEBUG
  fmap_sff_read_header_print(stderr, rh);
#endif

  if(rh->rheader_length != n) {
      fmap_error("SFF read header length did not match", Exit, ReadFileError);
  }

  return rh;
}

void
fmap_sff_read_header_destroy(fmap_sff_read_header_t *rh)
{
  if(NULL == rh) return;
  fmap_string_destroy(rh->name);
  free(rh);
}

static fmap_sff_read_header_t *
fmap_sff_read_header_clone(fmap_sff_read_header_t *rh)
{
  fmap_sff_read_header_t *ret = NULL;

  ret = fmap_calloc(1, sizeof(fmap_sff_read_header_t), "rh");

  ret->rheader_length = rh->rheader_length; 
  ret->name_length = rh->name_length; 
  ret->n_bases = rh->n_bases; 
  ret->clip_qual_left = rh->clip_qual_left; 
  ret->clip_qual_right = rh->clip_qual_left; 
  ret->clip_adapter_left = rh->clip_adapter_left; 
  ret->clip_adapter_right = rh->clip_adapter_left; 
  ret->name = fmap_string_clone(rh->name);

  return ret;
}

#ifdef FMAP_SFF_DEBUG
static void
fmap_sff_read_print(FILE *fp, fmap_sff_read_t *r, fmap_sff_header_t *gh, fmap_sff_read_header_t *rh)
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

fmap_sff_read_t *
fmap_sff_read_read(fmap_file_t *fp, fmap_sff_header_t *gh, fmap_sff_read_header_t *rh)
{
  fmap_sff_read_t *r = NULL;
  uint32_t i, n = 0;

  r = fmap_calloc(1, sizeof(fmap_sff_read_t), "r");

  r->flowgram = fmap_malloc(sizeof(uint16_t)*gh->flow_length, "r->flowgram");
  r->flow_index = fmap_malloc(sizeof(uint8_t)*rh->n_bases, "r->flow_index");

  r->bases = fmap_string_init(rh->n_bases+1);
  r->quality = fmap_string_init(rh->n_bases+1);

  if(gh->flow_length != fmap_file_fread(r->flowgram, sizeof(uint16_t), gh->flow_length, fp)
     || rh->n_bases != fmap_file_fread(r->flow_index, sizeof(uint8_t), rh->n_bases, fp)
     || rh->n_bases != fmap_file_fread(r->bases->s, sizeof(char), rh->n_bases, fp)
     || rh->n_bases != fmap_file_fread(r->quality->s, sizeof(char), rh->n_bases, fp)) {
      fmap_error("fmap_file_fread", Exit, ReadFileError);
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

  n += fmap_sff_read_padding(fp, n);

#ifdef FMAP_SFF_DEBUG
  fmap_sff_read_print(stderr, r, gh, rh);
#endif

  return r;
}

void
fmap_sff_read_destroy(fmap_sff_read_t *r)
{
  if(NULL == r) return;
  free(r->flowgram);
  free(r->flow_index);
  fmap_string_destroy(r->bases);
  fmap_string_destroy(r->quality);
  free(r);

}

static fmap_sff_read_t *
fmap_sff_read_clone(fmap_sff_read_t *r, fmap_sff_header_t *gh, fmap_sff_read_header_t *rh)
{
  fmap_sff_read_t *ret = NULL;
  int32_t i;

  ret = fmap_calloc(1, sizeof(fmap_sff_read_t), "r");

  ret->flowgram = fmap_malloc(sizeof(uint16_t)*gh->flow_length, "ret->flowgram");
  for(i=0;i<gh->flow_length;i++) {
      ret->flowgram[i] = r->flowgram[i];
  }

  ret->flow_index = fmap_malloc(sizeof(uint8_t)*rh->n_bases, "ret->flow_index");
  for(i=0;i<rh->n_bases;i++) {
      ret->flow_index[i] = r->flow_index[i];
  }

  ret->bases = fmap_string_clone(r->bases);
  ret->quality = fmap_string_clone(r->quality);

  return ret;
}

fmap_sff_t *
fmap_sff_init()
{
  fmap_sff_t *sff = NULL;

  sff = fmap_calloc(1, sizeof(fmap_sff_t), "sff");
  sff->gheader = NULL;
  sff->rheader = NULL;
  sff->read = NULL;

  return sff;
}

void 
fmap_sff_destroy(fmap_sff_t *sff)
{
  if(NULL == sff) return;
  fmap_sff_read_header_destroy(sff->rheader);
  fmap_sff_read_destroy(sff->read);
  free(sff);
}

fmap_sff_t *
fmap_sff_clone(fmap_sff_t *sff)
{
  fmap_sff_t *ret = NULL;

  ret = fmap_sff_init();

  ret->gheader = sff->gheader;
  ret->rheader = fmap_sff_read_header_clone(sff->rheader);
  ret->read = fmap_sff_read_clone(sff->read, sff->gheader, sff->rheader);

  return ret;

}

void
fmap_sff_reverse_compliment(fmap_sff_t *sff)
{
  int32_t i;

  // reverse flow index
  for(i=0;i<(sff->gheader->flow_length>>2);i++) {
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
  fmap_string_reverse_compliment(sff->read->bases, sff->is_int);
  // reverse the qualities
  fmap_string_reverse(sff->read->quality);
}

void
fmap_sff_to_int(fmap_sff_t *sff)
{
  int32_t i;
  if(1 == sff->is_int) return;
  for(i=0;i<sff->read->bases->l;i++) {
      sff->read->bases->s[i] = nt_char_to_int[(int)sff->read->bases->s[i]];
  }
  sff->is_int = 1;
}

void
fmap_sff_to_char(fmap_sff_t *sff)
{
  int32_t i;
  if(0 == sff->is_int) return;
  for(i=0;i<sff->read->bases->l;i++) {
      sff->read->bases->s[i] = "ACGTN"[(int)sff->read->bases->s[i]];
  }
  sff->is_int = 0;
}

inline fmap_string_t *
fmap_sff_get_bases(fmap_sff_t *sff)
{
  return sff->read->bases;
}

inline fmap_string_t *
fmap_sff_get_qualities(fmap_sff_t *sff)
{
  return sff->read->quality;
}

inline void
fmap_sff_remove_key_sequence(fmap_sff_t *sff)
{
  int32_t i;
  // remove the key sequence
  for(i=0;i<sff->rheader->n_bases - sff->gheader->key_length;i++) { 
      sff->read->bases->s[i] = sff->read->bases->s[i+sff->gheader->key_length];
      sff->read->quality->s[i] = sff->read->quality->s[i+sff->gheader->key_length];
  }
  sff->read->bases->l -= sff->gheader->key_length;
  sff->read->quality->l -= sff->gheader->key_length;
}
