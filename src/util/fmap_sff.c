#include <stdlib.h>
#include <string.h>

#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_definitions.h"
#include "../io/fmap_file.h"
#include "fmap_sff.h"

static void
fmap_sff_read_padding(fmap_file_t *fp, uint32_t n)
{
  char padding[8]="\0";
  n = n % 8;
  if(0 != n) {
      if(n != fmap_file_fread(padding, sizeof(char), n, fp)) {
          fmap_error("fmap_file_fread", Exit, ReadFileError);
      }
  }
}

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

  if(FMAP_SFF_MAGIC != h->magic) {
      fmap_error("SFF magic number did not match", Exit, ReadFileError);
  }
  if(0 != memcmp(&h->version, FMAP_SFF_VERSION, sizeof(uint32_t))) {
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

  if(h->gheader_length != n) {
      fmap_error("SFF global header did not match", Exit, ReadFileError);
  }

  fmap_sff_read_padding(fp, n);

  return h;
}

void
fmap_sff_header_destroy(fmap_sff_header_t *h)
{
  fmap_string_destroy(h->flow);
  fmap_string_destroy(h->key);
  free(h);
}

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

     rh->name = fmap_string_init(rh->name_length+1);

     if(rh->name_length != fmap_file_fread(rh->name->s, sizeof(char), rh->name_length, fp)) {
         fmap_error("fmap_file_fread", Exit, ReadFileError);
     }
     n += sizeof(char)*rh->name_length;
     
     // set read name length and null-terminator
     rh->name->l = rh->name_length;
     rh->name->s[rh->name->l]='\0';
  
     fmap_sff_read_padding(fp, n);

     return rh;
}

void
fmap_sff_read_header_destroy(fmap_sff_read_header_t *rh)
{
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

fmap_sff_read_t *
fmap_sff_read_read(fmap_file_t *fp, fmap_sff_header_t *gh, fmap_sff_read_header_t *rh)
{
  fmap_sff_read_t *r = NULL;
  uint32_t n = 0;

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

  // set length and null-terminators
  r->bases->l = rh->n_bases;
  r->quality->l = rh->n_bases+1;
  r->bases->s[r->bases->l]='\0';
  r->quality->s[r->quality->l]='\0';
     
  fmap_sff_read_padding(fp, n);

  return r;
}

void
fmap_sff_read_destroy(fmap_sff_read_t *r)
{
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
  fmap_error("Not implemented", Exit, OutOfRange);
}

void
fmap_sff_to_int(fmap_sff_t *sff)
{
  fmap_error("Not implemented", Exit, OutOfRange);
}
