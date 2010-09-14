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

  h->flow = fmap_malloc(sizeof(char)*h->flow_length, "h->flow_sequence");
  h->key = fmap_malloc(sizeof(char)*h->key_length, "h->key_sequence");
  
  if(h->flow_length != fmap_file_fread(h->flow, sizeof(char), h->flow_length, fp)
     || h->key_length != fmap_file_fread(h->key, sizeof(char), h->key_length, fp)) {
      fmap_error("fmap_file_fread", Exit, ReadFileError);
  }
  n += sizeof(char)*(h->flow_length + h->key_length);

  if(h->gheader_length != n) {
      fmap_error("SFF global header did not match", Exit, ReadFileError);
  }

  fmap_sff_read_padding(fp, n);

  return h;
}

void
fmap_sff_header_destroy(fmap_sff_header_t *h)
{
  free(h->key);
  free(h->flow);
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

     rh->name = fmap_malloc(sizeof(char)*(1+rh->name_length), "rh->name");

     if(rh->name_length != fmap_file_fread(rh->name, sizeof(char), rh->name_length, fp)) {
         fmap_error("fmap_file_fread", Exit, ReadFileError);
     }
     n += sizeof(char)*rh->name_length;
     rh->name[rh->name_length]='\0';
  
     fmap_sff_read_padding(fp, n);

     return rh;
}

void
fmap_sff_read_header_destroy(fmap_sff_read_header_t *rh)
{
  free(rh->name);
  free(rh);
}

fmap_sff_read_t *
fmap_sff_read_read(fmap_file_t *fp, fmap_sff_header_t *gh, fmap_sff_read_header_t *rh)
{
  fmap_sff_read_t *r = NULL;
  uint32_t n = 0;

  r = fmap_calloc(1, sizeof(fmap_sff_read_t), "r");
  
  r->flowgram = fmap_malloc(sizeof(uint16_t)*gh->flow_length, "r->flowgram");
  r->flow_index = fmap_malloc(sizeof(uint8_t)*rh->n_bases, "r->flow_index");
  r->bases = fmap_malloc(sizeof(char)*rh->n_bases, "r->bases");
  r->quality = fmap_malloc(sizeof(uint8_t)*rh->n_bases, "r->quality");

  if(gh->flow_length != fmap_file_fread(r->flowgram, sizeof(uint16_t), gh->flow_length, fp)
     || rh->n_bases != fmap_file_fread(r->flow_index, sizeof(uint8_t), rh->n_bases, fp)
     || rh->n_bases != fmap_file_fread(r->bases, sizeof(char), rh->n_bases, fp)
     || rh->n_bases != fmap_file_fread(r->quality, sizeof(uint8_t), rh->n_bases, fp)) {
      fmap_error("fmap_file_fread", Exit, ReadFileError);
  }
     
  fmap_sff_read_padding(fp, n);

  return r;
}

void
fmap_sff_read_destroy(fmap_sff_read_t *r)
{
  free(r->flowgram);
  free(r->flow_index);
  free(r->bases);
  free(r->quality);
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
