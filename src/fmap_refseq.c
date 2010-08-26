#include <stdio.h>
#include <stdlib.h>

#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_seq.h"
#include "fmap_refseq.h"

// Input: ASCII character
// Output: 2-bit DNA value
uint8_t nt_char_to_int[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline char *
fmap_refseq_get_file_name(const char *prefix, int32_t is_anno)
{
  char *fn = NULL;
  if(0 == is_anno) {
      fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_REFSEQ_PAC_FILE_EXTENSION)), "fn"); 
      strcpy(fn, prefix);
      strcat(fn, FMAP_REFSEQ_PAC_FILE_EXTENSION);
  }
  else {
      fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_REFSEQ_ANNO_FILE_EXTENSION)), "fn"); 
      strcpy(fn, prefix);
      strcat(fn, FMAP_REFSEQ_ANNO_FILE_EXTENSION);
  }
  return fn;
}

static inline void 
fmap_refseq_write_header(fmap_file_t *fp, fmap_refseq_t *refseq)
{
  if(1 != fmap_file_fwrite(&refseq->version_id, sizeof(uint64_t), 1, fp) 
     || 1 != fmap_file_fwrite(&refseq->seed, sizeof(uint32_t), 1, fp) 
     || 1 != fmap_file_fwrite(&refseq->num_annos, sizeof(uint32_t), 1, fp)
     || 1 != fmap_file_fwrite(&refseq->len, sizeof(uint64_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
}

static inline void
fmap_refseq_write_annos(fmap_file_t *fp, fmap_anno_t *anno) 
{
  uint32_t len = strlen(anno->name)+1; // include null terminator

  if(1 != fmap_file_fwrite(&len, sizeof(uint32_t), 1, fp) 
     || len != fmap_file_fwrite(anno->name, sizeof(char), len, fp)
     || 1 != fmap_file_fwrite(&anno->len, sizeof(uint64_t), 1, fp)
     || 1 != fmap_file_fwrite(&anno->offset, sizeof(uint64_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
}

static inline void
fmap_refseq_write_anno(fmap_file_t *fp, fmap_refseq_t *refseq)
{
  uint32_t i;
  // write annotation file
  fmap_refseq_write_header(fp, refseq); // write the header
  for(i=0;i<refseq->num_annos;i++) { // write the annotations
      fmap_refseq_write_annos(fp, &refseq->annos[i]);
  }
}

fmap_refseq_t *
fmap_refseq_read_fasta(const char *fn_fasta, int32_t compression)
{
  fmap_file_t *fp_fasta = NULL, *fp_pac = NULL, *fp_anno = NULL;;
  fmap_seq_t *seq = NULL;
  fmap_refseq_t *refseq = NULL;
  char *fn_pac = NULL, *fn_anno = NULL;
  uint8_t buffer[FMAP_REFSEQ_BUFFER_SIZE];
  int32_t i, l, buffer_length;
  uint8_t x = 0;

  refseq = fmap_calloc(1, sizeof(fmap_refseq_t), "refseq");

  refseq->version_id = FMAP_REFSEQ_VERSION_ID; 
  refseq->seed = FMAP_REFSEQ_SEED;
  srand48(refseq->seed);
  refseq->seq = buffer; // IMPORTANT: must nullify later
  refseq->annos = NULL;
  refseq->num_annos = 0;
  refseq->len = 0;
  memset(buffer, 0, FMAP_REFSEQ_BUFFER_SIZE);
  buffer_length = 0;

  // input files
  fp_fasta = fmap_file_fopen(fn_fasta, "rb", compression); 
  seq = fmap_seq_init(fp_fasta);

  // output files
  fn_pac = fmap_refseq_get_file_name(fn_fasta, 0);
  fp_pac = fmap_file_fopen(fn_pac, "wb", FMAP_REFSEQ_COMPRESSION);
  fn_anno = fmap_refseq_get_file_name(fn_fasta, 1);
  fp_anno = fmap_file_fopen(fn_anno, "wb", FMAP_REFSEQ_COMPRESSION);

  // read in sequences
  while(0 <= (l = fmap_seq_read(seq))) {
      refseq->num_annos++;
      refseq->annos = fmap_realloc(refseq->annos, sizeof(fmap_anno_t)*refseq->num_annos, "refseq->annos");
      refseq->annos[refseq->num_annos-1].name = fmap_strdup(seq->name.s);  // sequence name
      refseq->annos[refseq->num_annos-1].len = l;
      refseq->annos[refseq->num_annos-1].offset = (1 == refseq->num_annos) ? 0 : refseq->annos[refseq->num_annos-2].offset + refseq->annos[refseq->num_annos-2].len;
      if(refseq->num_annos <= 1) {
          refseq->annos[refseq->num_annos-1].offset = 0; 
      }
      else {
          refseq->annos[refseq->num_annos-1].offset = refseq->annos[refseq->num_annos-2].offset 
            + refseq->annos[refseq->num_annos-2].len;
      }

      // fill the buffer
      for(i=0;i<l;i++) {
          int c = nt_char_to_int[(int)seq->seq.s[i]];
          if(4 <= c) c = lrand48() & 0x3; // random base
          if(buffer_length == (FMAP_REFSEQ_BUFFER_SIZE << 2)) { // 2-bit
              if(fmap_refseq_seq_memory(buffer_length) != fmap_file_fwrite(buffer, sizeof(uint8_t), fmap_refseq_seq_memory(buffer_length), fp_pac)) {
                  fmap_error(fn_pac, Exit, WriteFileError);
              }
              memset(buffer, 0, FMAP_REFSEQ_BUFFER_SIZE);
              buffer_length = 0;
          }
          fmap_refseq_seq_store_i(refseq, buffer_length, c);
          buffer_length++;
      }
      refseq->len += l;
  }
  // write out the buffer
  if(fmap_refseq_seq_memory(buffer_length) != fmap_file_fwrite(buffer, sizeof(uint8_t), fmap_refseq_seq_memory(buffer_length), fp_pac)) {
      fmap_error(fn_pac, Exit, WriteFileError);
  }
  if(refseq->len % 4 == 0) { // add an extra byte if we completely filled all bits
      if(1 != fmap_file_fwrite(&x, sizeof(uint8_t), 1, fp_pac)) {
          fmap_error(fn_pac, Exit, WriteFileError);
      }
  }
  // store number of unused bits at the last byte
  x = refseq->len % 4;
  if(1 != fmap_file_fwrite(&x, sizeof(uint8_t), 1, fp_pac)) {
      fmap_error(fn_pac, Exit, WriteFileError);
  }
  refseq->seq = NULL; // IMPORTANT: nullify this

  // write annotation file
  fmap_refseq_write_anno(fp_anno, refseq); 

  fmap_file_fclose(fp_fasta);
  fmap_file_fclose(fp_pac);
  fmap_file_fclose(fp_anno);

  fmap_seq_destroy(seq);
  free(fn_pac);
  free(fn_anno);

  return refseq;
}

void 
fmap_refseq_write(fmap_refseq_t *refseq, const char *prefix)
{
  fmap_file_t *fp_pac = NULL, *fp_anno = NULL;;
  char *fn_pac = NULL, *fn_anno = NULL;

  // write annotation file
  fn_anno = fmap_refseq_get_file_name(prefix, 1);
  fp_anno = fmap_file_fopen(fn_anno, "wb", FMAP_REFSEQ_COMPRESSION);
  fmap_refseq_write_anno(fp_anno, refseq); 
  fmap_file_fclose(fp_anno);
  free(fn_anno);

  // write the sequence
  fn_pac = fmap_refseq_get_file_name(prefix, 0);
  fp_pac = fmap_file_fopen(fn_pac, "wb", FMAP_REFSEQ_COMPRESSION);
  if(fmap_refseq_seq_memory(refseq->len) != fmap_file_fwrite(refseq->seq, sizeof(uint8_t), fmap_refseq_seq_memory(refseq->len), fp_pac)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
  fmap_file_fclose(fp_pac);
  free(fn_pac);
}

static inline void 
fmap_refseq_read_header(fmap_file_t *fp, fmap_refseq_t *refseq)
{
  if(1 != fmap_file_fread(&refseq->version_id, sizeof(uint64_t), 1, fp) 
     || 1 != fmap_file_fread(&refseq->seed, sizeof(uint32_t), 1, fp) 
     || 1 != fmap_file_fread(&refseq->num_annos, sizeof(uint32_t), 1, fp)
     || 1 != fmap_file_fread(&refseq->len, sizeof(uint64_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
}

static inline void
fmap_refseq_read_annos(fmap_file_t *fp, fmap_anno_t *anno) 
{
  uint32_t len = 0;

  if(1 != fmap_file_fread(&len, sizeof(uint32_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }

  anno->name = fmap_malloc(sizeof(char)*len, "anno->name");

  if(len != fmap_file_fread(anno->name, sizeof(char), len, fp)
     || 1 != fmap_file_fread(&anno->len, sizeof(uint64_t), 1, fp)
     || 1 != fmap_file_fread(&anno->offset, sizeof(uint64_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
}

static inline void
fmap_refseq_read_anno(fmap_file_t *fp, fmap_refseq_t *refseq)
{
  uint32_t i;
  // read annotation file
  fmap_refseq_read_header(fp, refseq); // read the header
  refseq->annos = fmap_malloc(sizeof(fmap_anno_t)*refseq->num_annos, "refseq->annos"); // allocate memory
  for(i=0;i<refseq->num_annos;i++) { // read the annotations
      fmap_refseq_read_annos(fp, &refseq->annos[i]);
  }
}

fmap_refseq_t *
fmap_refseq_read(const char *prefix)
{
  fmap_file_t *fp_pac = NULL, *fp_anno = NULL;;
  char *fn_pac = NULL, *fn_anno = NULL;
  fmap_refseq_t *refseq = NULL;

  // allocate some memory 
  refseq = fmap_calloc(1, sizeof(fmap_refseq_t), "refseq");

  // read annotation file
  fn_anno = fmap_refseq_get_file_name(prefix, 1);
  fp_anno = fmap_file_fopen(fn_anno, "rb", FMAP_REFSEQ_COMPRESSION);
  fmap_refseq_read_anno(fp_anno, refseq); 
  fmap_file_fclose(fp_anno);
  free(fn_anno);

  // read the sequence
  fn_pac = fmap_refseq_get_file_name(prefix, 0);
  fp_pac = fmap_file_fopen(fn_pac, "rb", FMAP_REFSEQ_COMPRESSION);
  refseq->seq = fmap_malloc(sizeof(char)*fmap_refseq_seq_memory(refseq->len), "refseq->seq"); // allocate
  if(fmap_refseq_seq_memory(refseq->len) != fmap_file_fread(refseq->seq, sizeof(uint8_t), fmap_refseq_seq_memory(refseq->len), fp_pac)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
  fmap_file_fclose(fp_pac);
  free(fn_pac);

  return refseq;
}

static void
fmap_refseq_destroy_annos(fmap_anno_t *anno)
{
  free(anno->name);
  free(anno);
}

void
fmap_refseq_destroy(fmap_refseq_t *refseq)
{
  uint32_t i;

  for(i=0;i<refseq->num_annos;i++) {
      fmap_refseq_destroy_annos(&refseq->annos[i]);
  }
  free(refseq->seq);
  free(refseq);
}

size_t
fmap_refseq_size(fmap_refseq_t *refseq)
{
  size_t size = 0;
  uint32_t i;

  size += sizeof(fmap_refseq_t);
  size += sizeof(char)*fmap_refseq_seq_memory(refseq->len);
  size += sizeof(fmap_refseq_t)*refseq->num_annos;
  for(i=0;i<refseq->num_annos;i++) {
      size += (1+strlen(refseq->annos[i].name))*sizeof(char); // include null terminator
  }

  return size;
}
