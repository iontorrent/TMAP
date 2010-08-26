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

fmap_refseq_t *
fmap_refseq_read_fasta(const char *fn_fasta, int32_t compression)
{
  fmap_file_t *fp = NULL;
  fmap_seq_t *seq = NULL;
  fmap_refseq_t *refseq = NULL;
  int32_t i, l;

  refseq = fmap_calloc(1, sizeof(fmap_refseq_t), "refseq");

  refseq->version_id = FMAP_REFSEQ_VERSION_ID; 
  refseq->seed = FMAP_REFSEQ_SEED;
  srand48(refseq->seed);
  refseq->seq = NULL;
  refseq->annos = NULL;
  refseq->num_annos = 0;
  refseq->len = 0;

  fp = fmap_file_fopen(fn_fasta, "rb", compression); 
  seq = fmap_seq_init(fp);

  while(0 <= (l = fmap_seq_read(seq))) {
      refseq->seq = fmap_realloc(refseq->seq, sizeof(char)*FMAP_REFSEQ_BYTE_LENGTH(l + refseq->len), "refseq->seq");

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

      for(i=0;i<l;i++) {
          if(0 == (i & 3)) { // every fourth
              refseq->seq[FMAP_REFSEQ_BYTE_LENGTH(refseq->len + i)] = 0; // set to zero
          }
          int c = nt_char_to_int[(int)seq->seq.s[l]];
          if(4 <= c) {
              c = lrand48() & 0x3; // random base
          }
          FMAP_REFSEQ_STORE(refseq, refseq->len+i, c);
      }
      refseq->len += l;
  }

  fmap_file_fclose(fp);

  fmap_seq_destroy(seq);

  return refseq;
}

static inline void 
fmap_refseq_write_header(fmap_file_t *fp, fmap_refseq_t *refseq)
{
  if(sizeof(uint64_t) != fmap_file_fwrite(&refseq->version_id, sizeof(uint64_t), 1, fp) 
     || sizeof(uint32_t) != fmap_file_fwrite(&refseq->seed, sizeof(uint32_t), 1, fp) 
     || sizeof(uint32_t) != fmap_file_fwrite(&refseq->num_annos, sizeof(uint32_t), 1, fp)
     || sizeof(uint64_t) != fmap_file_fwrite(&refseq->len, sizeof(uint64_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
}

static inline void
fmap_refseq_write_annos(fmap_file_t *fp, fmap_anno_t *anno) 
{
  uint32_t len = strlen(anno->name)+1; // include null terminator

  if(sizeof(uint32_t) != fmap_file_fwrite(&len, sizeof(uint32_t), 1, fp) 
     || len*sizeof(char) != fmap_file_fwrite(anno->name, sizeof(char), len, fp)
     || sizeof(uint64_t) != fmap_file_fwrite(&anno->len, sizeof(uint64_t), 1, fp)
     || sizeof(uint64_t) != fmap_file_fwrite(&anno->offset, sizeof(uint64_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
}

void 
fmap_refseq_write(fmap_refseq_t *refseq, const char *prefix)
{
  int32_t i;
  char *fn = NULL;
  fmap_file_t *fp = NULL;

  // creat the output file name
  fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_REFSEQ_FILE_EXTENSION)), "fn"); 
  strcpy(fn, prefix);
  strcat(fn, FMAP_REFSEQ_FILE_EXTENSION);

  // open the file
  fp = fmap_file_fopen(fn, "wb", FMAP_REFSEQ_COMPRESSION); 
  
  // free the output file name memory
  free(fn);

  // write the header
  fmap_refseq_write_header(fp, refseq);

  // write the sequence
  if(sizeof(char)*FMAP_REFSEQ_BYTE_LENGTH(refseq->len) != fmap_file_fwrite(refseq->seq, sizeof(char), FMAP_REFSEQ_BYTE_LENGTH(refseq->len), fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }

  // write the annotations
  for(i=0;i<refseq->num_annos;i++) {
      fmap_refseq_write_annos(fp, &refseq->annos[i]);
  }

  // close the file
  fmap_file_fclose(fp);

}

static inline void 
fmap_refseq_read_header(fmap_file_t *fp, fmap_refseq_t *refseq)
{
  if(sizeof(uint64_t) != fmap_file_fread(&refseq->version_id, sizeof(uint64_t), 1, fp) 
     || sizeof(uint32_t) != fmap_file_fread(&refseq->seed, sizeof(uint32_t), 1, fp) 
     || sizeof(uint32_t) != fmap_file_fread(&refseq->num_annos, sizeof(uint32_t), 1, fp)
     || sizeof(uint64_t) != fmap_file_fread(&refseq->len, sizeof(uint64_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
}

static inline void
fmap_refseq_read_annos(fmap_file_t *fp, fmap_anno_t *anno) 
{
  uint32_t len = 0;

  if(sizeof(uint32_t) != fmap_file_fread(&len, sizeof(uint32_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }

  anno->name = fmap_malloc(sizeof(char)*len, "anno->name");
     
  if(len*sizeof(char) != fmap_file_fread(anno->name, sizeof(char), len, fp)
     || sizeof(uint64_t) != fmap_file_fread(&anno->len, sizeof(uint64_t), 1, fp)
     || sizeof(uint64_t) != fmap_file_fread(&anno->offset, sizeof(uint64_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
}

fmap_refseq_t *
fmap_refseq_read(const char *prefix)
{
  int32_t i;
  char *fn = NULL;
  fmap_file_t *fp = NULL;
  fmap_refseq_t *refseq = NULL;

  // allocate some memory 
  refseq = fmap_calloc(1, sizeof(fmap_refseq_t), "refseq");

  // create the file name
  fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_REFSEQ_FILE_EXTENSION)), "fn"); 
  strcpy(fn, prefix);
  strcat(fn, FMAP_REFSEQ_FILE_EXTENSION);

  // open the file
  fp = fmap_file_fopen(fn, "rb", FMAP_REFSEQ_COMPRESSION); 

  // read the header
  fmap_refseq_read_header(fp, refseq);

  // allocate some memory
  refseq->seq = fmap_malloc(sizeof(char)*FMAP_REFSEQ_BYTE_LENGTH(refseq->len), "refseq->seq");
  refseq->annos = fmap_malloc(sizeof(fmap_anno_t)*refseq->num_annos, "refseq->annos");

  // read the sequence
  if(sizeof(char)*FMAP_REFSEQ_BYTE_LENGTH(refseq->len) != fmap_file_fread(refseq->seq, sizeof(char), FMAP_REFSEQ_BYTE_LENGTH(refseq->len), fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }

  // read the annotations
  for(i=0;i<refseq->num_annos;i++) {
      fmap_refseq_read_annos(fp, &refseq->annos[i]);
  }

  // close the file
  fmap_file_fclose(fp);

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
  size += sizeof(char)*FMAP_REFSEQ_BYTE_LENGTH(refseq->len);
  size += sizeof(fmap_refseq_t)*refseq->num_annos;
  for(i=0;i<refseq->num_annos;i++) {
      size += (1+strlen(refseq->annos[i].name))*sizeof(char); // include null terminator
  }

  return size;
}
