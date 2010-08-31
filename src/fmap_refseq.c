#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_seq.h"
#include "fmap_io.h"
#include "fmap_definitions.h"
#include "fmap_progress.h"
#include "fmap_refseq.h"

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

uint64_t
fmap_refseq_fasta2pac(const char *fn_fasta, int32_t compression)
{
  fmap_file_t *fp_fasta = NULL, *fp_pac = NULL, *fp_anno = NULL;;
  fmap_seq_t *seq = NULL;
  fmap_refseq_t *refseq = NULL;
  char *fn_pac = NULL, *fn_anno = NULL;
  uint8_t buffer[FMAP_REFSEQ_BUFFER_SIZE];
  int32_t i, l, buffer_length;
  uint8_t x = 0;
  uint64_t ref_len;

  fmap_progress_set_start_time(clock());
  fmap_progress_print("packing the reference FASTA");

  refseq = fmap_calloc(1, sizeof(fmap_refseq_t), "refseq");

  refseq->version_id = FMAP_VERSION_ID; 
  refseq->seed = FMAP_REFSEQ_SEED;
  srand48(refseq->seed);
  refseq->seq = buffer; // IMPORTANT: must nullify later
  refseq->annos = NULL;
  refseq->num_annos = 0;
  refseq->len = 0;
  refseq->is_rev = 0;
  memset(buffer, 0, FMAP_REFSEQ_BUFFER_SIZE);
  buffer_length = 0;

  // input files
  fp_fasta = fmap_file_fopen(fn_fasta, "rb", compression); 
  seq = fmap_seq_init(fp_fasta);

  // output files
  fn_pac = fmap_get_file_name(fn_fasta, FMAP_PAC_FILE);
  fp_pac = fmap_file_fopen(fn_pac, "wb", FMAP_PAC_COMPRESSION);

  // read in sequences
  while(0 <= (l = fmap_seq_read(seq))) {
      fmap_progress_print2("packing contig [%s:1-%d]", seq->name.s, l);

      refseq->num_annos++;
      refseq->annos = fmap_realloc(refseq->annos, sizeof(fmap_anno_t)*refseq->num_annos, "refseq->annos");
      refseq->annos[refseq->num_annos-1].name = fmap_strdup(seq->name.s);  // sequence name
      refseq->annos[refseq->num_annos-1].len = l;
      refseq->annos[refseq->num_annos-1].offset = (1 == refseq->num_annos) ? 0 : refseq->annos[refseq->num_annos-2].offset + refseq->annos[refseq->num_annos-2].len;

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
  ref_len = refseq->len; // save for return
  
  // write annotation file
  fn_anno = fmap_get_file_name(fn_fasta, FMAP_ANNO_FILE);
  fp_anno = fmap_file_fopen(fn_anno, "wb", FMAP_ANNO_COMPRESSION);
  fmap_refseq_write_anno(fp_anno, refseq); 

  // close files
  fmap_file_fclose(fp_fasta);
  fmap_file_fclose(fp_pac);
  fmap_file_fclose(fp_anno);

  fmap_refseq_destroy(refseq); 
  fmap_seq_destroy(seq);
  free(fn_pac);
  free(fn_anno);

  fmap_progress_print2("packed the reference FASTA");

  fmap_refseq_pac2revpac(fn_fasta);

  return ref_len;
}

void
fmap_refseq_pac2revpac(const char *fn_fasta)
{
  int32_t i, j, c;
  fmap_refseq_t *refseq=NULL, *refseq_rev=NULL;
  
  fmap_progress_print("reversing the packed reference FASTA");
  fmap_progress_set_start_time(clock());

  refseq = fmap_refseq_read(fn_fasta, 0);

  // shallow copy
  refseq_rev = fmap_calloc(1, sizeof(fmap_refseq_t), "refseq_rev");
  (*refseq_rev) = (*refseq);
  
  // update sequence
  refseq_rev->seq = NULL;
  refseq_rev->seq = fmap_calloc(fmap_refseq_seq_memory(refseq->len), sizeof(uint8_t), "refseq_rev->seq");
  for(i=refseq->len-1;0<=i;i--) {
      c = fmap_refseq_seq_i(refseq, i);
      j = refseq->len - i - 1;
      if(j < 0) {
          fprintf(stderr, "i=%d j=%d c=%d\n", i, j, c);
      }
      fmap_refseq_seq_store_i(refseq_rev, j, c);
  }

  // write
  fmap_refseq_write(refseq_rev, fn_fasta, 1);

  // free
  free(refseq_rev->seq);
  free(refseq_rev);
  fmap_refseq_destroy(refseq);
}

void 
fmap_refseq_write(fmap_refseq_t *refseq, const char *fn_fasta, uint32_t is_rev)
{
  fmap_file_t *fp_pac = NULL, *fp_anno = NULL;;
  char *fn_pac = NULL, *fn_anno = NULL;
  uint8_t x = 0;

  // write annotation file
  if(0 == is_rev) {
      fn_anno = fmap_get_file_name(fn_fasta, FMAP_ANNO_FILE);
      fp_anno = fmap_file_fopen(fn_anno, "wb", FMAP_ANNO_COMPRESSION);
      fmap_refseq_write_anno(fp_anno, refseq); 
      fmap_file_fclose(fp_anno);
      free(fn_anno);
  }

  // write the sequence
  fn_pac = fmap_get_file_name(fn_fasta, (0 == is_rev) ? FMAP_PAC_FILE : FMAP_REV_PAC_FILE);
  fp_pac = fmap_file_fopen(fn_pac, "wb", (0 == is_rev) ? FMAP_PAC_COMPRESSION : FMAP_REV_PAC_COMPRESSION);
  if(fmap_refseq_seq_memory(refseq->len) != fmap_file_fwrite(refseq->seq, sizeof(uint8_t), fmap_refseq_seq_memory(refseq->len), fp_pac)) {
      fmap_error(NULL, Exit, WriteFileError);
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
  fmap_file_fclose(fp_pac);
  free(fn_pac);
  
  fmap_progress_print2("reversed the packed reference FASTA");
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
  refseq->annos = fmap_calloc(refseq->num_annos, sizeof(fmap_anno_t), "refseq->annos"); // allocate memory
  for(i=0;i<refseq->num_annos;i++) { // read the annotations
      fmap_refseq_read_annos(fp, &refseq->annos[i]);
  }
}

fmap_refseq_t *
fmap_refseq_read(const char *fn_fasta, uint32_t is_rev)
{
  fmap_file_t *fp_pac = NULL, *fp_anno = NULL;;
  char *fn_pac = NULL, *fn_anno = NULL;
  fmap_refseq_t *refseq = NULL;
  
  // allocate some memory 
  refseq = fmap_calloc(1, sizeof(fmap_refseq_t), "refseq");
  refseq->is_rev = is_rev;

  // read annotation file
  fn_anno = fmap_get_file_name(fn_fasta, FMAP_ANNO_FILE);
  fp_anno = fmap_file_fopen(fn_anno, "rb", FMAP_ANNO_COMPRESSION);
  fmap_refseq_read_anno(fp_anno, refseq); 
  fmap_file_fclose(fp_anno);
  free(fn_anno);

  // read the sequence
  fn_pac = fmap_get_file_name(fn_fasta, (0 == is_rev) ? FMAP_PAC_FILE : FMAP_REV_PAC_FILE);
  fp_pac = fmap_file_fopen(fn_pac, "rb", (0 == is_rev) ? FMAP_PAC_COMPRESSION : FMAP_REV_PAC_COMPRESSION);
  refseq->seq = fmap_malloc(sizeof(uint8_t)*fmap_refseq_seq_memory(refseq->len), "refseq->seq"); // allocate
  if(fmap_refseq_seq_memory(refseq->len) 
     != fmap_file_fread(refseq->seq, sizeof(uint8_t), fmap_refseq_seq_memory(refseq->len), fp_pac)) {
      fmap_error(NULL, Exit, WriteFileError);
  }
  fmap_file_fclose(fp_pac);
  free(fn_pac);


  return refseq;
}

void
fmap_refseq_destroy(fmap_refseq_t *refseq)
{
  uint32_t i;

  for(i=0;i<refseq->num_annos;i++) {
      free(refseq->annos[i].name);
  }
  free(refseq->annos);
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

// zero-based
static inline int32_t
fmap_refseq_get_seqid1(fmap_refseq_t *refseq, uint32_t pacpos)
{
  int32_t left, right, mid;

  if(refseq->len <= pacpos) {
      fmap_error("Coordinate was larger than the reference", Exit, OutOfRange);
  }

  left = 0; mid = 0; right = refseq->num_annos;
  while (left < right) {
      mid = (left + right) >> 1;
      if(refseq->annos[mid].offset <= pacpos) {
          if(mid == refseq->num_annos - 1) break;
          if(pacpos < refseq->annos[mid+1].offset) break;
          left = mid + 1;
      } else right = mid;
  }

  return mid;
}

// zero-based
static inline int32_t
fmap_refseq_get_seqid(fmap_refseq_t *refseq, uint32_t pacpos, uint32_t aln_len)
{
  int32_t seqid_left, seqid_right;

  seqid_left = fmap_refseq_get_seqid1(refseq, pacpos);
  seqid_right = fmap_refseq_get_seqid1(refseq, pacpos+aln_len-1);

  if(seqid_left != seqid_right) return -1; // overlaps two chromosomes

  return seqid_left;
}

// zero-based
static inline uint32_t
fmap_refseq_get_pos(fmap_refseq_t *refseq, uint32_t pacpos, uint32_t seqid)
{
  return pacpos - refseq->annos[seqid].offset;
}

inline int32_t
fmap_refseq_pac2real(fmap_refseq_t *refseq, uint32_t pacpos, uint32_t aln_length, uint32_t *seqid, uint32_t *pos)
{
  (*seqid) = fmap_refseq_get_seqid(refseq, pacpos, aln_length);
  if((*seqid) < 0) return -1;
  (*pos) = fmap_refseq_get_pos(refseq, pacpos, (*seqid));

  return (*pos);
}

int
fmap_refseq_fasta2pac_main(int argc, char *argv[])
{
  if(argc < 2) {
      fprintf(stderr, "Usage: %s %s <in.fasta>\n", PACKAGE, argv[0]);
      return 1;
  }
  fmap_refseq_fasta2pac(argv[1], FMAP_FILE_NO_COMPRESSION);
  return 0;
}
