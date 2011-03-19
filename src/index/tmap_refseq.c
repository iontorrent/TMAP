/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <config.h>
#include <unistd.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_string.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_definitions.h"
#include "../seq/tmap_seq.h"
#include "../io/tmap_file.h"
#include "../io/tmap_seq_io.h"
#include "tmap_refseq.h"

static inline int32_t
tmap_refseq_supported(tmap_refseq_t *refseq)
{
  int32_t i, j, k;
  int32_t refseq_n, tmap_n;
  char *refseq_v = refseq->package_version->s;
  char *tmap_v = TMAP_REFSEQ_MIN_VERSION;

  // sanity check on version names
  for(i=j=0;i<strlen(refseq_v);i++) {
      if('.' == refseq_v[i]) j++;
  }
  if(2 != j) {
      tmap_error("did not find three version numbers", Exit, OutOfRange);
  }
  for(i=j=0;i<strlen(tmap_v);i++) {
      if('.' == tmap_v[i]) j++;
  }
  if(2 != j) {
      tmap_error("did not find three version numbers", Exit, OutOfRange);
  }

  for(i=j=k=0;i<3;i++) { // three version numbers
      refseq_n = atoi(refseq_v + j);
      tmap_n = atoi(tmap_v + k);
      if(refseq_n < tmap_n) {
          return 0;
      }
      if(i < 2) {
          while(refseq_v[j] != '.') {
              j++;
          }
          while(tmap_v[k] != '.') {
              k++;
          }
      }
  } 

  return 1;
}

static inline void 
tmap_refseq_write_header(tmap_file_t *fp, tmap_refseq_t *refseq)
{
  if(1 != tmap_file_fwrite(&refseq->version_id, sizeof(uint64_t), 1, fp) 
     || 1 != tmap_file_fwrite(&refseq->package_version->l, sizeof(size_t), 1, fp)
     || refseq->package_version->l+1 != tmap_file_fwrite(refseq->package_version->s, sizeof(char), refseq->package_version->l+1, fp)
     || 1 != tmap_file_fwrite(&refseq->seed, sizeof(uint32_t), 1, fp) 
     || 1 != tmap_file_fwrite(&refseq->num_annos, sizeof(uint32_t), 1, fp)
     || 1 != tmap_file_fwrite(&refseq->len, sizeof(uint64_t), 1, fp)) {
      tmap_error(NULL, Exit, WriteFileError);
  }
}

static inline void 
tmap_refseq_print_header(tmap_file_t *fp, tmap_refseq_t *refseq)
{
  uint32_t i;
  tmap_file_fprintf(fp, "version_id:\t%llu\n", (unsigned long long int)refseq->version_id);
  tmap_file_fprintf(fp, "package version:\t%s\n", refseq->package_version->s);
  tmap_file_fprintf(fp, "seed:\t%u\n", refseq->seed);
  for(i=0;i<refseq->num_annos;i++) {
      tmap_file_fprintf(fp, "contig-%d:\t%s\t%u\n", i+1, refseq->annos[i].name->s, refseq->annos[i].len);
  }
  tmap_file_fprintf(fp, "length:\t%llu\n", (unsigned long long int)refseq->len);
}

static inline void
tmap_refseq_write_annos(tmap_file_t *fp, tmap_anno_t *anno) 
{
  uint32_t len = anno->name->l+1; // include null terminator

  if(1 != tmap_file_fwrite(&len, sizeof(uint32_t), 1, fp) 
     || len != tmap_file_fwrite(anno->name->s, sizeof(char), len, fp)
     || 1 != tmap_file_fwrite(&anno->len, sizeof(uint64_t), 1, fp)
     || 1 != tmap_file_fwrite(&anno->offset, sizeof(uint64_t), 1, fp)
     || 1 != tmap_file_fwrite(&anno->num_amb, sizeof(uint32_t), 1, fp)) {
      tmap_error(NULL, Exit, WriteFileError);
  }
  if(0 < anno->num_amb) {
      if(1 != tmap_file_fwrite(anno->amb_positions, sizeof(uint32_t), anno->num_amb, fp)
         || 1 != tmap_file_fwrite(anno->amb_bases, sizeof(uint8_t), anno->num_amb, fp)) {
          tmap_error(NULL, Exit, WriteFileError);
      }
  }
}

static inline void
tmap_refseq_write_anno(tmap_file_t *fp, tmap_refseq_t *refseq)
{
  uint32_t i;

  // write annotation file
  tmap_refseq_write_header(fp, refseq); // write the header
  for(i=0;i<refseq->num_annos;i++) { // write the annotations
      tmap_refseq_write_annos(fp, &refseq->annos[i]);
  }
}

uint64_t
tmap_refseq_fasta2pac(const char *fn_fasta, int32_t compression)
{
  tmap_file_t *fp_pac = NULL, *fp_anno = NULL;
  tmap_seq_io_t *seqio = NULL;
  tmap_seq_t *seq = NULL;
  tmap_refseq_t *refseq = NULL;
  char *fn_pac = NULL, *fn_anno = NULL;
  uint8_t buffer[TMAP_REFSEQ_BUFFER_SIZE];
  int32_t i, j, l, buffer_length;
  uint32_t num_IUPAC_found= 0, amb_bases_mem = 0;
  uint8_t x = 0;
  uint64_t ref_len;

  tmap_progress_print("packing the reference FASTA");

  refseq = tmap_calloc(1, sizeof(tmap_refseq_t), "refseq");

  refseq->version_id = TMAP_VERSION_ID; 
  refseq->package_version = tmap_string_clone2(PACKAGE_VERSION);
  refseq->seed = TMAP_REFSEQ_SEED;
  srand48(refseq->seed);
  refseq->seq = buffer; // IMPORTANT: must nullify later
  refseq->annos = NULL;
  refseq->num_annos = 0;
  refseq->len = 0;
  refseq->is_rev = 0;
  refseq->is_shm = 0;
  memset(buffer, 0, TMAP_REFSEQ_BUFFER_SIZE);
  buffer_length = 0;

  // input files
  seqio = tmap_seq_io_init(fn_fasta, TMAP_SEQ_TYPE_FQ, 0, compression);
  seq = tmap_seq_init(TMAP_SEQ_TYPE_FQ);

  // output files
  fn_pac = tmap_get_file_name(fn_fasta, TMAP_PAC_FILE);
  fp_pac = tmap_file_fopen(fn_pac, "wb", TMAP_PAC_COMPRESSION);

  // read in sequences
  while(0 <= (l = tmap_seq_io_read(seqio, seq))) {
      tmap_progress_print2("packing contig [%s:1-%d]", seq->data.fq->name->s, l);

      refseq->num_annos++;
      refseq->annos = tmap_realloc(refseq->annos, sizeof(tmap_anno_t)*refseq->num_annos, "refseq->annos");
      refseq->annos[refseq->num_annos-1].name = tmap_string_clone(seq->data.fq->name); 
      refseq->annos[refseq->num_annos-1].len = l;
      refseq->annos[refseq->num_annos-1].offset = (1 == refseq->num_annos) ? 0 : refseq->annos[refseq->num_annos-2].offset + refseq->annos[refseq->num_annos-2].len;
      refseq->annos[refseq->num_annos-1].amb_positions = NULL;
      refseq->annos[refseq->num_annos-1].amb_bases = NULL;
      refseq->annos[refseq->num_annos-1].num_amb = 0;
      amb_bases_mem = 0;

      // fill the buffer
      for(i=0;i<l;i++) {
          uint8_t c = tmap_nt_char_to_int[(int)seq->data.fq->seq->s[i]];
          // handle IUPAC codes 
          if(4 <= c) {
              // change it to a mismatched base than the IUPAC code
              c = tmap_iupac_char_to_bit_string[(int)seq->data.fq->seq->s[i]];

              // warn users about IUPAC codes
              if(0 == num_IUPAC_found) { 
                  tmap_error("IUPAC codes were found and will be converted to random DNA bases", Warn, OutOfRange);
              }
              num_IUPAC_found++;

              // store IUPAC bases
              refseq->annos[refseq->num_annos-1].num_amb++;
              if(amb_bases_mem < refseq->annos[refseq->num_annos-1].num_amb) { // allocate more memory if necessary
                  amb_bases_mem = refseq->annos[refseq->num_annos-1].num_amb;
                  tmap_roundup32(amb_bases_mem);
                  refseq->annos[refseq->num_annos-1].amb_positions = tmap_realloc(refseq->annos[refseq->num_annos-1].amb_positions, sizeof(uint32_t) * amb_bases_mem, "refseq->annos[refseq->num_annos-1].amb_positions");
                  refseq->annos[refseq->num_annos-1].amb_bases = tmap_realloc(refseq->annos[refseq->num_annos-1].amb_bases, sizeof(uint8_t) * amb_bases_mem, "refseq->annos[refseq->num_annos-1].amb_bases");
              }
              refseq->annos[refseq->num_annos-1].amb_positions[refseq->annos[refseq->num_annos-1].num_amb-1] = i+1; // one-based
              refseq->annos[refseq->num_annos-1].amb_bases[refseq->annos[refseq->num_annos-1].num_amb-1] = tmap_iupac_char_to_int[(int)seq->data.fq->seq->s[i]];
              
              // randomize 
              if(c < 15) { // not an N
                  int32_t k, bases[4];
                  // how many bases does this not represent?
                  for(j=k=0;j<4;j++) {
                      if(!(c & (0x1 << j))) {
                          bases[k] = j;
                          k++;
                      }
                  } 
                  if(0 == k) {
                      tmap_error("bug encountered", Exit, OutOfRange);
                  }
                  // choose a random non-mismatch base
                  k = (uint8_t)(drand48() * k); // get the ith base
                  if(c & (0x1 << bases[k])) {
                      tmap_error("bug encountered", Exit, OutOfRange);
                  }
                  c = bases[k];
              }
              else {
                  c = lrand48() & 0x3; // random base
              }
          }
          if(3 < c) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }
          if(buffer_length == (TMAP_REFSEQ_BUFFER_SIZE << 2)) { // 2-bit
              if(tmap_refseq_seq_memory(buffer_length) != tmap_file_fwrite(buffer, sizeof(uint8_t), tmap_refseq_seq_memory(buffer_length), fp_pac)) {
                  tmap_error(fn_pac, Exit, WriteFileError);
              }
              memset(buffer, 0, TMAP_REFSEQ_BUFFER_SIZE);
              buffer_length = 0;
          }
          tmap_refseq_seq_store_i(refseq, buffer_length, c);
          buffer_length++;
      }
      refseq->len += l;
      // re-size the amibiguous bases
      if(refseq->annos[refseq->num_annos-1].num_amb < amb_bases_mem) {
          amb_bases_mem = refseq->annos[refseq->num_annos-1].num_amb;
          refseq->annos[refseq->num_annos-1].amb_positions = tmap_realloc(refseq->annos[refseq->num_annos-1].amb_positions, sizeof(uint32_t) * amb_bases_mem, "refseq->annos[refseq->num_annos-1].amb_positions");
          refseq->annos[refseq->num_annos-1].amb_bases = tmap_realloc(refseq->annos[refseq->num_annos-1].amb_bases, sizeof(uint8_t) * amb_bases_mem, "refseq->annos[refseq->num_annos-1].amb_bases");
      }
  }
  // write out the buffer
  if(tmap_refseq_seq_memory(buffer_length) != tmap_file_fwrite(buffer, sizeof(uint8_t), tmap_refseq_seq_memory(buffer_length), fp_pac)) {
      tmap_error(fn_pac, Exit, WriteFileError);
  }
  if(refseq->len % 4 == 0) { // add an extra byte if we completely filled all bits
      if(1 != tmap_file_fwrite(&x, sizeof(uint8_t), 1, fp_pac)) {
          tmap_error(fn_pac, Exit, WriteFileError);
      }
  }
  // store number of unused bits at the last byte
  x = refseq->len % 4;
  if(1 != tmap_file_fwrite(&x, sizeof(uint8_t), 1, fp_pac)) {
      tmap_error(fn_pac, Exit, WriteFileError);
  }
  refseq->seq = NULL; // IMPORTANT: nullify this
  ref_len = refseq->len; // save for return
      
  tmap_progress_print2("total genome length [%u]", refseq->len);
  if(0 < num_IUPAC_found) {
      if(1 == num_IUPAC_found) {
          tmap_progress_print("%u IUPAC base was found and converted to random DNA base", num_IUPAC_found);
      }
      else {
          tmap_progress_print("%u IUPAC bases were found and converted to random DNA bases", num_IUPAC_found);
      }
  }

  // write annotation file
  fn_anno = tmap_get_file_name(fn_fasta, TMAP_ANNO_FILE);
  fp_anno = tmap_file_fopen(fn_anno, "wb", TMAP_ANNO_COMPRESSION);
  tmap_refseq_write_anno(fp_anno, refseq); 

  // close files
  tmap_file_fclose(fp_pac);
  tmap_file_fclose(fp_anno);

  // check sequence name uniqueness
  for(i=0;i<refseq->num_annos;i++) {
      for(j=i+1;j<refseq->num_annos;j++) {
          if(0 == strcmp(refseq->annos[i].name->s, refseq->annos[j].name->s)) {
              tmap_file_fprintf(tmap_file_stderr, "Contigs have the same name: #%d [%s] and #%d [%s]\n",
                                i+1, refseq->annos[i].name->s, 
                                j+1, refseq->annos[j].name->s); 
              tmap_error("Contig names must be unique", Exit, OutOfRange);
          }
      }
  }

  tmap_refseq_destroy(refseq); 
  tmap_seq_io_destroy(seqio);
  tmap_seq_destroy(seq);
  free(fn_pac);
  free(fn_anno);

  tmap_progress_print2("packed the reference FASTA");

  tmap_refseq_pac2revpac(fn_fasta);

  return ref_len;
}

void
tmap_refseq_pac2revpac(const char *fn_fasta)
{
  uint32_t i, j, c;
  tmap_refseq_t *refseq=NULL, *refseq_rev=NULL;

  tmap_progress_print("reversing the packed reference FASTA");

  refseq = tmap_refseq_read(fn_fasta, 0);

  // shallow copy
  refseq_rev = tmap_calloc(1, sizeof(tmap_refseq_t), "refseq_rev");
  (*refseq_rev) = (*refseq);

  // update sequence
  refseq_rev->seq = NULL;
  refseq_rev->seq = tmap_calloc(tmap_refseq_seq_memory(refseq->len), sizeof(uint8_t), "refseq_rev->seq");
  for(i=0;i<refseq->len;i++) {
      c = tmap_refseq_seq_i(refseq, i);
      j = refseq->len - i - 1;
      tmap_refseq_seq_store_i(refseq_rev, j, c);
  }

  // write
  tmap_refseq_write(refseq_rev, fn_fasta, 1);

  // free
  free(refseq_rev->seq);
  free(refseq_rev);
  tmap_refseq_destroy(refseq);

  tmap_progress_print2("reversed the packed reference FASTA");
}

void 
tmap_refseq_write(tmap_refseq_t *refseq, const char *fn_fasta, uint32_t is_rev)
{
  tmap_file_t *fp_pac = NULL, *fp_anno = NULL;
  char *fn_pac = NULL, *fn_anno = NULL;
  uint8_t x = 0;

  // write annotation file
  if(0 == is_rev) {
      fn_anno = tmap_get_file_name(fn_fasta, TMAP_ANNO_FILE);
      fp_anno = tmap_file_fopen(fn_anno, "wb", TMAP_ANNO_COMPRESSION);
      tmap_refseq_write_anno(fp_anno, refseq); 
      tmap_file_fclose(fp_anno);
      free(fn_anno);
  }

  // write the sequence
  fn_pac = tmap_get_file_name(fn_fasta, (0 == is_rev) ? TMAP_PAC_FILE : TMAP_REV_PAC_FILE);
  fp_pac = tmap_file_fopen(fn_pac, "wb", (0 == is_rev) ? TMAP_PAC_COMPRESSION : TMAP_REV_PAC_COMPRESSION);
  if(tmap_refseq_seq_memory(refseq->len) != tmap_file_fwrite(refseq->seq, sizeof(uint8_t), tmap_refseq_seq_memory(refseq->len), fp_pac)) {
      tmap_error(NULL, Exit, WriteFileError);
  }
  if(refseq->len % 4 == 0) { // add an extra byte if we completely filled all bits
      if(1 != tmap_file_fwrite(&x, sizeof(uint8_t), 1, fp_pac)) {
          tmap_error(fn_pac, Exit, WriteFileError);
      }
  }
  // store number of unused bits at the last byte
  x = refseq->len % 4;
  if(1 != tmap_file_fwrite(&x, sizeof(uint8_t), 1, fp_pac)) {
      tmap_error(fn_pac, Exit, WriteFileError);
  }
  tmap_file_fclose(fp_pac);
  free(fn_pac);
}

static inline void 
tmap_refseq_read_header(tmap_file_t *fp, tmap_refseq_t *refseq)
{
  size_t package_version_l;
  if(1 != tmap_file_fread(&refseq->version_id, sizeof(uint64_t), 1, fp) 
     || 1 != tmap_file_fread(&package_version_l, sizeof(size_t), 1, fp)) {
      tmap_error(NULL, Exit, ReadFileError);
  }
  if(refseq->version_id != TMAP_VERSION_ID) {
      tmap_error("version id did not match", Exit, ReadFileError);
  }

  refseq->package_version = tmap_string_init(package_version_l+1); // add one for the null terminator
  refseq->package_version->l = package_version_l;
  if(refseq->package_version->l+1 != tmap_file_fread(refseq->package_version->s, sizeof(char), refseq->package_version->l+1, fp)) {
      tmap_error(NULL, Exit, ReadFileError);
  }
  if(0 == tmap_refseq_supported(refseq)) {
      tmap_error("the reference index is not supported", Exit, ReadFileError);
  }
     
  if(1 != tmap_file_fread(&refseq->seed, sizeof(uint32_t), 1, fp) 
     || 1 != tmap_file_fread(&refseq->num_annos, sizeof(uint32_t), 1, fp)
     || 1 != tmap_file_fread(&refseq->len, sizeof(uint64_t), 1, fp)) {
      tmap_error(NULL, Exit, ReadFileError);
  }

}

static inline void
tmap_refseq_read_annos(tmap_file_t *fp, tmap_anno_t *anno) 
{
  uint32_t len = 0; // includes the null-terminator
  
  if(1 != tmap_file_fread(&len, sizeof(uint32_t), 1, fp)) {
      tmap_error(NULL, Exit, ReadFileError);
  }

  anno->name = tmap_string_init(len);

  if(len != tmap_file_fread(anno->name->s, sizeof(char), len, fp)
     || 1 != tmap_file_fread(&anno->len, sizeof(uint64_t), 1, fp)
     || 1 != tmap_file_fread(&anno->offset, sizeof(uint64_t), 1, fp)
     || 1 != tmap_file_fread(&anno->num_amb, sizeof(uint32_t), 1, fp)) {
      tmap_error(NULL, Exit, ReadFileError);
  }
  if(0 < anno->num_amb) {
      anno->amb_positions = tmap_malloc(sizeof(uint32_t) * anno->num_amb, "anno->amb_positions");
      anno->amb_bases = tmap_malloc(sizeof(uint32_t) * anno->num_amb, "anno->amb_bases");
      if(1 != tmap_file_fread(anno->amb_positions, sizeof(uint32_t), anno->num_amb, fp)
         || 1 != tmap_file_fread(anno->amb_bases, sizeof(uint8_t), anno->num_amb, fp)) {
          tmap_error(NULL, Exit, WriteFileError);
      }
  }
  else {
      anno->amb_positions = NULL;
      anno->amb_bases = NULL;
  }
  // set name length
  anno->name->l = len-1;
}

static inline void
tmap_refseq_read_anno(tmap_file_t *fp, tmap_refseq_t *refseq)
{
  uint32_t i;
  // read annotation file
  tmap_refseq_read_header(fp, refseq); // read the header
  refseq->annos = tmap_calloc(refseq->num_annos, sizeof(tmap_anno_t), "refseq->annos"); // allocate memory
  for(i=0;i<refseq->num_annos;i++) { // read the annotations
      tmap_refseq_read_annos(fp, &refseq->annos[i]);
  }
}

tmap_refseq_t *
tmap_refseq_read(const char *fn_fasta, uint32_t is_rev)
{
  tmap_file_t *fp_pac = NULL, *fp_anno = NULL;
  char *fn_pac = NULL, *fn_anno = NULL;
  tmap_refseq_t *refseq = NULL;

  // allocate some memory 
  refseq = tmap_calloc(1, sizeof(tmap_refseq_t), "refseq");
  refseq->is_rev = is_rev;
  refseq->is_shm = 0;

  // read annotation file
  fn_anno = tmap_get_file_name(fn_fasta, TMAP_ANNO_FILE);
  fp_anno = tmap_file_fopen(fn_anno, "rb", TMAP_ANNO_COMPRESSION);
  tmap_refseq_read_anno(fp_anno, refseq); 
  tmap_file_fclose(fp_anno);
  free(fn_anno);

  // read the sequence
  fn_pac = tmap_get_file_name(fn_fasta, (0 == is_rev) ? TMAP_PAC_FILE : TMAP_REV_PAC_FILE);
  fp_pac = tmap_file_fopen(fn_pac, "rb", (0 == is_rev) ? TMAP_PAC_COMPRESSION : TMAP_REV_PAC_COMPRESSION);
  refseq->seq = tmap_malloc(sizeof(uint8_t)*tmap_refseq_seq_memory(refseq->len), "refseq->seq"); // allocate
  if(tmap_refseq_seq_memory(refseq->len) 
     != tmap_file_fread(refseq->seq, sizeof(uint8_t), tmap_refseq_seq_memory(refseq->len), fp_pac)) {
      tmap_error(NULL, Exit, ReadFileError);
  }
  tmap_file_fclose(fp_pac);
  free(fn_pac);


  return refseq;
}

size_t
tmap_refseq_shm_num_bytes(tmap_refseq_t *refseq)
{
  // returns the number of bytes to allocate for shared memory
  int32_t i;
  size_t n = 0;

  n += sizeof(uint64_t); // version_id
  n += sizeof(size_t); // package_version->l
  n += sizeof(uint32_t); // seed
  n += sizeof(uint32_t); // annos
  n += sizeof(uint64_t); // len
  n += sizeof(uint32_t); // is_rev
  n += sizeof(char)*(refseq->package_version->l+1); // package_version->s
  n += sizeof(uint8_t)*tmap_refseq_seq_memory(refseq->len); // seq 
  for(i=0;i<refseq->num_annos;i++) {
      n += sizeof(uint64_t); // len
      n += sizeof(uint64_t); // offset
      n += sizeof(size_t); // annos[i].name->l
      n += sizeof(uint32_t); // annos[i].num_amb
      n += sizeof(char)*(refseq->annos[i].name->l+1); // annos[i].name->s
      n += sizeof(uint32_t)*refseq->annos[i].num_amb; // amb_positions
      n += sizeof(uint8_t)*refseq->annos[i].num_amb; // amb_bases
  }

  return n;
}

size_t
tmap_refseq_shm_read_num_bytes(const char *fn_fasta, uint32_t is_rev)
{
  size_t n = 0;
  tmap_file_t *fp_anno = NULL;
  char *fn_anno = NULL;
  tmap_refseq_t *refseq = NULL;

  // allocate some memory 
  refseq = tmap_calloc(1, sizeof(tmap_refseq_t), "refseq");
  refseq->is_rev = is_rev;
  refseq->is_shm = 0;

  // read the annotation file
  fn_anno = tmap_get_file_name(fn_fasta, TMAP_ANNO_FILE);
  fp_anno = tmap_file_fopen(fn_anno, "rb", TMAP_ANNO_COMPRESSION);
  tmap_refseq_read_anno(fp_anno, refseq);
  tmap_file_fclose(fp_anno);
  free(fn_anno);

  // No need to read in the pac
  refseq->seq = NULL;

  // get the number of bytes
  n = tmap_refseq_shm_num_bytes(refseq);

  // destroy
  tmap_refseq_destroy(refseq);

  return n;
}

uint8_t *
tmap_refseq_shm_pack(tmap_refseq_t *refseq, uint8_t *buf)
{
  int32_t i;

  // fixed length data
  memcpy(buf, &refseq->version_id, sizeof(uint64_t)); buf += sizeof(uint64_t);
  memcpy(buf, &refseq->package_version->l, sizeof(size_t)); buf += sizeof(size_t);
  memcpy(buf, &refseq->seed, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(buf, &refseq->num_annos, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(buf, &refseq->len, sizeof(uint64_t)); buf += sizeof(uint64_t);
  memcpy(buf, &refseq->is_rev, sizeof(uint32_t)); buf += sizeof(uint32_t);
  // variable length data
  memcpy(buf, refseq->package_version->s, sizeof(char)*(refseq->package_version->l+1));
  buf += sizeof(char)*(refseq->package_version->l+1); 
  memcpy(buf, refseq->seq, tmap_refseq_seq_memory(refseq->len)*sizeof(uint8_t)); 
  buf += tmap_refseq_seq_memory(refseq->len)*sizeof(uint8_t);

  for(i=0;i<refseq->num_annos;i++) {
      // fixed length data
      memcpy(buf, &refseq->annos[i].len, sizeof(uint64_t)); buf += sizeof(uint64_t); 
      memcpy(buf, &refseq->annos[i].offset, sizeof(uint64_t)); buf += sizeof(uint64_t); 
      memcpy(buf, &refseq->annos[i].name->l, sizeof(size_t)); buf += sizeof(size_t);
      memcpy(buf, &refseq->annos[i].num_amb, sizeof(uint32_t)); buf += sizeof(uint32_t);
      // variable length data
      memcpy(buf, refseq->annos[i].name->s, sizeof(char)*(refseq->annos[i].name->l+1));
      buf += sizeof(char)*(refseq->annos[i].name->l+1);
      if(0 < refseq->annos[i].num_amb) {
          memcpy(buf, refseq->annos[i].amb_positions, sizeof(uint32_t)*refseq->annos[i].num_amb);
          buf += sizeof(uint32_t)*refseq->annos[i].num_amb;
          memcpy(buf, refseq->annos[i].amb_bases, sizeof(uint8_t)*refseq->annos[i].num_amb);
          buf += sizeof(uint8_t)*refseq->annos[i].num_amb;
      }
  }

  return buf;
}

tmap_refseq_t *
tmap_refseq_shm_unpack(uint8_t *buf)
{
  int32_t i;
  tmap_refseq_t *refseq = NULL;
  
  if(NULL == buf) return NULL;

  refseq = tmap_calloc(1, sizeof(tmap_refseq_t), "refseq");

  // fixed length data
  memcpy(&refseq->version_id, buf, sizeof(uint64_t)) ; buf += sizeof(uint64_t);
  if(refseq->version_id != TMAP_VERSION_ID) {
      tmap_error("version id did not match", Exit, ReadFileError);
  }
      
  refseq->package_version = tmap_string_init(0);
  memcpy(&refseq->package_version->l, buf, sizeof(size_t)); buf += sizeof(size_t);
  memcpy(&refseq->seed, buf, sizeof(uint32_t)) ; buf += sizeof(uint32_t);
  memcpy(&refseq->num_annos, buf, sizeof(uint32_t)) ; buf += sizeof(uint32_t);
  memcpy(&refseq->len, buf, sizeof(uint64_t)) ; buf += sizeof(uint64_t);
  memcpy(&refseq->is_rev, buf, sizeof(uint32_t)) ; buf += sizeof(uint32_t);

  // variable length data
  refseq->package_version->s = (char*)buf;
  refseq->package_version->m = refseq->package_version->l+1;
  buf += sizeof(char)*(refseq->package_version->l+1); 
  if(0 == tmap_refseq_supported(refseq)) {
      tmap_error("the reference index is not supported", Exit, ReadFileError);
  }
  refseq->seq = (uint8_t*)buf;
  buf += tmap_refseq_seq_memory(refseq->len)*sizeof(uint8_t);
  refseq->annos = tmap_calloc(refseq->num_annos, sizeof(tmap_anno_t), "refseq->annos");
  for(i=0;i<refseq->num_annos;i++) {
      // fixed length data
      memcpy(&refseq->annos[i].len, buf, sizeof(uint64_t)); buf += sizeof(uint64_t); 
      memcpy(&refseq->annos[i].offset, buf, sizeof(uint64_t)); buf += sizeof(uint64_t); 
      refseq->annos[i].name = tmap_string_init(0);
      memcpy(&refseq->annos[i].name->l, buf, sizeof(size_t)); buf += sizeof(size_t);
      refseq->annos[i].name->m = refseq->annos[i].name->l+1;
      memcpy(&refseq->annos[i].num_amb, buf, sizeof(uint32_t)); buf += sizeof(uint32_t);
      // variable length data
      refseq->annos[i].name->s = (char*)buf;
      buf += sizeof(char)*refseq->annos[i].name->l+1;
      if(0 < refseq->annos[i].num_amb) {
          refseq->annos[i].amb_positions = (uint32_t*)buf;
          buf += sizeof(uint32_t)*refseq->annos[i].num_amb;
          refseq->annos[i].amb_bases = (uint8_t*)buf;
          buf += sizeof(uint8_t)*refseq->annos[i].num_amb;
      }
      else {
          refseq->annos[i].amb_positions = NULL;
          refseq->annos[i].amb_bases = NULL;
      }
  }

  refseq->is_shm = 1;

  return refseq;
}

void
tmap_refseq_destroy(tmap_refseq_t *refseq)
{
  uint32_t i;

  if(1 == refseq->is_shm) {
      free(refseq->package_version);
      for(i=0;i<refseq->num_annos;i++) {
          free(refseq->annos[i].name);
      }
      free(refseq->annos);
      free(refseq);
  }
  else {
      tmap_string_destroy(refseq->package_version);
      for(i=0;i<refseq->num_annos;i++) {
          tmap_string_destroy(refseq->annos[i].name);
          free(refseq->annos[i].amb_positions);
          free(refseq->annos[i].amb_bases);
      }
      free(refseq->annos);
      free(refseq->seq);
      free(refseq);
  }
}

// zero-based
static inline int32_t
tmap_refseq_get_seqid1(const tmap_refseq_t *refseq, uint32_t pacpos)
{
  int32_t left, right, mid;

  if(refseq->len < pacpos) {
      tmap_error("Coordinate was larger than the reference", Exit, OutOfRange);
  }

  left = 0; mid = 0; right = refseq->num_annos;
  while (left < right) {
      mid = (left + right) >> 1;
      if(refseq->annos[mid].offset < pacpos) {
          if(mid == refseq->num_annos - 1) break;
          if(pacpos <= refseq->annos[mid+1].offset) break;
          left = mid + 1;
      } else right = mid;
  }

  if(refseq->num_annos < mid) {
      return refseq->num_annos;
  }

  return mid;
}

// zero-based
static inline int32_t
tmap_refseq_get_seqid(const tmap_refseq_t *refseq, uint32_t pacpos, uint32_t aln_len)
{
  int32_t seqid_left, seqid_right;

  seqid_left = tmap_refseq_get_seqid1(refseq, pacpos);
  if(1 == aln_len) return seqid_left;

  seqid_right = tmap_refseq_get_seqid1(refseq, pacpos+aln_len-1);
  if(seqid_left != seqid_right) return -1; // overlaps two chromosomes

  return seqid_left;
}

// zero-based
static inline uint32_t
tmap_refseq_get_pos(const tmap_refseq_t *refseq, uint32_t pacpos, uint32_t seqid)
{
  // note: offset is zero-based
  return pacpos - refseq->annos[seqid].offset;
}

inline uint32_t
tmap_refseq_pac2real(const tmap_refseq_t *refseq, uint32_t pacpos, uint32_t aln_length, uint32_t *seqid, uint32_t *pos)
{
  (*seqid) = tmap_refseq_get_seqid(refseq, pacpos, aln_length);
  if((*seqid) == (uint32_t)-1) {
      (*pos) = (uint32_t)-1;
      return 0;
  }
  (*pos) = tmap_refseq_get_pos(refseq, pacpos, (*seqid));

  return (*pos);
}

inline int32_t
tmap_refseq_subseq(const tmap_refseq_t *refseq, uint32_t pacpos, uint32_t length, uint8_t *target)
{
  uint32_t k, l, pacpos_upper;
  if(0 == length) {
      return 0;
  }
  else if(pacpos + length - 1 < refseq->len) {
      pacpos_upper = pacpos + length - 1;
  }
  else {
      pacpos_upper = refseq->len;
  }
  for(k=pacpos,l=0;k<=pacpos_upper;k++,l++) {
      // k-1 since pacpos is one-based
      target[l] = tmap_refseq_seq_i(refseq, k-1);
  }
  return l;
}

inline int32_t
tmap_refseq_amb_bases(const tmap_refseq_t *refseq, uint32_t seqid, uint32_t start, uint32_t end)
{
  int64_t low, high, mid;
  if(0 == refseq->annos[seqid-1].num_amb) {
      return 0;
  }
  else if(1 == refseq->annos[seqid-1].num_amb) {
      if(start <= refseq->annos[seqid-1].amb_positions[0] &&
         refseq->annos[seqid-1].amb_positions[0] <= end) {
          return 1;
      }
      else {
          return 0;
      }
  }

  low = 0; 
  high = refseq->annos[seqid-1].num_amb - 1;
  while(low <= high) {
      mid = (low + high) / 2;
      if(refseq->annos[seqid-1].amb_positions[mid] == start) {
          return mid+1;
      }
      else if(refseq->annos[seqid-1].amb_positions[mid] < start) {
          if(end <= refseq->annos[seqid-1].amb_positions[mid]) {
              return mid+1;
          }
          low = mid + 1;
      }
      else if(start < refseq->annos[seqid-1].amb_positions[mid]) {
          if(refseq->annos[seqid-1].amb_positions[mid] <= end) {
              return mid+1;
          }
          high = mid - 1;
      }
  }
  return 0;
}

int
tmap_refseq_fasta2pac_main(int argc, char *argv[])
{
  int c, help=0;

  while((c = getopt(argc, argv, "vh")) >= 0) {
      switch(c) {
        case 'v': tmap_progress_set_verbosity(1); break;
        case 'h': help = 1; break;
        default: return 1;
      }
  }
  if(1 != argc - optind || 1 == help) {
      tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [-vh] <in.fasta>\n", PACKAGE, argv[0]);
      return 1;
  }

  tmap_refseq_fasta2pac(argv[optind], TMAP_FILE_NO_COMPRESSION);

  return 0;
}

int
tmap_refseq_refinfo_main(int argc, char *argv[])
{
  int c, help=0;
  tmap_refseq_t *refseq = NULL;
  tmap_file_t *fp_anno = NULL;
  char *fn_anno = NULL;
  char *fn_fasta = NULL;

  while((c = getopt(argc, argv, "vh")) >= 0) {
      switch(c) {
        case 'v': tmap_progress_set_verbosity(1); break;
        case 'h': help = 1; break;
        default: return 1;
      }
  }
  if(1 != argc - optind || 1 == help) {
      tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [-vh] <in.fasta>\n", PACKAGE, argv[0]);
      return 1;
  }

  fn_fasta = argv[optind];

  // Note: 'tmap_file_stdout' should not have been previously modified
  tmap_file_stdout = tmap_file_fdopen(fileno(stdout), "wb", TMAP_FILE_NO_COMPRESSION);

  // allocate some memory 
  refseq = tmap_calloc(1, sizeof(tmap_refseq_t), "refseq");
  refseq->is_rev = 0;
  refseq->is_shm = 0;

  // read the annotation file
  fn_anno = tmap_get_file_name(fn_fasta, TMAP_ANNO_FILE);
  fp_anno = tmap_file_fopen(fn_anno, "rb", TMAP_ANNO_COMPRESSION);
  tmap_refseq_read_anno(fp_anno, refseq);
  tmap_file_fclose(fp_anno);
  free(fn_anno);
  
  // no need to read in the pac
  refseq->seq = NULL;

  // print the header
  tmap_refseq_print_header(tmap_file_stdout, refseq);

  // destroy
  tmap_refseq_destroy(refseq);

  // close the output
  tmap_file_fclose(tmap_file_stdout);

  return 0;
}

int
tmap_refseq_pac2fasta_main(int argc, char *argv[])
{
  int c, help=0;
  uint32_t i, j;
  char *fn_fasta = NULL;
  tmap_refseq_t *refseq = NULL;

  while((c = getopt(argc, argv, "vh")) >= 0) {
      switch(c) {
        case 'v': tmap_progress_set_verbosity(1); break;
        case 'h': help = 1; break;
        default: return 1;
      }
  }
  if(1 != argc - optind || 1 == help) {
      tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [-vh] <in.fasta>\n", PACKAGE, argv[0]);
      return 1;
  }

  fn_fasta = argv[optind];

  // Note: 'tmap_file_stdout' should not have been previously modified
  tmap_file_stdout = tmap_file_fdopen(fileno(stdout), "wb", TMAP_FILE_NO_COMPRESSION);

  // read in the reference sequence
  refseq = tmap_refseq_read(fn_fasta, 0);

  for(i=0;i<refseq->num_annos;i++) {
      tmap_file_fprintf(tmap_file_stdout, ">%s", refseq->annos[i].name->s); // new line handled later
      for(j=0;j<refseq->annos[i].len;j++) {
          if(0 == (j % TMAP_REFSEQ_FASTA_LINE_LENGTH)) {
              tmap_file_fprintf(tmap_file_stdout, "\n");
          }
          tmap_file_fprintf(tmap_file_stdout, "%c", "ACGTN"[(int)tmap_refseq_seq_i(refseq, j + refseq->annos[i].offset)]);
      }
      tmap_file_fprintf(tmap_file_stdout, "\n");
  }

  // destroy
  tmap_refseq_destroy(refseq);

  // close the output
  tmap_file_fclose(tmap_file_stdout);

  return 0;
}
