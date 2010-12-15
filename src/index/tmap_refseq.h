/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_REFSEQ_H_
#define TMAP_REFSEQ_H_

#include <stdint.h>
#include "../util/tmap_string.h"
#include "../io/tmap_file.h"

/*! 
  DNA Reference Sequence Library
  */

// seed for our random number generator
#define TMAP_REFSEQ_SEED 13
// buffer size for reading in from a FASTA file
#define TMAP_REFSEQ_BUFFER_SIZE 0x10000

/*! 
  @param  _len  the number of bases stored 
  @return       the number of bytes allocated
  */
#define tmap_refseq_seq_memory(_len) ((((_len)-1) >> 2) + 1)

/*! 
  @param  _i  the 0-based base position 
  @return     the 0-based byte index
  */
#define tmap_refseq_seq_byte_i(_i) ((_i) >> 2)

/*! 
  @param  _i  the 0-based base position 
  @return     the number of bits the base is shifted (returns a multiple of two)
  @details    the reference is stored in a 2-bit format
  */
#define tmap_refseq_seq_byte_shift(_i) ((0x3 - ((_i) & 0x3)) << 1)

/*! 
  @param  _refseq  pointer to the structure holding the reference sequence
  @param  _i       the 0-based base position to store
  @param  _c       the base's 2-bit integer representation
  */
#define tmap_refseq_seq_store_i(_refseq, _i, _c) (_refseq->seq[tmap_refseq_seq_byte_i(_i)] |= _c << tmap_refseq_seq_byte_shift(_i))

/*! 
  @param  _refseq  pointer to the structure holding the reference sequence
  @param  _i       the 0-based base position to retrieve
  @return          the base's 2-bit integer representation
  */
#define tmap_refseq_seq_i(_refseq, _i) ((_refseq->seq[tmap_refseq_seq_byte_i(_i)] >> tmap_refseq_seq_byte_shift(_i)) & 0x3)

/*! 
  */
typedef struct {
    tmap_string_t *name;  /*!< the name of the contig */
    uint64_t len;  /*!< the length of the current contig  */
    uint64_t offset;  /*!< the offset from the start of the reference (zero-based) */
} tmap_anno_t;

/*! 
  */
typedef struct {
    uint64_t version_id;  /*!< the version id of this file */
    uint32_t seed;  /*!< the random base generator seed */
    uint8_t *seq;  /*!< the packed nucleotide sequence, with contigs concatenated */
    tmap_anno_t *annos;  /*!< the annotations about the contigs */
    uint32_t num_annos;  /*!< the number of contigs (and annotations) */
    uint64_t len;  /*!< the total length of the reference sequence */
    uint32_t is_rev;  /*!< 1 if the reference sequence was reversed, 0 otherwise */
    uint32_t is_shm;  /*!< 1 if loaded from shared memory, 0 otherwise */
} tmap_refseq_t;

/*! 
  @param  fn_fasta     the file name of the fasta file
  @param  compression  the type of compression, if any to be used
  @return              the length of the reference sequence
  */
uint64_t
tmap_refseq_fasta2pac(const char *fn_fasta, int32_t compression);

/*! 
  @param  fn_fasta     the file name of the fasta file
  */
void
tmap_refseq_pac2revpac(const char *fn_fasta);

/*! 
  @param  refseq    pointer to the structure in which to store the data 
  @param  fn_fasta  the fn_fasta of the file to be written, usually the fasta file name 
  @param  is_rev    0 if to write the reverse packed sequence, 1 otherwise
  @details          this will not overwrite the annotation file if "is_rev" is 1
  */
void
tmap_refseq_write(tmap_refseq_t *refseq, const char *fn_fasta, uint32_t is_rev);

/*! 
  @param  fn_fasta  the fn_fasta of the file to be read, usually the fasta file name 
  @param  is_rev    0 if to read the reverse packed sequence, 1 otherwise
  @return           a pointer to the initialized memory
  */
tmap_refseq_t *
tmap_refseq_read(const char *fn_fasta, uint32_t is_rev);

/*! 
  @param  refseq  the refseq structure 
  @return         the number of bytes required for this bwt in shared memory
  */
size_t
tmap_refseq_shm_num_bytes(tmap_refseq_t *refseq);

/*! 
  @param  fn_fasta  the fn_fasta of the file to be read, usually the fasta file name 
  @param  is_rev    0 if to read the reverse packed sequence, 1 otherwise
  @return           the number of bytes required for this bwt in shared memory
  */
size_t
tmap_refseq_shm_read_num_bytes(const char *fn_fasta, uint32_t is_rev);

/*! 
  @param  refseq  the refseq structure to pack 
  @param  buf     the byte array in which to pack the refseq data
  @return         a pointer to the next unused byte in memory
  */
uint8_t *
tmap_refseq_shm_pack(tmap_refseq_t *refseq, uint8_t *buf);

/*! 
  @param  buf  the byte array in which to unpack the refseq data
  @return      a pointer to the initialized refseq structure
  */
tmap_refseq_t *
tmap_refseq_shm_unpack(uint8_t *buf);

/*! 
  @param  refseq  pointer to the structure in which the data is stored
  */
void
tmap_refseq_destroy(tmap_refseq_t *refseq);

/*! 
  @param  refseq      pointer to the structure in which the data is stored
  @param  pacpos      the packed FASTA position
  @param  aln_length  the alignment length
  @param  seqid       the zero-based sequence index to be returned
  @param  pos         the one-based position to be returned
  @return             the one-based position, 0 if not found (i.e. overlaps two chromosomes)
  */
inline uint32_t
tmap_refseq_pac2real(const tmap_refseq_t *refseq, uint32_t pacpos, uint32_t aln_length, uint32_t *seqid, uint32_t *pos);

/*! 
  main-like function for 'tmap fasta2pac'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
tmap_refseq_fasta2pac_main(int argc, char *argv[]);
#endif
