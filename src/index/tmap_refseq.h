/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_REFSEQ_H
#define TMAP_REFSEQ_H

#include <stdint.h>
#include "../util/tmap_string.h"
#include "../util/tmap_definitions.h"
#include "../io/tmap_file.h"

/*! 
  DNA Reference Sequence Library
  */

// buffer size for reading in from a FASTA file
#define TMAP_REFSEQ_BUFFER_SIZE 0x10000

// the number of bases on a given line for "tmap pac2refseq"
#define TMAP_REFSEQ_FASTA_LINE_LENGTH 72

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
    uint32_t *amb_positions_start;  /*!< start positions of ambiguous bases (one-based) */
    uint32_t *amb_positions_end;  /*!< end positions of ambiguous bases (one-based) */
    uint8_t *amb_bases;  /*!< the ambiguous bases (IUPAC code) */
    uint32_t num_amb;  /*!< the number of ambiguous bases */
} tmap_anno_t;

/*! 
  */
typedef struct {
    uint64_t version_id;  /*!< the version id of this file */
    tmap_string_t *package_version;  /*!< the package version */
    uint8_t *seq;  /*!< the packed nucleotide sequence, with contigs concatenated */
    tmap_anno_t *annos;  /*!< the annotations about the contigs */
    int32_t num_annos;  /*!< the number of contigs (and annotations) */
    uint64_t len;  /*!< the total length of the reference sequence */
    uint32_t is_shm;  /*!< 1 if loaded from shared memory, 0 otherwise */
} tmap_refseq_t;

/*!
  returns the index version format given a package version
  @param  v  the package version string
  @return    the index version format string
  */
const char * 
tmap_refseq_get_version_format(const char *v);

/*! 
  @param  fn_fasta     the file name of the fasta file
  @param  compression  the type of compression, if any to be used
  @param  fwd_only     1 if to pack the forward sequence only, 0 otherwise
  @return              the length of the reference sequence
  */
uint64_t
tmap_refseq_fasta2pac(const char *fn_fasta, int32_t compression, int32_t fwd_only);

/*! 
  @param  fn_fasta     the file name of the fasta file
  */
void
tmap_refseq_pac2revpac(const char *fn_fasta);

/*! 
  @param  refseq    pointer to the structure in which to store the data 
  @param  fn_fasta  the fn_fasta of the file to be written, usually the fasta file name 
  */
void
tmap_refseq_write(tmap_refseq_t *refseq, const char *fn_fasta);

/*! 
  @param  fn_fasta  the fn_fasta of the file to be read, usually the fasta file name 
  @return           a pointer to the initialized memory
  */
tmap_refseq_t *
tmap_refseq_read(const char *fn_fasta);

/*! 
  @param  refseq  the refseq structure 
  @return         the number of bytes required for this bwt in shared memory
  */
size_t
tmap_refseq_shm_num_bytes(tmap_refseq_t *refseq);

/*! 
  @param  fn_fasta  the fn_fasta of the file to be read, usually the fasta file name 
  @return           the number of bytes required for this bwt in shared memory
  */
size_t
tmap_refseq_shm_read_num_bytes(const char *fn_fasta);

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
  @param  pacpos      the packed FASTA position (one-based)
  @param  aln_length  the alignment length
  @param  seqid       the zero-based sequence index to be returned
  @param  pos         the one-based position to be returned
  @param  strand      the strand (0 - forward, 1 - reverse)
  @return             the one-based position, 0 if not found (i.e. overlaps two chromosomes)
  */
inline tmap_bwt_int_t
tmap_refseq_pac2real(const tmap_refseq_t *refseq, tmap_bwt_int_t pacpos, uint32_t aln_length, uint32_t *seqid, uint32_t *pos, uint8_t *strand);

/*! 
  Retrieves a subsequence of the reference in 2-bit format
  @param  refseq  pointer to the structure in which the data is stored
  @param  pacpos  the packed FASTA position (one-based)
  @param  length  the subsequence length retrieve
  @param  target  the target in which to store the data (must be allocated with enough memory)
  @return         the length retrieved
  */
inline int32_t
tmap_refseq_subseq(const tmap_refseq_t *refseq, tmap_bwt_int_t pacpos, uint32_t length, uint8_t *target);

/*! 
  Retrieves a subsequence of the reference in 2-bit format
  @param  refseq  pointer to the structure in which the data is stored
  @param  seqid   the sequence id (one-based)
  @param  start   the start position (one-based)
  @param  end     the end position (one-based)
  @param  target  pre-allocated memory for the target
  @param  to_n    change all ambiguous bases to N, otherwise they will be returned as the correct code
  @param  conv    the number of bases converted to ambiguity bases
  @return         the target sequence if successful, NULL otherwise
  */
inline uint8_t*
tmap_refseq_subseq2(const tmap_refseq_t *refseq, uint32_t seqid, uint32_t start, uint32_t end, uint8_t *target, int32_t to_n, int32_t *conv);

/*! 
  Checks if the given reference range has ambiguous bases
  @param  refseq  pointer to the structure in which the data is stored
  @param  seqid   the sequence index (one-based)
  @param  start   the start position (one-based and inclusive)
  @param  end     the end position (one-based and inclusive)
  @return         zero if none were found, otherwise the one-based index into "amb_bases" array
  */
inline int32_t
tmap_refseq_amb_bases(const tmap_refseq_t *refseq, uint32_t seqid, uint32_t start, uint32_t end);

/*! 
  main-like function for 'tmap fasta2pac'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
tmap_refseq_fasta2pac_main(int argc, char *argv[]);

/*! 
  main-like function for 'tmap refinfo'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
tmap_refseq_refinfo_main(int argc, char *argv[]);

/*! 
  main-like function for 'tmap pac2fasta'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
tmap_refseq_pac2fasta_main(int argc, char *argv[]);
#endif
