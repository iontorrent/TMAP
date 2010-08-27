#ifndef FMAP_REFSEQ_H_
#define FMAP_REFSEQ_H_

#include <stdint.h>
#include "fmap_io.h"

// version id of the indexes
#define FMAP_REFSEQ_VERSION_ID 0
// seed for our random number generator
#define FMAP_REFSEQ_SEED 13
// refseq file extension
#define FMAP_REFSEQ_ANNO_FILE_EXTENSION ".fmap.refseq.anno"
#define FMAP_REFSEQ_PAC_FILE_EXTENSION ".fmap.refseq.pac"
#define FMAP_REFSEQ_BWT_FILE_EXTENSION ".fmap.refseq.bwt"
#define FMAP_REFSEQ_SA_FILE_EXTENSION ".fmap.refseq.sa"
// the implementation relies on no compression
#define FMAP_REFSEQ_ANNO_COMPRESSION FMAP_FILE_NO_COMPRESSION 
#define FMAP_REFSEQ_PAC_COMPRESSION FMAP_FILE_NO_COMPRESSION 
#define FMAP_REFSEQ_BWT_COMPRESSION FMAP_FILE_NO_COMPRESSION 
#define FMAP_REFSEQ_SA_COMPRESSION FMAP_FILE_NO_COMPRESSION 
// buffer size for reading in from a FASTA file
#define FMAP_REFSEQ_BUFFER_SIZE 0x10000

enum {
    FMAP_REFSEQ_ANNO_FILE = 0,
    FMAP_REFSEQ_PAC_FILE = 1,
    FMAP_REFSEQ_BWT_FILE = 2,
    FMAP_REFSEQ_SA_FILE = 3,
};

/*! @macro
  @param  _len  the number of bases stored 
  @return       the number of bytes allocated
  */
#define fmap_refseq_seq_memory(_len) ((((_len)-1) >> 2) + 1)

/*! @macro
  @param  _i  the 0-based base position 
  @return     the 0-based byte index
  */
#define fmap_refseq_seq_byte_i(_i) ((_i) >> 2)

/*! @macro
  @param  _i  the 0-based base position 
  @return     the number of bits the base is shifted (returns a multiple of two)
  @dicsussion the reference is stored in a 2-bit format
  */
#define fmap_refseq_seq_byte_shift(_i) ((3 - ((_i) & 3)) << 1)

/*! @macro
  @param  _refseq  pointer to the structure holding the reference sequence
  @param  _i       the 0-based base position to store
  @param  _c       the base's 2-bit integer representation
  */
#define fmap_refseq_seq_store_i(_refseq, _i, _c) (_refseq->seq[fmap_refseq_seq_byte_i(_i)] |= _c << fmap_refseq_seq_byte_shift(_i))

/*! @macro
  @param  _refseq  pointer to the structure holding the reference sequence
  @param  _i       the 0-based base position to retrieve
  @return          the base's 2-bit integer representation
  */
#define fmap_refseq_seq_i(_refseq, _i) ((_refseq->seq[fmap_refseq_seq_byte_i(_i)] >> fmap_refseq_seq_byte_shift(_i)) & 3)

/*! @typedef
  @abstract  
  @field  name    the name of the contig
  @field  len     the length of the current contig 
  @field  offset  the offset from the start of the reference
  */
typedef struct {
    char *name;
    uint64_t len;
    uint64_t offset;
} fmap_anno_t;

/*! @typedef
  @abstract  
  @field  version_id    the version id of this file
  @field  seed          the random base generator seed
  @field  seq           the packed nucleotide sequence, with contigs concatenated
  @field  annos         the annotations about the contigs
  @field  len           the total length of the reference sequence
  */
typedef struct {
    uint64_t version_id;
    uint32_t seed;
    uint8_t *seq;
    fmap_anno_t *annos;
    uint32_t num_annos;
    uint64_t len;
} fmap_refseq_t;

/*! @typedef
  @abstract       gets the name of a specific file based on the reference sequence
  @param  prefix   the prefix of the file to be written, usually the fasta file name 
  @param  type    the type of file based on this reference sequence
  @return         a pointer to the file name string
  */
inline char *
fmap_refseq_get_file_name(const char *prefix, int32_t type);

/*! @function
  @abstract
  @field  fn_fasta     the file name of the fasta file
  @field  compression  the type of compression, if any to be used
  @return             the length of the reference sequence
  */
uint64_t
fmap_refseq_fasta2pac(const char *fn_fasta, int32_t compression);

/*! @function
  @abstract
  @field  refseq  pointer to the structure in which to store the data 
  @field  prefix   the prefix of the file to be written, usually the fasta file name 
  */
void
fmap_refseq_write(fmap_refseq_t *refseq, const char *prefix);

/*! @function
  @abstract
  @field  prefix   the prefix of the file to be read, usually the fasta file name 
  @return        a pointer to the initialized memory
  */
fmap_refseq_t *
fmap_refseq_read(const char *prefix);

/*! @function
  @abstract
  @field  refseq  pointer to the structure in which the data is stored
  */
void
fmap_refseq_destroy(fmap_refseq_t *refseq);

/*! @function
  @abstract
  @field  refseq  pointer to the structure in which the data is stored
  @return        the number of bytes used by this data structure
  */
size_t
fmap_refseq_size(fmap_refseq_t *refseq);

/*! @function
  @abstract
  @field  fn_fasta     the file name of the fasta file
  @return             the number of bytes used by the sequence 
  */
int64_t 
fmap_refseq_seq_len(const char *fn_fasta);

/*! @function
  @abstract
  @field  refseq      pointer to the structure in which the data is stored
  @field  pacpos      the packed FASTA position
  @field  aln_length  the alignment length
  @field  seqid       the zero-based sequence index to be returned
  @field  pos         the zero-based position to be returned
  @return        the zero-based position, -1 if not found (overlaps two chromosomes)
  */
inline int32_t
fmap_refseq_pac2real(fmap_refseq_t *refseq, uint32_t pacpos, uint32_t aln_length, uint32_t *seqid, uint32_t *pos);

/*! @function
  @abstract     main-like function for 'fmap fasta2pac'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
fmap_refseq_fasta2pac_main(int argc, char *argv[]);
#endif
