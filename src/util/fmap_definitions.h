#ifndef FMAP_DEFINITIONS_H_
#define FMAP_DEFINITIONS_H_

#include <stdint.h>

// TODO: incorporate this into all static files
#define FMAP_VERSION_ID ('f' + 'm' + 'a' + 'p')

// refseq file extension
#define FMAP_ANNO_FILE_EXTENSION ".fmap.anno"
#define FMAP_PAC_FILE_EXTENSION ".fmap.pac"
#define FMAP_REV_PAC_FILE_EXTENSION ".fmap.rpac"
#define FMAP_BWT_FILE_EXTENSION ".fmap.bwt"
#define FMAP_REV_BWT_FILE_EXTENSION ".fmap.rbwt"
#define FMAP_SA_FILE_EXTENSION ".fmap.sa"
#define FMAP_REV_SA_FILE_EXTENSION ".fmap.rsa"
// the implementation relies on no compression
#define FMAP_ANNO_COMPRESSION FMAP_FILE_NO_COMPRESSION 
#define FMAP_PAC_COMPRESSION FMAP_FILE_NO_COMPRESSION 
#define FMAP_REV_PAC_COMPRESSION FMAP_FILE_NO_COMPRESSION 
#define FMAP_BWT_COMPRESSION FMAP_FILE_NO_COMPRESSION 
#define FMAP_REV_BWT_COMPRESSION FMAP_FILE_NO_COMPRESSION 
#define FMAP_SA_COMPRESSION FMAP_FILE_NO_COMPRESSION
#define FMAP_REV_SA_COMPRESSION FMAP_FILE_NO_COMPRESSION

/*
   CIGAR operations.
   */
/*! @abstract CIGAR: match */
#define BAM_CMATCH      0
/*! @abstract CIGAR: insertion to the reference */
#define BAM_CINS        1
/*! @abstract CIGAR: deletion from the reference */
#define BAM_CDEL        2
/*! @abstract CIGAR: skip on the reference (e.g. spliced alignment) */
#define BAM_CREF_SKIP   3
/*! @abstract CIGAR: clip on the read with clipped sequence present in qseq */
#define BAM_CSOFT_CLIP  4
/*! @abstract CIGAR: clip on the read with clipped sequence trimmed off */
#define BAM_CHARD_CLIP  5
/*! @abstract CIGAR: padding */
#define BAM_CPAD        6

// TODO: document
enum {
    FMAP_ANNO_FILE     = 0,
    FMAP_PAC_FILE      = 1,
    FMAP_REV_PAC_FILE  = 2,
    FMAP_BWT_FILE      = 3,
    FMAP_REV_BWT_FILE  = 4,
    FMAP_SA_FILE       = 5,
    FMAP_REV_SA_FILE   = 6
};

// TODO: document
enum {
    FMAP_READS_FORMAT_UNKNOWN  = -1,
    FMAP_READS_FORMAT_FASTA    = 0,
    FMAP_READS_FORMAT_FASTQ    = 1,
    FMAP_READS_FORMAT_SFF      = 2
};

//TODO: document
extern uint8_t nt_char_to_int[256];
//TODO: document
extern uint8_t nt_char_to_rc_char[256];
//TODO: document
#define CHAR2QUAL(c) ((uint8_t)c-33)
//TODO: document
#define QUAL2CHAR(q) (char)(((q<=93)?q:93)+33)

//TODO: document
#ifndef htonll
#define htonll(x) ((((uint64_t)htonl(x)) << 32) + htonl(x >> 32))
#endif
//TODO: document
#ifndef ntohll
#define ntohll(x) ((((uint64_t)ntohl(x)) << 32) + ntohl(x >> 32))
#endif

/*! @function
  @abstract       gets the name of a specific file based on the reference sequence
  @param  prefix   the prefix of the file to be written, usually the fasta file name 
  @param  type    the type of file based on this reference sequence
  @return         a pointer to the file name string
  */
inline char *
fmap_get_file_name(const char *prefix, int32_t type);

/*! @function
  @abstract       
  @param  optarg  the string of the file format
  @return         the format type
  */
int 
fmap_get_reads_file_format_int(char *optarg);

/*! @function
  @abstract   checks the extension of the file to recognize its format     
  @param  fn  the file name 
  @return     the format type
  */
int 
fmap_get_reads_file_format_from_fn_int(char *fn);


/*! @function
  @abstract       
  @param  format  the interger file format specifier
  @return         the format type (string)
  */
char *
fmap_get_reads_file_format_string(int format);

#endif
