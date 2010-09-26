#ifndef FMAP_DEFINITIONS_H_
#define FMAP_DEFINITIONS_H_

#include <stdint.h>

/*! 
  details  Generic Functions
  */

/*! d FMAP_VERSION_ID
  the magic id for fmap
  */
#define FMAP_VERSION_ID ('f' + 'm' + 'a' + 'p')

/* 
 * File extensions
 */
/*! d FMAP_ANNO_FILE_EXTENSION
  the file extension for the reference sequence annotations
  */
#define FMAP_ANNO_FILE_EXTENSION ".fmap.anno"
/*! d FMAP_PAC_FILE_EXTENSION
  the file extension for the packed forward reference sequence
  */
#define FMAP_PAC_FILE_EXTENSION ".fmap.pac"
/*! d FMAP_REV_PAC_FILE_EXTENSION
  the file extension for the packed reverse reference sequence
  */
#define FMAP_REV_PAC_FILE_EXTENSION ".fmap.rpac"
/*! d FMAP_BWT_FILE_EXTENSION
  the file extension for the forward BWT structure
  */
#define FMAP_BWT_FILE_EXTENSION ".fmap.bwt"
/*! d FMAP_REV_BWT_FILE_EXTENSION
  the file extension for the reverse BWT structure
  */
#define FMAP_REV_BWT_FILE_EXTENSION ".fmap.rbwt"
/*! d FMAP_SA_FILE_EXTENSION
  the file extension for the forward SA structure
  */
#define FMAP_SA_FILE_EXTENSION ".fmap.sa"
/*! d FMAP_REV_SA_FILE_EXTENSION
  the file extension for the reverse SA structure
  */
#define FMAP_REV_SA_FILE_EXTENSION ".fmap.rsa"

// The default compression types for each file
// Note: the implementation relies on no compression
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
#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6

/*! 
  for each type of file, the integer id associated with this file
  details  can be used with 'fmap_get_file_name' 
  */
enum {
    FMAP_ANNO_FILE     = 0, /*!< the reference sequence annotation file */
    FMAP_PAC_FILE      = 1, /*!< the packed forward reference sequence file */
    FMAP_REV_PAC_FILE  = 2, /*!< the packed reverse reference sequence file */
    FMAP_BWT_FILE      = 3, /*!< the packed forward BWT file */
    FMAP_REV_BWT_FILE  = 4, /*!< the packed reverse BWT file */
    FMAP_SA_FILE       = 5, /*!< the packed forward SA file */
    FMAP_REV_SA_FILE   = 6 /*!< the packed reverse SA file */
};

/*! 
  */
enum {
    FMAP_READS_FORMAT_UNKNOWN  = -1, /*!< the reads format is unrecognized */
    FMAP_READS_FORMAT_FASTA    = 0, /*!< the reads are in FASTA format */
    FMAP_READS_FORMAT_FASTQ    = 1, /*!< the reads are in FASTQ format */
    FMAP_READS_FORMAT_SFF      = 2 /*!< the reads are in SFF format */
};

/*! @var  nt_char_to_int
  details  converts a DNA base in ASCII format to its 2-bit format [0-4]. 
  */
extern uint8_t nt_char_to_int[256];

/*! @var  nt_char_to_rc_char
  details  converts a DNA base in ASCII format to reverse compliment in ASCII format.
  */
extern uint8_t nt_char_to_rc_char[256];

/*! 
  @param  c  the quality value in ASCII format
  @return    the quality value in integer format
  */
#define CHAR2QUAL(c) ((uint8_t)c-33)

/*! 
  @param  q  the quality value in integer format
  @return    the quality value in ASCII format
  */
#define QUAL2CHAR(q) (char)(((((unsigned char)q)<=93)?q:93)+33)

#ifndef htonll
/*! 
  converts a 64-bit value to network order
  @param  x  the 64-bit value to convert
  @return    the converted 64-bit value
  */
#define htonll(x) ((((uint64_t)htonl(x)) << 32) + htonl(x >> 32))
#endif

#ifndef ntohll
/*! 
  converts a 64-bit value to host order
  @param  x  the 64-bit value to convert
  @return    the converted 64-bit value
  */
#define ntohll(x) ((((uint64_t)ntohl(x)) << 32) + ntohl(x >> 32))
#endif

#ifndef fmap_roundup32
/*! 
  rounds up to the nearest power of two integer
  @param  x  the integer to round up
  @return    the smallest integer greater than x that is a power of two 
  */
#define fmap_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/*! 
  @param  v  the value to take the log 2
  @return    log of the value, base two
  */
inline uint32_t 
fmap_log2(uint32_t v);

/*! 
       gets the name of a specific file based on the reference sequence
  @param  prefix   the prefix of the file to be written, usually the fasta file name 
  @param  type    the type of file based on this reference sequence
  @return         a pointer to the file name string
  */
inline char *
fmap_get_file_name(const char *prefix, int32_t type);

/*! 
       
  @param  optarg  the string of the file format
  @return         the format type
  */
int 
fmap_get_reads_file_format_int(char *optarg);

/*! 
             checks the extension of the file to recognize its format     
  @param  fn            the file name 
  @param  reads_format  pointer to the reads format, if any (unknown|fastq|fq|fasta|fa|sff)
  @param  compr_type    pointer the type of compression used, if any (none|gz|bz2)
  details           if the reads_format is unknown, it will be populated; similarly for compr_type.
  */
void
fmap_get_reads_file_format_from_fn_int(char *fn, int32_t *reads_format, int32_t *compr_type);

/*! 
       
  @param  format  the interger file format specifier
  @return         the format type (string)
  */
char *
fmap_get_reads_file_format_string(int format);

#endif
