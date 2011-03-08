/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_DEFINITIONS_H
#define TMAP_DEFINITIONS_H

#include <stdint.h>
#include <config.h>

/*! 
  Generic Functions
  */

/*! d TMAP_VERSION_ID
  the magic id for tmap
  */
#define TMAP_VERSION_ID ('f' + 'm' + 'a' + 'p')

/* 
 * File extensions
 */
/*! d TMAP_ANNO_FILE_EXTENSION
  the file extension for the reference sequence annotations
  */
#define TMAP_ANNO_FILE_EXTENSION ".tmap.anno"
/*! d TMAP_PAC_FILE_EXTENSION
  the file extension for the packed forward reference sequence
  */
#define TMAP_PAC_FILE_EXTENSION ".tmap.pac"
/*! d TMAP_REV_PAC_FILE_EXTENSION
  the file extension for the packed reverse reference sequence
  */
#define TMAP_REV_PAC_FILE_EXTENSION ".tmap.rpac"
/*! d TMAP_BWT_FILE_EXTENSION
  the file extension for the forward BWT structure
  */
#define TMAP_BWT_FILE_EXTENSION ".tmap.bwt"
/*! d TMAP_REV_BWT_FILE_EXTENSION
  the file extension for the reverse BWT structure
  */
#define TMAP_REV_BWT_FILE_EXTENSION ".tmap.rbwt"
/*! d TMAP_SA_FILE_EXTENSION
  the file extension for the forward SA structure
  */
#define TMAP_SA_FILE_EXTENSION ".tmap.sa"
/*! d TMAP_REV_SA_FILE_EXTENSION
  the file extension for the reverse SA structure
  */
#define TMAP_REV_SA_FILE_EXTENSION ".tmap.rsa"

// The default compression types for each file
// Note: the implementation relies on no compression
#define TMAP_ANNO_COMPRESSION TMAP_FILE_NO_COMPRESSION 
#define TMAP_PAC_COMPRESSION TMAP_FILE_NO_COMPRESSION 
#define TMAP_REV_PAC_COMPRESSION TMAP_FILE_NO_COMPRESSION 
#define TMAP_BWT_COMPRESSION TMAP_FILE_NO_COMPRESSION 
#define TMAP_REV_BWT_COMPRESSION TMAP_FILE_NO_COMPRESSION 
#define TMAP_SA_COMPRESSION TMAP_FILE_NO_COMPRESSION
#define TMAP_REV_SA_COMPRESSION TMAP_FILE_NO_COMPRESSION

/*
   CIGAR operations, from samtools.
   */
#ifndef BAM_BAM_H
#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#endif

/*! 
  for each type of file, the integer id associated with this file
  @details  can be used with 'tmap_get_file_name' 
  */
enum {
    TMAP_ANNO_FILE     = 0, /*!< the reference sequence annotation file */
    TMAP_PAC_FILE      = 1, /*!< the packed forward reference sequence file */
    TMAP_REV_PAC_FILE  = 2, /*!< the packed reverse reference sequence file */
    TMAP_BWT_FILE      = 3, /*!< the packed forward BWT file */
    TMAP_REV_BWT_FILE  = 4, /*!< the packed reverse BWT file */
    TMAP_SA_FILE       = 5, /*!< the packed forward SA file */
    TMAP_REV_SA_FILE   = 6 /*!< the packed reverse SA file */
};

/*! 
*/
enum {
    TMAP_READS_FORMAT_UNKNOWN  = -1, /*!< the reads format is unrecognized */
    TMAP_READS_FORMAT_FASTA    = 0, /*!< the reads are in FASTA format */
    TMAP_READS_FORMAT_FASTQ    = 1, /*!< the reads are in FASTQ format */
    TMAP_READS_FORMAT_SFF      = 2, /*!< the reads are in SFF format */
#ifdef HAVE_SAMTOOLS
    TMAP_READS_FORMAT_SAM      = 3, /*!< the reads are in SAM format */
    TMAP_READS_FORMAT_BAM      = 4 /*!< the reads are in BAM format */
#endif
};

/*! 
  @param  algo_id  the algorithm identifier
  @return          algorithm name
  */
char *
tmap_algo_id_to_name(uint16_t algo_id);

/*! nt_char_to_int
  @details  converts a DNA base in ASCII format to its 2-bit format [0-4]. 
  */
extern uint8_t tmap_nt_char_to_int[256];

/*! nt_char_to_rc_char
  @details  converts a DNA base in ASCII format to reverse compliment in ASCII format.
  */
extern uint8_t tmap_nt_char_to_rc_char[256];

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

#ifndef tmap_roundup32
/*! 
  rounds up to the nearest power of two integer
  @param  x  the integer to round up
  @return    the smallest integer greater than x that is a power of two 
  */
#define tmap_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/*!
  @param  reads_format  the reads format
  @return               the sequence format (for tmap_seq_t)
  */
int32_t 
tmap_reads_format_to_seq_type(int32_t reads_format);

/*! 
  @param  v  the value to take the log 2
  @return    log of the value, base two
  */
inline uint32_t 
tmap_log2(uint32_t v);

/*! 
  gets the name of a specific file based on the reference sequence
  @param  prefix   the prefix of the file to be written, usually the fasta file name 
  @param  type    the type of file based on this reference sequence
  @return         a pointer to the file name string
  */
inline char *
tmap_get_file_name(const char *prefix, int32_t type);

/*! 
  @param  optarg  the string of the file format
  @return         the format type
  */
int 
tmap_get_reads_file_format_int(char *optarg);

/*! 
  checks the extension of the file to recognize its format     
  @param  fn            the file name 
  @param  reads_format  pointer to the reads format, if any (unknown|fastq|fq|fasta|fa|sff)
  @param  compr_type    pointer the type of compression used, if any (none|gz|bz2)
  @details              if the reads_format is unknown, it will be populated; similarly for compr_type.
  */
void
tmap_get_reads_file_format_from_fn_int(char *fn, int32_t *reads_format, int32_t *compr_type);

/*! 
  @param  format  the interger file format specifier
  @return         the format type (string)
  */
char *
tmap_get_reads_file_format_string(int format);

/*!
  reverses a given string
  @param  seq  the string to reverse
  @param  len  the length of the string
  */
inline void
tmap_reverse(char *seq, int32_t len);

/*!
  reverse compliments a given string
  @param  seq  the character DNA sequence
  @param  len  the length of the DNA sequence
  */
inline void
tmap_reverse_compliment(char *seq, int32_t len); 

/*!
  compliments a given string
  @param  seq  the character DNA sequence
  @param  len  the length of the DNA sequence
  */
inline void
tmap_compliment(char *seq, int32_t len); 

#endif
