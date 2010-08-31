#ifndef FMAP_DEFINTIONS_H_
#define FMAP_DEFINTIONS_H_

#include <stdint.h>

// TODO: incorporate this into all static files
#define FMAP_VERSION_ID 0

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

enum {
    FMAP_ANNO_FILE     = 0,
    FMAP_PAC_FILE      = 1,
    FMAP_REV_PAC_FILE  = 2,
    FMAP_BWT_FILE      = 3,
    FMAP_REV_BWT_FILE  = 4,
    FMAP_SA_FILE       = 5,
    FMAP_REV_SA_FILE   = 6
};

extern uint8_t nt_char_to_int[256];

/*! @typedef
  @abstract       gets the name of a specific file based on the reference sequence
  @param  prefix   the prefix of the file to be written, usually the fasta file name 
  @param  type    the type of file based on this reference sequence
  @return         a pointer to the file name string
  */
inline char *
fmap_get_file_name(const char *prefix, int32_t type);

#endif
