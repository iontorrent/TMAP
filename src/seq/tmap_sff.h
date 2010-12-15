/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SFF_H_
#define TMAP_SFF_H_

#define TMAP_SFF_MAGIC 0x2E736666
#define TMAP_SFF_VERSION 1

// uncomment this to allow for some SFF debuggin
//#define TMAP_SFF_DEBUG 1

#include <stdint.h>
#include "../util/tmap_string.h"
#include "../io/tmap_file.h"

/*! 
  A Library for SFF data
  */

/*! 
  */
typedef struct {
    uint32_t magic;  /*!< the magic number for this file */
    uint32_t version;  /*!< the version number */
    uint64_t index_offset;  /*!< not currently used (value is 0) */
    uint32_t index_length;  /*!< not currently used (value is 0) */
    uint32_t n_reads;  /*!< the number of reads in the file */
    uint32_t gheader_length;  /*!< the number of bytes in the global header including padding */
    uint16_t key_length;  /*!< the length of the key sequence used with these reads */
    uint16_t flow_length;  /*!< the number of nucleotide flows used in this experiment */
    uint8_t flowgram_format;  /*!< the manner in which signal values are encoded (value is 1) */
    tmap_string_t *flow;  /*!< the string specifying the ith nucleotide flowed  */
    tmap_string_t *key;  /*!< the string specifying the ith nucleotide of the sequence key */
} tmap_sff_header_t;

/*! 
  */
typedef struct {
    uint16_t rheader_length;  /*!< the number of bytes in the  */
    uint16_t name_length;  /*!< the number of characters in the name of the read (not including the null-terminator) */
    uint32_t n_bases;  /*!< the number of bases in the read */
    uint16_t clip_qual_left;  /*!< the 1-based coordinate of the first base after the (quality) left clipped region (zero if no clipping has been applied) */
    uint16_t clip_qual_right;  /*!< the 1-based coordinate of the first base after the (quality) right clipped region (zero if no clipping has been applied) */
    uint16_t clip_adapter_left;  /*!< the 1-based coordinate of the first base after the (adapter) left clipped region (zero if no clipping has been applied) */
    uint16_t clip_adapter_right;  /*!< the 1-based coordinate of the first base after the (adapter) right clipped region (zero if no clipping has been applied) */
    tmap_string_t *name;  /*!< the read name  */
} tmap_sff_read_header_t;

/*! 
  */
typedef struct {
    uint16_t *flowgram;  /*!< the flowgram  */
    uint8_t *flow_index;  /*!< the 1-based flow index for each base called */
    tmap_string_t *bases;  /*!< the called bases */
    tmap_string_t *quality;  /*!< the quality score for each base call */
} tmap_sff_read_t;

/*! 
  */
typedef struct {
    tmap_sff_header_t *gheader;  /*!< pointer to the global header */
    tmap_sff_read_header_t *rheader;  /*!< pointer to the read header */
    tmap_sff_read_t *read;  /*!< pointer to the read */
    int32_t is_int;  /*!< 1 if the bases are integer values, 0 otherwise */
} tmap_sff_t;

/*! 
  @param  fp  the file pointer from which to read
  @return     a pointer to the sff header read in
  */
tmap_sff_header_t *
tmap_sff_header_read(tmap_file_t *fp);

/*! 
  @param  h  a pointer to the sff header to destroy
  */
void
tmap_sff_header_destroy(tmap_sff_header_t *h);

/*! 
  @param  fp  the file pointer from which to read
  @return     a pointer to the sff read header read in
  */
tmap_sff_read_header_t *
tmap_sff_read_header_read(tmap_file_t *fp);

/*! 
  @param  rh  a pointer to the sff read header to destroy
  */
void
tmap_sff_read_header_destroy(tmap_sff_read_header_t *rh);

/*! 
  @param  fp  the file pointer from which to read
  @param  gh  the sff global header
  @param  rh  the sff read header
  @return     a pointer to the sff read to read in
  */
tmap_sff_read_t *
tmap_sff_read_read(tmap_file_t *fp, tmap_sff_header_t *gh, tmap_sff_read_header_t *rh);

/*! 
  @param  r  a pointer to the sff read to destroy
  */
void
tmap_sff_read_destroy(tmap_sff_read_t *r);

/*! 
  @return a pointer to the empty sff 
  */
tmap_sff_t *
tmap_sff_init();

/*! 
  @param  sff  a pointer to the sff to destroy
  */
void
tmap_sff_destroy(tmap_sff_t *sff);

/*! 
  @param  sff  a pointer to the sff to clone
  @return a pointer to the cloned sff
*/
tmap_sff_t *
tmap_sff_clone(tmap_sff_t *sff);

/*! 
  @param  sff  a pointer to the sff 
*/
void
tmap_sff_reverse_compliment(tmap_sff_t *sff);

/*! 
  @param  sff  a pointer to the sff 
*/
void
tmap_sff_to_int(tmap_sff_t *sff);

/*! 
  @param  sff  a pointer to the sff 
*/
void
tmap_sff_to_char(tmap_sff_t *sff);
/*!
  gets the read's bases
  @param  sff  a pointer to a sequence structure
  @details     this will include the key sequence qualities
 */
inline tmap_string_t *
tmap_sff_get_bases(tmap_sff_t *sff);

/*!
  gets the read's qualities
  @param  sff  a pointer to a sequence structure
  @details     this will include the key sequence qualities
 */
inline tmap_string_t *
tmap_sff_get_qualities(tmap_sff_t *sff);

/*! 
  removes the key sequence from the read and quality fields
  @param  sff  pointer to the structure to convert
  @details     this will only remove the key sequence from the SFF
  structure, and then only the read and quality (not the read header etc.)
  */
inline void
tmap_sff_remove_key_sequence(tmap_sff_t *sff);

#endif
