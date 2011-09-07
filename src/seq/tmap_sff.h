/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SFF_H
#define TMAP_SFF_H

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
  The SFF Header.
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
  The SFF Read Header.
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
    uint16_t clip_left;  /*< the number of bases hard-clipped on the right of the sequence, not including the key sequence */
    uint16_t clip_right;  /*< the number of bases hard-clipped on the right of the sequence */
} tmap_sff_read_header_t;

/*! 
  The SFF Read.
  */
typedef struct {
    uint16_t *flowgram;  /*!< the flowgram  */
    uint8_t *flow_index;  /*!< the 1-based flow index for each base called */
    tmap_string_t *bases;  /*!< the called bases */
    tmap_string_t *quality;  /*!< the quality score for each base call */
} tmap_sff_read_t;

/*! 
  The SFF Entry.
  */
typedef struct {
    tmap_sff_header_t *gheader;  /*!< pointer to the global header */
    tmap_sff_read_header_t *rheader;  /*!< pointer to the read header */
    tmap_sff_read_t *read;  /*!< pointer to the read */
    uint8_t is_int;  /*!< 1 if the bases are integer values, 0 otherwise */
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
tmap_sff_reverse(tmap_sff_t *sff);

/*! 
  @param  sff  a pointer to the sff 
*/
void
tmap_sff_reverse_compliment(tmap_sff_t *sff);

/*! 
  @param  sff  a pointer to the sff 
*/
void
tmap_sff_compliment(tmap_sff_t *sff);

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
  removes the key sequence and clipped sequence from the read and quality fields
  @param  sff              pointer to the structure to convert
  @param  remove_clipping  1 if we are to remove clipped sequence, 0 otherwise
  @details                 this will not modify the header
  */
inline void
tmap_sff_remove_key_sequence(tmap_sff_t *sff, int32_t remove_clipping);

/*!
  @param  sff        pointer to the structure to convert
  @param  flow_order a pointer to where the flow order should be stored 
  @return            the length of the flow order
 */
int32_t
tmap_sff_get_flow_order_int(tmap_sff_t *sff, uint8_t **flow_order);

/*!
  @param  sff      pointer to the structure to convert
  @param  key_seq  pointer to where the key sequence should be stored 
  @return          the length of the key sequence
 */
int32_t
tmap_sff_get_key_seq_int(tmap_sff_t *sff, uint8_t **key_seq);

/*!
  @param  sff      pointer to the structure to convert
  @param  flowgram  pionter to where the flowgram should be stored
  @param  mem      memory size already allocated to flowgram
  @return          the flowgram length
 */
int32_t
tmap_sff_get_flowgram(tmap_sff_t *sff, uint16_t **flowgram, int32_t mem);

#endif
