#ifndef FMAP_SFF_H_
#define FMAP_SFF_H_

#define FMAP_SFF_MAGIC 0x2E736666
#define FMAP_SFF_VERSION 1

// uncomment this to allow for some SFF debuggin
//#define FMAP_SFF_DEBUG 1

#include <stdint.h>
#include "../util/fmap_string.h"
#include "../io/fmap_file.h"

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
    fmap_string_t *flow;  /*!< the string specifying the ith nucleotide flowed  */
    fmap_string_t *key;  /*!< the string specifying the ith nucleotide of the sequence key */
} fmap_sff_header_t;

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
    fmap_string_t *name;  /*!< the read name  */
} fmap_sff_read_header_t;

/*! 
  */
typedef struct {
    uint16_t *flowgram;  /*!< the flowgram  */
    uint8_t *flow_index;  /*!< the 1-based flow index for each base called */
    fmap_string_t *bases;  /*!< the called bases */
    fmap_string_t *quality;  /*!< the quality score for each base call */
} fmap_sff_read_t;

/*! 
  */
typedef struct {
    fmap_sff_header_t *gheader;  /*!< pointer to the global header */
    fmap_sff_read_header_t *rheader;  /*!< pointer to the read header */
    fmap_sff_read_t *read;  /*!< pointer to the read */
    int32_t is_int;  /*!< 1 if the bases are integer values, 0 otherwise */
} fmap_sff_t;

/*! 
  @param  fp  the file pointer from which to read
  @return     a pointer to the sff header read in
  */
fmap_sff_header_t *
fmap_sff_header_read(fmap_file_t *fp);

/*! 
  @param  h  a pointer to the sff header to destroy
  */
void
fmap_sff_header_destroy(fmap_sff_header_t *h);

/*! 
  @param  fp  the file pointer from which to read
  @return     a pointer to the sff read header read in
  */
fmap_sff_read_header_t *
fmap_sff_read_header_read(fmap_file_t *fp);

/*! 
  @param  rh  a pointer to the sff read header to destroy
  */
void
fmap_sff_read_header_destroy(fmap_sff_read_header_t *rh);

/*! 
  @param  fp  the file pointer from which to read
  @param  gh  the sff global header
  @param  rh  the sff read header
  @return     a pointer to the sff read to read in
  */
fmap_sff_read_t *
fmap_sff_read_read(fmap_file_t *fp, fmap_sff_header_t *gh, fmap_sff_read_header_t *rh);

/*! 
  @param  r  a pointer to the sff read to destroy
  */
void
fmap_sff_read_destroy(fmap_sff_read_t *r);

/*! 
  @return a pointer to the empty sff 
  */
fmap_sff_t *
fmap_sff_init();

/*! 
  @param  sff  a pointer to the sff to destroy
  */
void
fmap_sff_destroy(fmap_sff_t *sff);

/*! 
  @param  sff  a pointer to the sff to clone
  @return a pointer to the cloned sff
*/
fmap_sff_t *
fmap_sff_clone(fmap_sff_t *sff);

/*! 
  @param  sff  a pointer to the sff 
*/
void
fmap_sff_reverse_compliment(fmap_sff_t *sff);

/*! 
  @param  sff  a pointer to the sff 
*/
void
fmap_sff_to_int(fmap_sff_t *sff);

/*! 
  @param  sff  a pointer to the sff 
*/
void
fmap_sff_to_char(fmap_sff_t *sff);
/*!
  gets the read's bases
  @param  sff  a pointer to a sequence structure
  @details     this will include the key sequence qualities
 */
inline fmap_string_t *
fmap_sff_get_bases(fmap_sff_t *sff);

/*!
  gets the read's qualities
  @param  sff  a pointer to a sequence structure
  @details     this will include the key sequence qualities
 */
inline fmap_string_t *
fmap_sff_get_qualities(fmap_sff_t *sff);

/*! 
  @param  sff  pointer to the structure to convert
  @details     this will only remove the key sequence from the SFF
  structure, and then only the read and quality (not the read header etc.)
  */
void
fmap_sff_remove_key_sffuence(fmap_sff_t *sff);

#endif
