#ifndef FMAP_SFF_H_
#define FMAP_SFF_H_

#define FMAP_SFF_MAGIC 0x2E736666
#define FMAP_SFF_VERSION "\0\0\0\1"

#include <stdint.h>
#include "../io/fmap_file.h"
#include "fmap_string.h"

/*! @typedef
  @field  magic           the magic number for this file
  @field  version         the version number
  @field  index_offset    not currently used (value is 0)
  @field  index_length    not currently used (value is 0)
  @field  n_reads         the number of reads in the file
  @field  gheader_length  the number of bytes in the global header including padding
  @field  key_length      the length of the key sequence used with these reads
  @field  flow_length      the number of nucleotide flows used in this experiment
  @field  flowgram_format  the manner in which signal values are encoded (value is 1)
  @field  flow             the string specifying the ith nucleotide flowed 
  @field  key             the string specifying the ith nucleotide of the sequence key
  */
typedef struct {
    uint32_t magic;
    uint32_t version;
    uint64_t index_offset;
    uint64_t index_length;
    uint32_t n_reads;
    uint32_t gheader_length;
    uint16_t key_length;
    uint16_t flow_length;
    uint8_t flowgram_format;
    fmap_string_t *flow;
    fmap_string_t *key;
} fmap_sff_header_t;

/*! @typedef
  @field  rheader_length      the number of bytes in the 
  @field  name_length         the number of characters in the name of the read (not including the null-terminator)
  @field  n_bases             the number of bases in the read
  @field  clip_qual_left      the 1-based coordinate of the first base after the (quality) left clipped region (zero if no clipping has been applied)
  @field  clip_qual_right     the 1-based coordinate of the first base after the (quality) right clipped region (zero if no clipping has been applied)
  @field  clip_adapter_left   the 1-based coordinate of the first base after the (adapter) left clipped region (zero if no clipping has been applied)
  @field  clip_adapter_right  the 1-based coordinate of the first base after the (adapter) right clipped region (zero if no clipping has been applied)
  @field  name                the read name 
  */
typedef struct {
    uint16_t  rheader_length;
    uint16_t  name_length;
    uint32_t  n_bases;
    uint16_t  clip_qual_left;
    uint16_t  clip_qual_right;
    uint16_t  clip_adapter_left;
    uint16_t  clip_adapter_right;
    fmap_string_t *name;
} fmap_sff_read_header_t;

/*! @typedef
  @field  flowgram    the flowgram 
  @field  flow_index  the 1-based flow index for each base called
  @field  bases      the called bases
  @field  quality    the quality score for each base call
  */
typedef struct {
    uint16_t *flowgram;   
    uint8_t *flow_index; 
    fmap_string_t *bases;
    fmap_string_t *quality;
} fmap_sff_read_t;

/*! @typedef
  @field  gheader  pointer to the global header
  @field  rheader  pointer to the read header
  @field  read     pointer to the read
  */
typedef struct {
    fmap_sff_header_t *gheader;
    fmap_sff_read_header_t *rheader;
    fmap_sff_read_t *read;
} fmap_sff_t;

/*! @function
  @abstract
  @param  fp  the file pointer from which to read
  @return     a pointer to the sff header read in
  */
fmap_sff_header_t *
fmap_sff_header_read(fmap_file_t *fp);

/*! @function
  @abstract
  @param  h  a pointer to the sff header to destroy
  */
void
fmap_sff_header_destroy(fmap_sff_header_t *h);

/*! @function
  @abstract
  @param  fp  the file pointer from which to read
  @return     a pointer to the sff read header read in
  */
fmap_sff_read_header_t *
fmap_sff_read_header_read(fmap_file_t *fp);

/*! @function
  @abstract
  @param  h  a pointer to the sff read header to destroy
  */
void
fmap_sff_read_header_destroy(fmap_sff_read_header_t *rh);

/*! @function
  @abstract
  @param  fp  the file pointer from which to read
  @return     a pointer to the sff read to read in
  */
fmap_sff_read_t *
fmap_sff_read_read(fmap_file_t *fp, fmap_sff_header_t *gh, fmap_sff_read_header_t *rh);

/*! @function
  @abstract
  @param  h  a pointer to the sff read to destroy
  */
void
fmap_sff_read_destroy(fmap_sff_read_t *r);

/*! @function
  @abstract
  @return a pointer to the empty sff 
  */
fmap_sff_t *
fmap_sff_init();

/*! @function
  @abstract
  @param  sff  a pointer to the sff to destroy
  */
void
fmap_sff_destroy(fmap_sff_t *sff);

/*! @function
  @abstract
  @return a pointer to the empty sff 
  */
fmap_sff_t *
fmap_sff_init();

/*! @function
  @abstract
  @param  sff  a pointer to the sff to clone
  @return a pointer to the cloned sff
*/
fmap_sff_t *
fmap_sff_clone(fmap_sff_t *sff);

/*! @function
  @abstract
  @param  sff  a pointer to the sff 
*/
void
fmap_sff_reverse_compliment(fmap_sff_t *sff);

/*! @function
  @abstract
  @param  sff  a pointer to the sff 
*/
void
fmap_sff_to_int(fmap_sff_t *sff);

#endif
