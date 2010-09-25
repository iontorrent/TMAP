#ifndef FMAP_SEQ_H_
#define FMAP_SEQ_H_

#include "fmap_fq.h"
#include "fmap_sff.h"

/*! @header
  @abstract  An Abstract Library for DNA Sequence Data
  */

/*! @enum  Sequence Type 
  @constant  FMAP_SEQ_TYPE_NOTYPE    unknown type
  @constant  FMAP_SEQ_TYPE_FQ        FASTA/FASTQ input/output
  @constant  FMAP_SEQ_TYPE_SFF       SFF input/output
  @discussion  the type of DNA sequence data
  */
enum {
    FMAP_SEQ_TYPE_NOTYPE = -1,
    FMAP_SEQ_TYPE_FQ = 0,
    FMAP_SEQ_TYPE_SFF = 1
};

/*! @typedef 
  @field  type  the type associated with this structure
  @field  data  pointer to the particular read data structure
  */
typedef struct {
    int8_t type;
    union {
        fmap_fq_t *fq;
        fmap_sff_t *sff;
    } data;
} fmap_seq_t;

/*! @function
  @abstract
  @param  type  the type associated with this structure
  @return       pointer to the initialized memory 
  */
fmap_seq_t *
fmap_seq_init(int8_t type);

/*! @function
  @abstract
  @param  seq  pointer to the structure
  */
void
fmap_seq_destroy(fmap_seq_t *seq);

/*! @function
  @abstract
  @param  seq  pointer to the structure to clone
  @return      pointer to the initialized memory 
  */
fmap_seq_t *
fmap_seq_clone(fmap_seq_t *seq);

/*! @function
  @abstract
  @param  seq  pointer to the structure to clone
  */
void
fmap_seq_reverse_compliment(fmap_seq_t *seq);

/*! @function
  @abstract
  @param  seq  pointer to the structure to clone
  */
void
fmap_seq_to_int(fmap_seq_t *seq);

/*! @function
  @abstract
  @param  seq  pointer to the structure to clone
  @return      a pointer to the name string
  */
fmap_string_t *
fmap_seq_get_name(fmap_seq_t *seq);

/*! @function
  @abstract
  @param  seq  pointer to the structure to clone
  @return      a pointer to the base sequence string
  */
fmap_string_t *
fmap_seq_get_bases(fmap_seq_t *seq);

/*! @function
  @abstract
  @param  seq  pointer to the structure to clone
  @return      a pointer to the quality string
  */
fmap_string_t *
fmap_seq_get_qualities(fmap_seq_t *seq);

/*! @function
  @abstract
  @param  seq  pointer to the structure to convert
  @return      a pointer to the fq structure
  */
fmap_seq_t *
fmap_seq_sff2fq(fmap_seq_t *seq);

#endif