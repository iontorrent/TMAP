#ifndef FMAP_SEQ_H_
#define FMAP_SEQ_H_

#include "fmap_fq.h"
#include "fmap_sff.h"

/*! 
  An Abstract Library for DNA Sequence Data
  */

/*! 
  details  the type of DNA sequence data
  */
enum {
    FMAP_SEQ_TYPE_NOTYPE = -1, /*!< unknown type */
    FMAP_SEQ_TYPE_FQ = 0, /*!< FASTA/FASTQ input/output */
    FMAP_SEQ_TYPE_SFF = 1 /*!< SFF input/output */
};

/*! 
  */
typedef struct {
    int8_t type;  /*!< the type associated with this structure */
    union {
        fmap_fq_t *fq;  /*!< the pointer to the fastq structure */
        fmap_sff_t *sff;  /*!< the pointer to the sff structure */
    } data;
} fmap_seq_t;

/*! 
  @param  type  the type associated with this structure
  @return       pointer to the initialized memory 
  */
fmap_seq_t *
fmap_seq_init(int8_t type);

/*! 
  @param  seq  pointer to the structure
  */
void
fmap_seq_destroy(fmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  @return      pointer to the initialized memory 
  */
fmap_seq_t *
fmap_seq_clone(fmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  */
void
fmap_seq_reverse_compliment(fmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  */
void
fmap_seq_to_int(fmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  @return      a pointer to the name string
  */
fmap_string_t *
fmap_seq_get_name(fmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  @return      a pointer to the base sequence string
  */
fmap_string_t *
fmap_seq_get_bases(fmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  @return      a pointer to the quality string
  */
fmap_string_t *
fmap_seq_get_qualities(fmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to convert
  @return      a pointer to the fq structure
  */
fmap_seq_t *
fmap_seq_sff2fq(fmap_seq_t *seq);

#endif
