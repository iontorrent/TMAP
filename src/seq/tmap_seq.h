/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SEQ_H_
#define TMAP_SEQ_H_

#include "tmap_fq.h"
#include "tmap_sff.h"

/*! 
  An Abstract Library for DNA Sequence Data
  */

/*! 
  @details  the type of DNA sequence data
  */
enum {
    TMAP_SEQ_TYPE_NOTYPE = -1, /*!< unknown type */
    TMAP_SEQ_TYPE_FQ = 0, /*!< FASTA/FASTQ input/output */
    TMAP_SEQ_TYPE_SFF = 1 /*!< SFF input/output */
};

/*! 
  */
typedef struct {
    int8_t type;  /*!< the type associated with this structure */
    union {
        tmap_fq_t *fq;  /*!< the pointer to the fastq structure */
        tmap_sff_t *sff;  /*!< the pointer to the sff structure */
    } data;
} tmap_seq_t;

/*! 
  @param  type  the type associated with this structure
  @return       pointer to the initialized memory 
  */
tmap_seq_t *
tmap_seq_init(int8_t type);

/*! 
  @param  seq  pointer to the structure
  */
void
tmap_seq_destroy(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  @return      pointer to the initialized memory 
  */
tmap_seq_t *
tmap_seq_clone(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  */
void
tmap_seq_reverse_compliment(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  */
void
tmap_seq_to_int(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  @return      a pointer to the name string
  */
tmap_string_t *
tmap_seq_get_name(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  @return      a pointer to the base sequence string
  */
inline tmap_string_t *
tmap_seq_get_bases(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to clone
  @return      a pointer to the quality string
  */
inline tmap_string_t *
tmap_seq_get_qualities(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to convert
  @details     this will only remove the key sequence from a SFF 
  structure, and then only the read and quality (not the read header etc.)
  */
void
tmap_seq_remove_key_sequence(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to convert
  @return      a pointer to the fq structure
  */
tmap_seq_t *
tmap_seq_sff2fq(tmap_seq_t *seq);

#endif
