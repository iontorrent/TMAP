/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SEQ_H
#define TMAP_SEQ_H

#include <config.h>
#include "tmap_fq.h"
#include "tmap_sff.h"
#ifdef HAVE_SAMTOOLS
#include "tmap_sam.h"
#endif

/*! 
  An Abstract Library for DNA Sequence Data
  */

// TODO: we could turn this into a macro library

/*! 
  the type of DNA sequence data
  */
enum {
    TMAP_SEQ_TYPE_NOTYPE = -1, /*!< unknown type */
    TMAP_SEQ_TYPE_FQ = 0, /*!< FASTA/FASTQ input/output */
    TMAP_SEQ_TYPE_SFF = 1, /*!< SFF input/output */
#ifdef HAVE_SAMTOOLS
    TMAP_SEQ_TYPE_SAM = 2, /*!< SAM input/output */
    TMAP_SEQ_TYPE_BAM = 3 /*!< BAM input/output */
#endif
};

/*! 
  */
typedef struct {
    int8_t type;  /*!< the type associated with this structure */
    union {
        tmap_fq_t *fq;  /*!< the pointer to the fastq structure */
        tmap_sff_t *sff;  /*!< the pointer to the sff structure */
#ifdef HAVE_SAMTOOLS
        tmap_sam_t *sam;  /*!< the pointer to the SAM/BAM structure */
#endif
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
  @param  seq  pointer to the structure to reverse 
  */
void
tmap_seq_reverse(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to reverse compliment
  */
void
tmap_seq_reverse_compliment(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure to compliment
  */
void
tmap_seq_compliment(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure 
  */
void
tmap_seq_to_int(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure 
  @return      a pointer to the name string
  */
tmap_string_t *
tmap_seq_get_name(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure 
  @return      a pointer to the base sequence string
  */
inline tmap_string_t *
tmap_seq_get_bases(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure 
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
