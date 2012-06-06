/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SEQ_H
#define TMAP_SEQ_H

#include <config.h>
#include "tmap_fq.h"
#include "tmap_sff.h"
#include "tmap_sam.h"
#include "../samtools/sam_header.h"

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
    TMAP_SEQ_TYPE_SAM = 2, /*!< SAM input/output */
    TMAP_SEQ_TYPE_BAM = 3 /*!< BAM input/output */
};

/*! 
  */
typedef struct {
    int8_t type;  /*!< the type associated with this structure */
    union {
        tmap_fq_t *fq;  /*!< the pointer to the fastq structure */
        tmap_sff_t *sff;  /*!< the pointer to the sff structure */
        tmap_sam_t *sam;  /*!< the pointer to the SAM/BAM structure */
    } data;
    const char *ks; /*!< the key sequence associated with this structure, NULL if none */
    const char *fo; /*!< the flow order associated with this structure, NULL if none */
    const sam_header_record_t *rg_record; /*!< the read group record associated with this structure, NULL if none */
    const sam_header_record_t *pg_record; /*!< the program group record associated with this structure, NULL if none */
    int32_t fo_start_idx; /*!< the flow order start index associated with this structure, -1 if none */
    uint16_t *flowgram; /*!< the flowgram */
    int32_t flowgram_len; /*!< the flowgram length */
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
  converts bases to character values
  @param  seq  a pointer to a sequence structure
  */
void
tmap_seq_to_char(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure 
  return       0 if the sequence is character format, 1 otherwise
  */
int32_t
tmap_seq_is_int(tmap_seq_t *seq);

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
  @return      the sequence length
  */
inline int32_t
tmap_seq_get_bases_length(tmap_seq_t *seq);

/*! 
  @param  seq  pointer to the structure 
  @return      a pointer to the quality string
  */
inline tmap_string_t *
tmap_seq_get_qualities(tmap_seq_t *seq);

/*! 
  @param                   seq  pointer to the structure 
  @param  remove_clipping  1 if we are to remove clipped sequence, 0 otherwise
  @param  key_seq          the key sequence (integer format) to enforce, NULL otherwise
  @param  key_seq_len      the key sequence length, 0 otherwise
  @return                  0 if the key sequence did not match, 1 otherwise 
  @details                 this will not modify the header
  */
int32_t
tmap_seq_remove_key_sequence(tmap_seq_t *seq, int32_t remove_clipping, uint8_t *key_seq, int32_t key_seq_len);

/*! 
  @param  seq  pointer to the structure to convert
  @return      a pointer to the fq structure
  */
tmap_seq_t *
tmap_seq_sff2fq(tmap_seq_t *seq);

/*!
  Updates the Read Group record pointer, and flow space information.
  @param  seq  pointer to the structure to update
  @param  idx  the file index from which the read was read, corresponding to the ith read group
  @param  header  the SAM Header
  */
void
tmap_seq_update(tmap_seq_t *seq, int32_t idx, sam_header_t *header);

#endif
