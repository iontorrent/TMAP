/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SAM_H
#define TMAP_SAM_H

#include <config.h>
#ifdef HAVE_SAMTOOLS
#include <bam.h>
#include <sam.h>
#endif
#include "../util/tmap_definitions.h"
#include "../index/tmap_refseq.h"
#include "../io/tmap_file.h"

#ifdef HAVE_SAMTOOLS

/*!
 Structure for holding SAM/BAM records
 */
typedef struct {
    bam1_t *b; /*!< the SAM/BAM structure */
    struct _tmap_sam_io_t *io; /*!< the file from which this SAM was read */
    tmap_string_t *name;  /*!< the name string */
    tmap_string_t *seq;  /*!< the sequence string */
    tmap_string_t *qual;  /*!< the quality string */
    int32_t is_int;  /*!< 1 if the sequence is in integer format, 0 otherwise  */
} tmap_sam_t;

/*! 
  initializes sequence read structure
  @return     pointer to the initialized memory 
  */
inline tmap_sam_t *
tmap_sam_init();

/*! 
  @param  sam  a pointer to the sequence structure
  */
inline void
tmap_sam_destroy(tmap_sam_t *sam);

/*! 
  clones the given sequence read structure
  @param  sam  pointer to the sequence read structure to be copied  
  @return     pointer to the initialized memory 
  */
inline tmap_sam_t *
tmap_sam_clone(tmap_sam_t *sam);

/*! 
  reverse the sequence and reverse the qualities
  @param  sam  a pointer to a sequence structure
  */
void
tmap_sam_reverse(tmap_sam_t *sam);

/*! 
  reverse compliments the sequence and reverse the qualities
  @param  sam  a pointer to a sequence structure
  */
void
tmap_sam_reverse_compliment(tmap_sam_t *sam);

/*! 
  compliments the sequence 
  @param  sam  a pointer to a sequence structure
  */
void
tmap_sam_compliment(tmap_sam_t *sam);

/*! 
  converts bases to integer values
  @param  sam  a pointer to a sequence structure
  */
void
tmap_sam_to_int(tmap_sam_t *sam);

/*! 
  converts bases to character values
  @param  sam  a pointer to a sequence structure
  */
void
tmap_sam_to_char(tmap_sam_t *sam);

/*!
  gets the read's bases
  @param  sam  a pointer to a sequence structure
 */
inline tmap_string_t *
tmap_sam_get_bases(tmap_sam_t *sam);

/*!
  gets the read's qualities
  @param  sam  a pointer to a sequence structure
 */
inline tmap_string_t *
tmap_sam_get_qualities(tmap_sam_t *sam);

/*!
  @param  sam        pointer to the structure to convert
  @param  flow_order a pointer to where the flow order should be stored 
  @return            the length of the flow order
  @details         this is currently not implemented
 */
int32_t
tmap_sam_get_flow_order_int(tmap_sam_t *sam, uint8_t **flow_order);

/*!
  @param  sam      pointer to the structure to convert
  @param  key_seq  pointer to where the key sequence should be stored 
  @return          the length of the key sequence
  @details         this is currently not implemented
 */
int32_t
tmap_sam_get_key_seq_int(tmap_sam_t *sam, uint8_t **key_seq);

/*!
  @param  sam      pointer to the structure to convert
  @param  flowgram  pionter to where the flowgram should be stored
  @param  mem      memory size already allocated to flowgram
  @return          the flowgram length
  @details         this is currently not implemented
 */
int32_t
tmap_sam_get_flowgram(tmap_sam_t *sam, uint16_t **flowgram, int32_t mem);

#endif

#endif
