/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SEQS_H
#define TMAP_SEQS_H

#include <config.h>
#include "tmap_fq.h"
#include "tmap_sff.h"
#include "tmap_sam.h"
#include "tmap_seq.h"

/*! 
  An Abstract Library for DNA Sequence Data
  */


/*! 
  */
typedef struct {
    tmap_seq_t **seqs; // TODO
    int32_t n; // TODO
    int32_t m; // TODO
    int32_t type; // TODO
} tmap_seqs_t;

/*! 
  @param  type  the type associated with this structure
  @return       pointer to the initialized memory 
  */
tmap_seqs_t *
tmap_seqs_init(int8_t type);

/*! 
  @param  seqs  pointer to the structure
  */
void
tmap_seqs_destroy(tmap_seqs_t *seqs);

/*! 
  @param  seqs  pointer to the structure to clone
  @return      pointer to the initialized memory 
  */
tmap_seqs_t *
tmap_seqs_clone(tmap_seqs_t *seqs);

/*!
  @param  seqs  pointer to the structure
  @details     performs a shallow copy
  */
void
tmap_seqs_add(tmap_seqs_t *seqs, tmap_seq_t *seq);

/*!
  @param  seqs  pointer to the structure
  @param  i     the ith end to to get (memory bound)
  @return  the sequence structure, or NULL if none available
  */
tmap_seq_t *
tmap_seqs_get(tmap_seqs_t *seqs, int32_t i);

#endif
