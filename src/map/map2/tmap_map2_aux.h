/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP2_AUX_H
#define TMAP_MAP2_AUX_H

#include "../../util/tmap_string.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_sa.h"
#include "tmap_map2_mempool.h"
#include "tmap_map2.h"

#define TMAP_MAP2_MASK_LEVEL 0.90f

/*! 
  Auxiliary Functions for BWT-like (long-read) Algorithm
  */

/*! 
  stores an alignment hit
  */
typedef struct {
    uint32_t k;  /*!< the lower suffix array interval, or suffix array position  */
    uint32_t l;  /*!< the upper suffix array interval, or 0 when k is the suffix array position */
    uint32_t flag:18;  /*!< records the origin of the hit (forward/reverse bwt in the 17th/18th bit respectively); the strand in the 5th bit; the first bit stores if the hit was repetitive */
    uint32_t n_seeds:14;  /*!< the number of seeds used in the forward alignment */
    int32_t qlen;  /*!< the length of the query in the alignment */
    int32_t tlen;  /*!< the length of the target in the alignment */
    int32_t G;  /*!< the alignment score */
    int32_t G2;  /*!< the sub-optimal alignment score */
    int32_t beg;  /*!< the beginning of the alignment in the query (0-based) */
    int32_t end;  /*!< the end of the alignment in the query (0-based) */
} tmap_map2_hit_t;

/*! 
  stores alignment hits
  */
typedef struct {
    int32_t n;  /*!< the number of hits */
    int32_t max;  /*!< the memory allocateed for the hits */
    tmap_map2_hit_t *hits;  /*!< the hits */
} tmap_map2_aln_t;

/*!
  allocates an alignment
  @return an initialized alignment
  */
tmap_map2_aln_t*
tmap_map2_aln_init();

/*!
  reallocates an alignment
  @param  a  the alignment to reallocate
  @param  n  the new number of alignments to store
  */
void
tmap_map2_aln_realloc(tmap_map2_aln_t *a, int32_t n);

/*! 
  destroys an alignment
  @param  a  pointer to the alignment
  */
void
tmap_map2_aln_destroy(tmap_map2_aln_t *a);

/*! 
  performs the  BWA-like (long-read) algorithm 
  @param  _opt    pointer to the program parameters
  @param  seqs    pointer to the query sequences (forward, reverse compliment, reverse, compliment) 
  @param  refseq  pointer to the reference sequence structure
  @param  bwt     pointer to the bwt structure
  @param  sa      pointer to the SA structure
  @param  rand    the random number generator
  @param  pool    pointer to a global memory pool
  @return         pointer to the alignment
  */
tmap_map_sams_t *
tmap_map2_aux_core(tmap_map_opt_t *_opt,
                   tmap_seq_t *seqs[4],
                   tmap_refseq_t *refseq,
                   tmap_bwt_t *bwt[2],
                   tmap_sa_t *sa[2],
                   tmap_rand_t *rand,
                   tmap_map2_global_mempool_t *pool);
#endif
