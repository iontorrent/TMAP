/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP2_AUX_H
#define TMAP_MAP2_AUX_H

#include "../util/tmap_string.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_sa.h"
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
    int32_t len;  /*!< the length of the alignment */
    int32_t G;  /*!< the alignment score */
    int32_t G2;  /*!< the sub-optimal alignment score */
    int32_t beg;  /*!< the beginning of the alignment (0-based) */
    int32_t end;  /*!< the end of the alignment (0-based) */
    int32_t n_cigar;  /*!< the number of cigar operations per hit */
    uint32_t *cigar;  /*!< the cigar operations per hit */
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
  resolves duplicate hits
  @param  bwt     pointer to the bwt structure
  @param  sa      pointer to the suffix array
  @param  b       pointer to the alignment
  @param  IS      the maximum occurrence interval for seeding
  @param  min_as  the minimum alignment score to accept a hit
  */ 
int32_t
tmap_map2_aux_resolve_duphits(const tmap_bwt_t *bwt, const tmap_sa_t *sa, tmap_map2_aln_t *b, 
                              int32_t IS, int32_t min_as);

/*! 
  performs the  BWA-like (long-read) algorithm 
  @param  _opt    pointer to the program parameters
  @param  query   pointer to the query sequence 
  @param  refseq  pointer to the reference sequence structure
  @param  bwt     pointer to the bwt structure
  @param  sa      pointer to the SA structure
  @param  pool    pointer to a global memory pool
  @return         pointer to the alignment
  */
tmap_map_sams_t *
tmap_map2_aux_core(tmap_map_opt_t *_opt,
                   tmap_seq_t *query,
                   tmap_refseq_t *refseq,
                   tmap_bwt_t *bwt[2],
                   tmap_sa_t *sa[2],
                   tmap_map2_global_mempool_t *pool);

#endif
