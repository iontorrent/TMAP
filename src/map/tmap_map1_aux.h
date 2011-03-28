/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP1_AUX_H
#define TMAP_MAP1_AUX_H

#include "../util/tmap_fibheap.h"
#include "tmap_map1.h"

/*! 
  Auxiliary Functions for the BWA-like (short-read) Mapping Algorithm
  */

/*! 
  */
typedef struct {
    uint32_t score;  /*!< the current alignment score */
    uint16_t n_mm;  /*!< the current number of mismatches  */
    int16_t n_gapo;  /*!< the current number of gap opens */
    int16_t n_gape;  /*!< the current number of gap extensions */
    uint8_t state:7;  /*!< the current state (match/mismatch/insertion/deletion) */
    uint8_t strand:1;  /*!< the strand of the alignment */
    int16_t offset;  /*!< the number of (read) bases used (one-based) */
    int16_t last_diff_offset;  /*!< the last offset of a base difference (mismatch/insertion/deletion) (zero-based) */
    tmap_bwt_match_occ_t match_sa;  /*!< the current SA interval information */
    uint32_t i;  /*!< the zero-based index of this element in the memory pool  */
    int32_t prev_i;  /*!< the zero-based index of the previous element (in the alignment) in the memory pool */
} tmap_map1_aux_stack_entry_t;

typedef struct {
    int32_t n_entries, m_entries;
    tmap_map1_aux_stack_entry_t **entries;
} tmap_map1_aux_bin_t;

/*! 
 Entry stack for searching.
  */
typedef struct {
    tmap_map1_aux_stack_entry_t **entry_pool;  /*!<  the memory pool of entries */
    int32_t entry_pool_length;  /*!< the memory pool length */ 
    int32_t entry_pool_i;  /*!< the next available entry in the memory pool */
    int32_t best_score;  /*!< the best score for any entry in this stack */
    int32_t n_bins;
    tmap_map1_aux_bin_t *bins;
    int32_t n_entries;
} tmap_map1_aux_stack_t;

/*
  Initializes the entry stack
  @return  a pointer to the initialized stack
 */
tmap_map1_aux_stack_t *
tmap_map1_aux_stack_init();

/*
  Destroys the entry stack
  @param  stack  the stack structure to destroy
 */
void
tmap_map1_aux_stack_destroy(tmap_map1_aux_stack_t *stack);

/*! 
  @param  seq         the base sequences (forward/reverse-complimented)
  @param  refseq      the reference sequence structure 
  @param  bwt         the BWT structure 
  @param  sa          the SA structure 
  @param  width       the bounds within the read (forward/reverse)
  @param  seed_width  the bounds within the seed (forward/reverse)
  @param  opt         the program parameters structure
  @param  stack       the stack structure
  @param  seed2_len   the secondary seed length (overides the parameter in opt) 
  @return             pointer to the alignments
  */
tmap_map_sams_t *
tmap_map1_aux_core(tmap_seq_t *seq[2], tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2],
                   tmap_bwt_match_width_t *width[2], tmap_bwt_match_width_t *seed_width[2], tmap_map_opt_t *opt,
                   tmap_map1_aux_stack_t *stack, int32_t seed2_len);

#endif
