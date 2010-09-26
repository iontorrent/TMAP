#ifndef FMAP_MAP1_AUX_H_
#define FMAP_MAP1_AUX_H_

#include "../util/fmap_fibheap.h"
#include "fmap_map1.h"

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
    int16_t last_diff_offset;  /*!< the last offset of a base difference (mismatch/insertion/deletion) */
    fmap_bwt_match_occ_t match_sa;  /*!< the current SA interval information */
    uint32_t i;  /*!< the zero-based index of this element in the memory pool  */
    int32_t prev_i;  /*!< the zero-based index of the previous element (in the alignment) in the memory pool */
} fmap_map1_aux_stack_entry_t;

/*! 
  */
typedef struct {
    fmap_map1_aux_stack_entry_t **entry_pool;  /*!<  the memory pool of entries */
    int32_t entry_pool_length;  /*!< the memory pool length */ 
    int32_t entry_pool_i;  /*!< the next available entry in the memory pool */
    fmap_fibheap_t *heap;  /*!< the entry heap for this stack*/
    int32_t best_score;  /*!< the best score for any entry in this stack */
} fmap_map1_aux_stack_t;

// TODO
fmap_map1_aux_stack_t *
fmap_map1_aux_stack_init();

// TODO
void
fmap_map1_aux_stack_destroy(fmap_map1_aux_stack_t *stack);

/*! 
  @param  seq         the base sequences (forward/reverse-complimented)
  @param  bwt         the BWT structure (reversed)
  @param  width       the bounds within the read (forward/reverse)
  @param  seed_width  the bounds within the seed (forward/reverse)
  @param  opt         the program parameters structure
  @param  stack       the stack structure
  @param  n_alns      the number of alignments returned
  @return             pointer to the alignments, NULL terminated
  */
fmap_map1_aln_t **
fmap_map1_aux_core(fmap_seq_t *seq[2], fmap_bwt_t *bwt,
                   fmap_bwt_match_width_t *width[2], fmap_bwt_match_width_t *seed_width[2], fmap_map1_opt_t *opt,
                   fmap_map1_aux_stack_t *stack, int32_t *n_alns);

#endif
