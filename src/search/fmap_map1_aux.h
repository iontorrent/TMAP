#ifndef FMAP_MAP1_AUX_H_
#define FMAP_MAP1_AUX_H_

#include "../util/fmap_fibheap.h"
#include "fmap_map1.h"

/*! @typedef 
  @field  score             the current alignment score
  @field  n_mm              the current number of mismatches 
  @field  n_gapo            the current number of gap opens
  @field  n_gape            the current number of gap extensions
  @field  state             the current state (match/mismatch/insertion/deletion)
  @field  strand            the strand of the alignment
  @field  offset            the zero-based offset from the start of the read
  @field  last_diff_offset  the last offset of a base difference (mismatch/insertion/deletion)
  @field  k                 the lower range of the SA interval
  @field  l                 the upper range of the SA interval
  @field  i                 the zero-based index of this element in the memory pool 
  @field  prev_i            the zero-based index of the previous element (in the alignment) in the memory pool
*/
typedef struct {
    uint32_t score;
    uint32_t n_mm:9, n_gapo:10, n_gape:10, state:2, strand:1;
    int16_t offset, last_diff_offset;
    uint32_t k, l; // SA interval
    uint32_t i, prev_i;
} fmap_map1_aux_stack_entry_t;

/*! @typedef
  @field  entry_pool
  @field  entry_pool_length
  @field  entry_pool_iu
  @field  heap
  @field  best_score
 */
typedef struct {
    fmap_map1_aux_stack_entry_t *entry_pool; // memory pool
    int32_t entry_pool_length; // number of entries in the memory pool
    int32_t entry_pool_i; // 0-based index into the memory pool
    fmap_fibheap_t *heap;
    int32_t best_score;
} fmap_map1_aux_stack_t;

// TODO
fmap_map1_aux_stack_t *
fmap_map1_aux_stack_init();

// TODO
void
fmap_map1_aux_stack_destroy(fmap_map1_aux_stack_t *stack);

/*! @function
  @abstract
  @param  seq         the base sequences (forward/reverse-complimented)
  @param  bwt         the BWT structure (reversed)
  @param  width       the bounds within the read (forward/reverse)
  @param  seed_width  the bounds within the seed (forward/reverse)
  @param  opt         the program parameters structure
  @param  stack       the stack structure
  @return             pointer to the alignments, NULL terminated
  */
fmap_map1_aln_t **
fmap_map1_aux_core(fmap_seq_t *seq[2], fmap_bwt_t *bwt,
                   fmap_bwt_match_width_t *width[2], fmap_bwt_match_width_t *seed_width[2], fmap_map1_opt_t *opt,
                   fmap_map1_aux_stack_t *stack);

#endif
