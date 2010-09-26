#ifndef FMAP_MAP2_CHAIN_H_
#define FMAP_MAP2_CHAIN_H_

#include "fmap_map2_aux.h"

/*! @header
  @abstract  Chaining Functions for the BWA-like (long-read) Algorithm
  */

/*! @typedef
  @abstract      structure to resolve chaining for Smith-Waterman extension
  @param  tbeg   the lower suffix array interval for the target
  @param  tend   the upper suffix array interval for the target
  @param  qbeg   the lower suffix array interval for the query
  @param  qend   the upper suffix array interval for the query
  @param  flag    the origin of the chain (forward/reverse bwt)
  @param  idx    0-based index within the originating hits
  @param  chain  the chain index; also used as a counter
  */
typedef struct {
    uint32_t tbeg, tend;
    int qbeg, qend;
    uint32_t flag:1, idx:31;
    int chain; // also reuse as a counter
} fmap_map2_chain_t;

/*! @function
  @abstract    filters multiple seeds within a given band
  @param  opt  the function options
  @param  len  the sequence length
  @param  b    pointer to the alignment
  */
void 
fmap_map2_chain_filter(const fmap_map2_opt_t *opt, int len, fmap_map2_aln_t *b[2]);

#endif
