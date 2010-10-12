#ifndef FMAP_MAP2_CHAIN_H_
#define FMAP_MAP2_CHAIN_H_

#include "fmap_map2_aux.h"

/*! 
  Chaining Functions for the BWA-like (long-read) Algorithm
  */

/*! 
  structure to resolve chaining for Smith-Waterman extension
  */
typedef struct {
    uint32_t tbeg;  /*!< the lower suffix array interval for the target */
    uint32_t tend;  /*!< the upper suffix array interval for the target */
    int qbeg;  /*!< the lower suffix array interval for the query */
    int qend;  /*!< the upper suffix array interval for the query */
    uint32_t flag:1;  /*!< the origin of the chain (forward/reverse bwt) */
    uint32_t idx:31;  /*!< 0-based index within the originating hits */
    int chain;  /*!< the chain index; also used as a counter */
} fmap_map2_chain_t;

/*! 
  filters multiple seeds within a given band
  @param  opt  the function options
  @param  len  the sequence length
  @param  b    pointer to the alignment
  */
void 
fmap_map2_chain_filter(const fmap_map2_opt_t *opt, int len, fmap_map2_aln_t *b[2]);

#endif
