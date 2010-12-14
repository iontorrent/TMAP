#ifndef FMAP_MAP3_AUX_H_
#define FMAP_MAP3_AUX_H_

#include "fmap_map3.h"

/*!
  Holds the seed matches
  */
typedef struct {
    uint32_t k; /*!< the lower SA interval */
    uint32_t l; /*!< the upper SA interval */
    uint16_t start; /*!< the # of bases from the start of the read (0-based) */
    int16_t offset; /*!< the # of bases inserted or deleted from the read */
} fmap_map3_aux_seed_t;

/*! 
  Holds the reference co-ordinate seed matches
 */
typedef struct {
    uint32_t seqid; /*!< the sequence index (0-based) */
    uint32_t pos; /*!< the position (0-based) */
    uint16_t start; /*!< the # of bases from the start of the read (0-based) */
    // int16_t offset; /*!< the # of bases inserted (+) or deleted (-) from the seed */
} fmap_map3_aux_hit_t;

/*!
  Core mapping routine
  @param  seq     the sequence to align (forward/reverse-compliment)
  @param  flow     the flow order (forward/reverse-compliment)
  @param  refseq  the reference sequence structure (forward)
  @param  bwt     the BWT structure (reverse)
  @param  sa      the SA structure (reverse)
  @param  opt     the program options
  @details        the sequences should be in 2-bit format
  @reutrn         the alignments
  */
fmap_map_sams_t *
fmap_map3_aux_core(fmap_seq_t *seq[2],
                   uint8_t *flow[2],
                   fmap_refseq_t *refseq,
                   fmap_bwt_t *bwt,
                   fmap_sa_t *sa,
                   fmap_map3_opt_t *opt);

#endif
