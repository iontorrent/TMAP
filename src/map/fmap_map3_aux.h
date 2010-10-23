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
  Initializes an alignment structure
  @return  a pointer to the initialized alignment
  */
fmap_map3_aln_t *
fmap_map3_aln_init();

/*!
  Destroys this alignment
  @param  aln  a pointer to the alignment to destroy
  */
void
fmap_map3_aln_destroy(fmap_map3_aln_t *aln);

/*!
  Reallocates the number hits within this alignment
  @param  aln  a pointer to the alignment to reallocate
  @param  n    the new number of hits
  @return  a pointer to the re-initialized alignment
  */
void
fmap_map3_aln_realloc(fmap_map3_aln_t *aln, int32_t n);

/*!
  Core mapping routine
  @param  seq     the sequence to align (forward/reverse-compliment)
  @param  refseq  the reference sequence structure (forward)
  @param  bwt     the BWT structure (reverse)
  @param  sa      the SA structure (reverse)
  @param  opt     the program options
  @details        the sequences should be in 2-bit format
  */
fmap_map3_aln_t *
fmap_map3_aux_core(fmap_seq_t *seq[2],
                   fmap_refseq_t *refseq,
                   fmap_bwt_t *bwt,
                   fmap_sa_t *sa,
                   fmap_map3_opt_t *opt);

#endif
