#ifndef FMAP_MAP2_AUX_H_
#define FMAP_MAP2_AUX_H_

#include "../util/fmap_string.h"
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_sa.h"
#include "fmap_map2_mempool.h"
#include "fmap_map2.h"

#define FMAP_MAP2_MASK_LEVEL 0.90f

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
} fmap_map2_hit_t;

/*! 
  stores alignment hits
  */
typedef struct {
    int32_t n;  /*!< the number of hits */
    int32_t max;  /*!< the maximum memory for the number of hits */
    fmap_map2_hit_t *hits;  /*!< the hits */
    int32_t *n_cigar;  /*!< the number of cigar operations per hit */
    uint32_t **cigar;  /*!< the cigar operations per hit */
} fmap_map2_aln_t;

/*! 
  */
typedef struct {
    uint8_t strand:1;  /*!< the strand */
    uint32_t seqid;  /*!< the zero-based reference contig index */
    uint32_t pos;  /*!< the zero-based reference position */
    uint8_t mapq;  /*!< the mapping quality */
    int32_t n_cigar;  /*!< the number of cigar operaters */
    uint32_t *cigar;  /*!< the cigar operators */
    int32_t AS;  /*!< the alignment score */
    int32_t XS;  /*!< the sub-optimal alignment score */
    uint16_t XF:2;  /*!< support for the forward/reverse alignment (1-forward 2-reverse 3-both) */
    uint16_t XE:14;  /*!< the number of supporting seeds */
    int32_t XI;  /*!< the suffix interval size */
} fmap_map2_sam_entry_t;

/*! 
  stores sam entries to be printed
  @field  num_entries  the number of entries
  @field  entries      the array of entries
  */
typedef struct __fmap_map2_sam_t {
    int32_t num_entries;
    fmap_map2_sam_entry_t *entries;
} fmap_map2_sam_t;

/*! 
  destroys an alignment
  @param  a  pointer to the alignment
  */
void
fmap_map2_aln_destroy(fmap_map2_aln_t *a);

/*! 
  resolves duplicate hits
  @param  bwt  pointer to the bwt structure
  @param  sa   pointer to the suffix array
  @param  b    pointer to the alignment
  @param  IS   the maximum occurrence interval for seeding
  */ 
int32_t
fmap_map2_aux_resolve_duphits(const fmap_bwt_t *bwt, const fmap_sa_t *sa, fmap_map2_aln_t *b, int32_t IS);

/*! 
  initializes a container for sam entries
  @param  n  the number of entries to initialize
  @return    a pointer to the initialized memory
  */
fmap_map2_sam_t *
fmap_map2_sam_init(int32_t n);

/*! 
  resizes a container for sam entries
  @param  sam  pointer to the sam entries structure
  @param  n    the new number of entries 
  @return      a pointer to the re-initialized memory
  */
fmap_map2_sam_t *
fmap_map2_sam_realloc(fmap_map2_sam_t *sam, int32_t n);

/*! 
  destroys a container for sam entries
  @param  sam  pointer to the sam entries structure
  */
void
fmap_map2_sam_destroy(fmap_map2_sam_t *sam);

/*! 
  performs the  BWA-like (long-read) algorithm 
  @param  _opt    pointer to the program parameters
  @param  query   pointer to the query sequence 
  @param  refseq  pointer to the reference sequence structure
  @param  bwt     pointer to the bwt structure
  @param  sa      pointer to the SA structure
  @param  pool    pointer to a global memory pool
  */
fmap_map2_sam_t *
fmap_map2_aux_core(fmap_map2_opt_t *_opt,
                   fmap_seq_t *query,
                   fmap_refseq_t *refseq,
                   fmap_bwt_t *bwt[2],
                   fmap_sa_t *sa[2],
                   fmap_map2_global_mempool_t *pool);

#endif
