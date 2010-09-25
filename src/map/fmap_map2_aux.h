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

/*! @header
  @abstract Auxiliary Functions for the BWA-like (long-read) Algorithm
  */

/*! @typedef
  @abstract        stores an alignment hit
  @param  k        the lower suffix array interval, or suffix array position 
  @param  l        the upper suffix array interval, or 0 when k is the suffix array position
  @param  flag      records the origin of the hit (forward/reverse bwt in the 17th/18th bit respectively); the strand in the 5th bit; the first bit stores if the hit was repetitive
  @param  n_seeds  the number of seeds used in the forward alignment
  @param  len      the length of the alignment
  @param  G        the alignment score
  @param  G2       the sub-optimal alignment score
  @param  beg      the beginning of the alignment
  @param  end      the end of the alignment
  */
typedef struct {
    uint32_t k, l, flag:18, n_seeds:14;
    int32_t len, G, G2;
    int32_t beg, end;
} fmap_map2_hit_t;

/*! @typedef
  @abstract        stores alignment hits
  @param  n        the number of hits
  @param  max      the maximum memory for the number of hits
  @param  hits     the hits
  @param  n_cigar  the number of cigar operations per hit
  @param  cigar    the cigar operations per hit
  */
typedef struct {
    int32_t n, max;
    fmap_map2_hit_t *hits;
    int32_t *n_cigar;
    uint32_t **cigar;
} fmap_map2_aln_t;

/*! @typedef
  @abstract  
  @param  strand   the strand
  @param  seqid    the zero-based reference contig index
  @param  pos      the zero-based reference position
  @param  mapq     the mapping quality
  @param  n_cigar  the number of cigar operaters
  @param  cigar    the cigar operators
  @param  AS       the alignment score
  @param  XS       the sub-optimal alignment score
  @param  XF       support for the forward/reverse alignment (1-forward 2-reverse 3-both)
  @param  XE       the number of supporting seeds
  @param  XI       the suffix interval size
  */
typedef struct {
    uint8_t strand:1; // 1-bit
    uint32_t seqid, pos; // zero-based
    uint8_t mapq;
    int32_t n_cigar;
    uint32_t *cigar;
    int32_t AS;
    int32_t XS;
    uint16_t XF:2, XE:14;
    int32_t XI;
} fmap_map2_sam_entry_t;

/*! @typedef
  @abstract            stores sam entries to be printed
  @param  num_entries  the number of entries
  @param  entries      the array of entries
  */
typedef struct __fmap_map2_sam_t {
    int32_t num_entries;
    fmap_map2_sam_entry_t *entries;
} fmap_map2_sam_t;

/*! @function
  @abstract  destroys an alignment
  @param  a  pointer to the alignment
  */
void
fmap_map2_aln_destroy(fmap_map2_aln_t *a);

/*! @function
  @abstract  resolves duplicate hits
  @param  bwt  pointer to the bwt structure
  @param  sa   pointer to the suffix array
  @param  b    pointer to the alignment
  @param  IS   the maximum occurrence interval for seeding
  */ 
int32_t
fmap_map2_aux_resolve_duphits(const fmap_bwt_t *bwt, const fmap_sa_t *sa, fmap_map2_aln_t *b, int32_t IS);

/*! @function
  @abstract  initializes a container for sam entries
  @param  n  the number of entries to initialize
  @return    a pointer to the initialized memory
  */
fmap_map2_sam_t *
fmap_map2_sam_init(int32_t n);

/*! @function
  @abstract    resizes a container for sam entries
  @param  sam  pointer to the sam entries structure
  @param  n    the new number of entries 
  @return      a pointer to the re-initialized memory
  */
fmap_map2_sam_t *
fmap_map2_sam_realloc(fmap_map2_sam_t *sam, int32_t n);

/*! @function
  @abstract    destroys a container for sam entries
  @param  sam  pointer to the sam entries structure
  */
void
fmap_map2_sam_destroy(fmap_map2_sam_t *sam);

/*! @function
  @abstract       performs the  BWA-like (long-read) algorithm 
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
