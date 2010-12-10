#ifndef FMAP_MAP_UTIL_H_
#define FMAP_MAP_UTIL_H_

#include "fmap_map2_aux.h"
#include "fmap_map3.h"

#define FMAP_MAP_UTIL_FSW_OFFSET 2
#define FMAP_MAP_UTIL_SCORE_MATCH 5
#define FMAP_MAP_UTIL_PEN_MM 3
#define FMAP_MAP_UTIL_PEN_GAPO 3
#define FMAP_MAP_UTIL_PEN_GAPE 1
#define FMAP_MAP_UTIL_FSCORE 7 

/*!
  The various modes to modify the alignment score
  */
enum {
    FMAP_MAP_UTIL_ALN_MODE_UNIQ_BEST      = 0,  /*!< output an alignment only if its score is better than all other alignments */
    FMAP_MAP_UTIL_ALN_MODE_RAND_BEST      = 1,  /*!< output a random best scoring alignment */
    FMAP_MAP_UTIL_ALN_MODE_ALL_BEST       = 2,  /*!< output all the alignments with the best score */
    FMAP_MAP_UTIL_ALN_MODE_ALL            = 3   /*!< Output all alignments */
};

typedef struct {
    uint16_t strand:1; /*!< the strand */
    uint32_t seqid;  /*!< the sequence index (0-based) */
    uint32_t pos; /*!< the position (0-based) */
    int32_t score; /*!< the alignment score */
    int32_t score_subo; /*!< the alignment score of the sub-optimal hit */
    int32_t n_cigar; /*!< the number of cigar operators */
    uint32_t *cigar; /*!< the cigar operator array */
} fmap_map_util_hit_t;

typedef struct {
    int32_t n; /*!< the number of hits */
    fmap_map_util_hit_t *hits; /*!< array of hits */
} fmap_map_util_fsw_aln_t;

inline void
fmap_map_util_map2_fsw(fmap_sff_t *sff, fmap_map2_sam_t *sam, fmap_refseq_t *refseq, fmap_map2_opt_t *opt);

inline void
fmap_map_util_map3_fsw(fmap_sff_t *sff, fmap_map3_aln_t *aln, fmap_refseq_t *refseq, fmap_map3_opt_t *opt);

#endif
