#ifndef FMAP_MAP_UTIL_H_
#define FMAP_MAP_UTIL_H_

/*
#include "fmap_map1.h"
#include "fmap_map2_aux.h"
#include "fmap_map3.h"
*/

#define FMAP_MAP_UTIL_FSW_OFFSET 2
#define FMAP_MAP_UTIL_SCORE_MATCH 5
#define FMAP_MAP_UTIL_PEN_MM 3
#define FMAP_MAP_UTIL_PEN_GAPO 3
#define FMAP_MAP_UTIL_PEN_GAPE 1
#define FMAP_MAP_UTIL_FSCORE 7 

#define FMAP_MAP_UTIL_FLOW_ORDER "TACG"

enum {
    FMAP_MAP_ALGO_NONE = 0x0,  /*!< dummy algorithm */
    FMAP_MAP_ALGO_MAP1 = 0x1,  /*!< the map1 algorithm */
    FMAP_MAP_ALGO_MAP2 = 0x2,  /*!< the map2 algorithm */
    FMAP_MAP_ALGO_MAP3 = 0x4,  /*!< the map3 algorithm */
};

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
    uint16_t n_mm;  /*!< the current number of mismatches  */
    uint16_t n_gapo;  /*!< the current number of gap opens */
    uint16_t n_gape;  /*!< the current number of gap extensions */
} fmap_map_map1_aux_t;

typedef struct {
    uint16_t XF:2;  /*!< support for the forward/reverse alignment (1-forward 2-reverse 3-both) */
    uint16_t XE:14;  /*!< the number of supporting seeds */
    int32_t XI;  /*!< the suffix interval size */
} fmap_map_map2_aux_t;

typedef struct {
    uint16_t n_seeds:15; /*!< the number seeds in this hit */
} fmap_map_map3_aux_t;

typedef struct {
    uint32_t algo_id:16; /*!< the algorithm id used to obtain this hit */
    uint32_t algo_stage:16; /*!< the algorithm id used to obtain this hit */
    uint16_t strand:1; /*!< the strand */
    uint32_t seqid;  /*!< the sequence index (0-based) */
    uint32_t pos; /*!< the position (0-based) */
    int16_t mapq; /*!< the mapping quality */
    int32_t score; /*!< the alignment score */
    int32_t score_subo; /*!< the alignment score of the sub-optimal hit */
    int32_t n_cigar; /*!< the number of cigar operators */
    uint32_t *cigar; /*!< the cigar operator array */
    union {
        fmap_map_map1_aux_t *map1_aux; /*!< auxiliary data for map1 */
        fmap_map_map2_aux_t *map2_aux; /*!< auxiliary data for map2 */
        fmap_map_map3_aux_t *map3_aux; /*!< auxiliary data for map3 */
    } aux;
} fmap_map_sam_t;

typedef struct {
    int32_t n; /*!< the number of hits */
    fmap_map_sam_t *sams; /*!< array of hits */
} fmap_map_sams_t;

void
fmap_map_sam_malloc_aux(fmap_map_sam_t *s, int32_t algo_id);

inline void
fmap_map_sam_destroy_aux(fmap_map_sam_t *s);

void
fmap_map_sam_destroy(fmap_map_sam_t *s);

fmap_map_sams_t *
fmap_map_sams_init();

void
fmap_map_sams_realloc(fmap_map_sams_t *s, int32_t n);

void
fmap_map_sams_destroy(fmap_map_sams_t *s);

void
fmap_map_sam_copy_and_nullify(fmap_map_sam_t *dest, fmap_map_sam_t *src);

void
fmap_map_sams_print(fmap_seq_t *seq, fmap_refseq_t *refseq, fmap_map_sams_t *sams, int32_t sam_sff_tags);

void
fmap_map_sams_filter(fmap_map_sams_t *sams, int32_t aln_output_mode);

void
fmap_map_sams_filter1(fmap_map_sams_t *sams, int32_t aln_output_mode, int32_t algo_id);

void
fmap_map_util_map1_adjust_score(fmap_map_sams_t *sams, int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape);

void
fmap_map_util_fsw(fmap_sff_t *sff, 
                  fmap_map_sams_t *sams, fmap_refseq_t *refseq,
                  int32_t bw, int32_t aln_global, int32_t score_thr,
                  int32_t score_match, int32_t pen_mm, int32_t pen_gapo, 
                  int32_t pen_gape, int32_t fscore);
#endif
