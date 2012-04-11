/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_UTIL_H
#define TMAP_MAP_UTIL_H

#include <sys/types.h>
#include "../../util/tmap_rand.h"
#include "../../sw/tmap_fsw.h"
#include "../../sw/tmap_vsw.h"
#include "tmap_map_opt.h"

#define __map_util_gen_ap(par, opt) do { \
    int32_t i; \
    for(i=0;i<25;i++) { \
        (par).matrix[i] = -(opt)->pen_mm; \
    } \
    for(i=0;i<4;i++) { \
        (par).matrix[i*5+i] = (opt)->score_match; \
    } \
    (par).gap_open = (opt)->pen_gapo; (par).gap_ext = (opt)->pen_gape; \
    (par).gap_end = (opt)->pen_gape; \
    (par).row = 5; \
    (par).band_width = (opt)->bw; \
} while(0)

#define __tmap_map_util_reverse_soft_clipping(_sc) \
  (((_sc) == TMAP_MAP_OPT_SOFT_CLIP_LEFT) ? \
   TMAP_MAP_OPT_SOFT_CLIP_RIGHT : \
   (((_sc) == TMAP_MAP_OPT_SOFT_CLIP_RIGHT) ? \
    TMAP_MAP_OPT_SOFT_CLIP_LEFT : (_sc)))

/*! 
  Auxiliary data for map1
  */
typedef struct {
    uint16_t n_mm;  /*!< the current number of mismatches  */
    uint16_t n_gapo;  /*!< the current number of gap opens */
    uint16_t n_gape;  /*!< the current number of gap extensions */
    uint16_t aln_ref;  /*!< the number of reference bases in the alignment */
    uint32_t num_all_sa;  /*!< the number of hits produced by map1 (though fewer may be reported due to -b) */
} tmap_map_map1_aux_t;

/*! 
  Auxiliary data for map2
  */
typedef struct {
    uint16_t XF:2;  /*!< support for the forward/reverse alignment (1-forward 2-reverse 3-both) */
    uint16_t flag:1; /*!< 0 for non-repetitive hit, 1 for repetitive hit */
    uint16_t XE:13;  /*!< the number of supporting seeds */
    int32_t XI;  /*!< the suffix interval size */
} tmap_map_map2_aux_t;

/*! 
  Auxiliary data for map3
  */
typedef struct {
    void *ptr; // NULL
} tmap_map_map3_aux_t;

/*! 
  Auxiliary data for map4
  */
typedef struct {
    void *ptr; // NULL
} tmap_map_map4_aux_t;

/*! 
  Auxiliary data for map3
  */
typedef struct {
    void *ptr; // NULL
} tmap_map_map_vsw_aux_t;

/*!
  General data structure for holding a mapping; for easy outputting to the SAM format
  */
typedef struct {
    uint32_t algo_id:16; /*!< the algorithm id used to obtain this hit */
    uint32_t algo_stage:16; /*!< the algorithm id used to obtain this hit */
    uint8_t strand; /*!< the strand */
    uint32_t seqid;  /*!< the sequence index (0-based) */
    uint32_t pos; /*!< the position (0-based) */
    int16_t mapq; /*!< the mapping quality */
    int32_t score; /*!< the alignment score */
    int32_t ascore;  /*!< the base alignment score (SFF only) */
    int32_t pscore;  /*!< the pairing base alignment score (pairing only) */
    uint8_t proper_pair:1;  /*!< 0 - if not a proper pair, 1 otherwise */
    double num_stds;  /*!< the number of standard deviations from the mean insert size */
    int16_t pmapq; /*!< the pairing mapping quality */
    int32_t score_subo; /*!< the alignment score of the sub-optimal hit */
    int32_t n_cigar; /*!< the number of cigar operators */
    uint32_t *cigar; /*!< the cigar operator array */
    uint16_t target_len; /*!< internal variable, the target length estimated by the seeding step */ 
    uint16_t n_seeds; /*!< the number seeds in this hit */
    union {
        tmap_map_map1_aux_t *map1_aux; /*!< auxiliary data for map1 */
        tmap_map_map2_aux_t *map2_aux; /*!< auxiliary data for map2 */
        tmap_map_map3_aux_t *map3_aux; /*!< auxiliary data for map3 */
        tmap_map_map4_aux_t *map4_aux; /*!< auxiliary data for map4 */
        tmap_map_map_vsw_aux_t *map_vsw_aux; /*!< auxiliary data for map_vsw */
    } aux;
    // for bounding the alignment with vectorized SW
    tmap_vsw_result_t result; /*!< the VSW boundaries (query/target start/end and scores) */
    uint32_t seed_start; /*!< the start of the seed in genomic coordinates used to map this read */
    uint32_t seed_end; /*!< the end of the seed in genomic coordinates used to map this read */
} tmap_map_sam_t;

/*!
  Stores multiple mappings for a given read
  */
typedef struct {
    int32_t n; /*!< the number of hits */
    int32_t max; /*!< the number of hits before filtering */
    tmap_map_sam_t *sams; /*!< array of hits */
} tmap_map_sams_t;

/*!
  The multi-end record structure
  */
typedef struct {
    tmap_map_sams_t **sams; /*!< the sam records */
    int32_t n; /*!< the number of records (multi-end) */
} tmap_map_record_t;

/*!
  make a copy of src and store it in dest
  @param  dest  the destination record
  @param  src   the source record
  */
inline void
tmap_map_sam_copy(tmap_map_sam_t *dest, tmap_map_sam_t *src);

/*!
  allocates memory for the auxiliary data specific to the algorithm specified by algo_id
  @param  s        the mapping structurem
  */
void
tmap_map_sam_malloc_aux(tmap_map_sam_t *s);

/*!
  destroys auxiliary data for the given mapping structure
  @param  s  the mapping structure
  */
inline void
tmap_map_sam_destroy_aux(tmap_map_sam_t *s);

/*!
  destroys the given mapping structure, including auxiliary data
  @param  s  the mapping structure
  */ 
void
tmap_map_sam_destroy(tmap_map_sam_t *s);

/*!
  allocate memory for an empty mapping structure, with no auxiliary data
  @param   prev  copies over the max from prev
  @return        a pointer to the initialized memory
  */
tmap_map_sams_t *
tmap_map_sams_init(tmap_map_sams_t *prev);

/*!
  reallocate memory for mapping structures; does not allocate auxiliary data
  @param  s  the mapping structure
  @param  n  the new number of mappings
  */
void
tmap_map_sams_realloc(tmap_map_sams_t *s, int32_t n);

/*!
  destroys memory associated with these mappings
  @param  s  the mapping structure
  */
void
tmap_map_sams_destroy(tmap_map_sams_t *s);

/*!
  Initializes a new multi-end mapping structure
  @param  num_ends  the number of ends in this record
  @return  the new multi-end mapping structure
 */
tmap_map_record_t*
tmap_map_record_init(int32_t num_ends);

/*!
  Clones a new multi-end mapping structure
  @param  src  the multi-end mapping structure to clone
  @return  the new multi-end mapping structure
 */
tmap_map_record_t*
tmap_map_record_clone(tmap_map_record_t *src);

/*!
  Merges the mappings of two multi-end mappings 
  @param  src   the multi-end mapping structure destination
  @param  dest  the multi-end mapping structure to merge from
 */
void
tmap_map_record_merge(tmap_map_record_t *dest, tmap_map_record_t *src);

/*!
  Destroys a record structure
  @param  record  the mapping structure
 */
void 
tmap_map_record_destroy(tmap_map_record_t *record);

/*!
  merges src into dest
  @param  dest  the destination mapping structure
  @param  src   the source mapping structure
  */
void
tmap_map_sams_merge(tmap_map_sams_t *dest, tmap_map_sams_t *src);

/*!
  clones src
  @param  src   the source mapping structure
  */
tmap_map_sams_t *
tmap_map_sams_clone(tmap_map_sams_t *src);

/*!
  copies the source into the destination, nullifying the source
  @param  dest  the destination mapping structure
  @param  src   the source mapping structure
  */
void
tmap_map_sam_copy_and_nullify(tmap_map_sam_t *dest, tmap_map_sam_t *src);

/*!
  prints the SAM records
  @param  seq           the original read sequence
  @param  refseq        the reference sequence
  @param  sams          the mappings to print
  @param  end_num       0 if there is no mate, 1 if this is the first fragment, 2 if the this is the last fragment
  @param  mates         the mate's mappings, NULL if there is no mate
  @param  sam_flowspace_tags  1 if SFF specific SAM tags are to be outputted, 0 otherwise
  @param  bidirectional  1 if a bidirectional SAM tag is to be added, 0 otherwise
  @param  seq_eq        1 if the SEQ field is to use '=' symbols, 0 otherwise
  */
void
tmap_map_sams_print(tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_map_sams_t *sams, int32_t end_num, 
                    tmap_map_sams_t *mates, int32_t sam_flowspace_tags, int32_t bidirectional, int32_t seq_eq);

/*!
  keep only the mappings with the given score 
  @param  sams     the mappings to keep
  @param  algo_id  the algorithm identifier
  @param  score    the score to keep
  */
void
tmap_map_util_keep_score(tmap_map_sams_t *sams, int32_t algo_id, int32_t score);

/*!
  filters mappings based on the output mode
  @param  sams             the mappings to filter
  @param  aln_output_mode  the output mode
  @param  rand             the random number generator
  */
void
tmap_map_sams_filter(tmap_map_sams_t *sams, int32_t aln_output_mode, tmap_rand_t *rand);

/*!
  filters mappings based on the output mode
  @param  sams             the mappings to filter
  @param  aln_output_mode  the output mode
  @param  algo_id          the algorithm identifier
  @param  rand             the random number generator
  only filters mappings based on the algorithm id (none process all)
  */
void
tmap_map_sams_filter1(tmap_map_sams_t *sams, int32_t aln_output_mode, int32_t algo_id, tmap_rand_t *rand);

/*!
  filters mappings that pass both the scoring and mapping quality thresholds
  @param  sams       the mappings to filter
  @param  score_thr  the score threshold
  @param  mapq_thr   the mapping quality threshold
  */
void
tmap_map_sams_filter2(tmap_map_sams_t *sams, int32_t score_thr, int32_t mapq_thr);

/*!
  removes duplicate alignments that fall within a given window
  @param  sams        the mappings to adjust 
  @param  dup_window  the window size to cluster mappings
  @param  rand        the random number generator
  */
void
tmap_map_util_remove_duplicates(tmap_map_sams_t *sams, int32_t dup_window, tmap_rand_t *rand);

/*!
 Computes the mapping quality score from a small set of summary statistics.
 @param  seq_len          the sequence length
 @param  n_best           the number of best scores
 @param  best_score       the best score
 @param  n_best_subo      the number of best suboptimal scores
 @param  best_subo_score  the best suboptimal score
 @param  opt              the program parameters
 @return                  the mapping quality
 */
int32_t
tmap_map_util_mapq_score(int32_t seq_len, int32_t n_best, int32_t best_score, int32_t n_best_subo, int32_t best_subo_score, tmap_map_opt_t *opt);

/*!
 Computes the mapping quality from the mappings of multiple algorithms
 @param  sams     the sams to update
 @param  seq_len  the sequence length
 @param  opt      the program parameters
 @return          0 upon success, non-zero otherwise
 */
int32_t
tmap_map_util_mapq(tmap_map_sams_t *sams, int32_t seq_len, tmap_map_opt_t *opt);

/*!
  perform local alignment
  @details              only fills in the score, start and end of the alignments
  @param  refseq        the reference sequence
  @param  sams          the seeded sams
  @param  seqs          the query sequence (forward, reverse compliment, reverse, and compliment)
  @param  rand          the random number generator
  @param  opt           the program parameters
  @return               the locally aligned sams
  */
tmap_map_sams_t *
tmap_map_util_sw_gen_score(tmap_refseq_t *refseq,
                 tmap_map_sams_t *sams,
                 tmap_seq_t **seqs,
                 tmap_rand_t *rand,
                 tmap_map_opt_t *opt);

/*!
  perform local alignment
  @details              generates the cigar after tmap_map_util_sw_gen_score has been called
  @param  refseq        the reference sequence
  @param  sams          the seeded sams
  @param  seqs          the query sequence (forward, reverse compliment, reverse, and compliment)
  @param  opt           the program parameters
  @return               the locally aligned sams
  */
tmap_map_sams_t *
tmap_map_util_sw_gen_cigar(tmap_refseq_t *refseq,
                 tmap_map_sams_t *sams, 
                 tmap_seq_t **seqs,
                 tmap_map_opt_t *opt);

/*!
  re-aligns mappings in flow space
  @param  fs             the flow sequence structure to re-use, NULL otherwise
  @param  seq            the seq read sequence
  @param  flow_order      the flow order
  @param  flow_order_len  the flow order length
  @param  key_seq        the key sequence
  @param  key_seq_len    the key sequence length
  @param  sams           the mappings to adjust 
  @param  refseq         the reference sequence
  @param  bw             the band width
  @param  softclip_type  the soft clip type
  @param  score_thr      the alignment score threshold
  @param  score_match    the match score
  @param  pen_mm         the mismatch penalty
  @param  pen_gapo       the gap open penalty
  @param  pen_gape       the gap extension penalty
  @param  fscore         the flow penalty
  @param  use_flowgram   1 to use the flowgram if available, 0 otherwise
  */
void
tmap_map_util_fsw(tmap_fsw_flowseq_t *fs, tmap_seq_t *seq, 
                  uint8_t *flow_order, int32_t flow_order_len,
                  uint8_t *key_seq, int32_t key_seq_len,
                  tmap_map_sams_t *sams, tmap_refseq_t *refseq,
                  int32_t bw, int32_t softclip_type, int32_t score_thr,
                  int32_t score_match, int32_t pen_mm, int32_t pen_gapo, 
                  int32_t pen_gape, int32_t fscore, int32_t use_flowgram);
#endif
