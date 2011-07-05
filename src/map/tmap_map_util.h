/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_UTIL_H
#define TMAP_MAP_UTIL_H

#include <sys/types.h>
#include "../sw/tmap_vsw.h"

#define TMAP_MAP_UTIL_FSW_OFFSET 2
#define TMAP_MAP_UTIL_SCORE_MATCH 1
#define TMAP_MAP_UTIL_PEN_MM 3
#define TMAP_MAP_UTIL_PEN_GAPO 5
#define TMAP_MAP_UTIL_PEN_GAPE 2
#define TMAP_MAP_UTIL_FSCORE 2

#define TMAP_MAP_UTIL_FLOW_ORDER "TACG"
#define TMAP_MAP_UTIL_MAX_DIFF_READ_LENGTH 250

#define __tmap_map_print_compression(_type) switch(_type) { \
  case TMAP_FILE_NO_COMPRESSION: \
                                 tmap_file_fprintf(tmap_file_stderr, " [none]\n"); \
    break; \
  case TMAP_FILE_GZ_COMPRESSION: \
                                 tmap_file_fprintf(tmap_file_stderr, " [gz]\n"); \
    break; \
  case TMAP_FILE_BZ2_COMPRESSION: \
                                  tmap_file_fprintf(tmap_file_stderr, " [bz2]\n"); \
    break; \
  default: \
           tmap_file_fprintf(tmap_file_stderr, " [?]\n"); \
    break; \
}

#define __gen_ap(par, opt) do { \
    int32_t i; \
    for(i=0;i<25;i++) { \
        (par).matrix[i] = -(opt)->pen_mm; \
    } \
    for(i=0;i<4;i++) { \
        (par).matrix[i*5+i] = (opt)->score_match; \
    } \
    (par).gap_open = (opt)->pen_gapo; (par).gap_ext = (opt)->pen_gape; \
    (par).gap_end = (opt)->pen_gape; \
    (par).row = 5; (par).band_width = opt->bw; \
} while(0)

#define __tmap_map_util_reverse_soft_clipping(_sc) \
  (((_sc) == TMAP_MAP_UTIL_SOFT_CLIP_LEFT) ? \
   TMAP_MAP_UTIL_SOFT_CLIP_RIGHT : \
   (((_sc) == TMAP_MAP_UTIL_SOFT_CLIP_RIGHT) ? \
    TMAP_MAP_UTIL_SOFT_CLIP_LEFT : (_sc)))

/*!
  The various algorithm types (flags)
  */
enum {
    TMAP_MAP_ALGO_NONE = 0x0,  /*!< dummy algorithm */
    TMAP_MAP_ALGO_MAP1 = 0x1,  /*!< the map1 algorithm */
    TMAP_MAP_ALGO_MAP2 = 0x2,  /*!< the map2 algorithm */
    TMAP_MAP_ALGO_MAP3 = 0x4,  /*!< the map3 algorithm */
    TMAP_MAP_ALGO_MAP_VSW = 0x4000,  /*!< the mapvsw algorithm */
    TMAP_MAP_ALGO_MAPALL = 0x8000, /*!< the mapall algorithm */
};

/*!
  The various soft-clipping types
  */
enum {
    TMAP_MAP_UTIL_SOFT_CLIP_ALL = 0,  /*!< allow soft-clipping on the right and left portion of the read */
    TMAP_MAP_UTIL_SOFT_CLIP_LEFT = 1,  /*!< allow soft-clipping on the left portion of the read */
    TMAP_MAP_UTIL_SOFT_CLIP_RIGHT = 2,  /*!< allow soft-clipping on the right portion of the read */
    TMAP_MAP_UTIL_SOFT_CLIP_NONE = 3,  /*!< do not soft-clip the read */
};

/*!
  The various modes to modify the alignment score
  */
enum {
    TMAP_MAP_UTIL_ALN_MODE_UNIQ_BEST      = 0,  /*!< output an alignment only if its score is better than all other alignments */
    TMAP_MAP_UTIL_ALN_MODE_RAND_BEST      = 1,  /*!< output a random best scoring alignment */
    TMAP_MAP_UTIL_ALN_MODE_ALL_BEST       = 2,  /*!< output all the alignments with the best score */
    TMAP_MAP_UTIL_ALN_MODE_ALL            = 3   /*!< Output all alignments */
};

typedef struct __tmap_map_opt_t {
    int32_t algo_id;
    int32_t algo_stage;

    // global options
    char **argv;  /*!< the command line argv structure */
    int argc;  /*!< the number of command line arguments passed */
    char *fn_fasta;  /*!< the fasta reference file name (-f) */
    char *fn_reads;  /*!< the reads file name (-r) */
    int32_t reads_format;  /*!< the reads file format (-F)  */
    int32_t score_match;  /*!< the match score (-A) */
    int32_t pen_mm;  /*!< the mismatch penalty (-M) */
    int32_t pen_gapo;  /*!< the indel open penalty (-O) */
    int32_t pen_gape;  /*!< the indel extension penalty (-E) */
    int32_t fscore;  /*!< the flow score penalty (-X) */
    char *flow_order; /*!< the flow order (-x) */
    int32_t flow_order_use_sff; /*!< the flow order should be from the sff (-x) */
    int32_t bw; /*!< the extra bases to add before and after the target during Smith-Waterman (-w) */
    int32_t softclip_type; /*!< soft clip type (-g) */
    int32_t dup_window; /*!< remove duplicate alignments from different algorithms within this bp window (-W) */
    int32_t max_seed_band; /*!< the band to group seeds (-B)*/
    int32_t score_thr;  /*!< the score threshold (match-score-scaled) (-T) */
    int32_t reads_queue_size;  /*!< the reads queue size (-q) */
    int32_t num_threads;  /*!< the number of threads (-n) */
    int32_t aln_output_mode;  /*!< specifies how to choose alignments (-a)  */
    char *sam_rg;  /*!< specifies the RG line in the SAM header (-R) */
    int32_t sam_sff_tags;  /*!< specifies to output SFF specific SAM tags (-Y) */
    int32_t remove_sff_clipping; /*!< removes sff cliping (-G) */
    int32_t input_compr;  /*!< the input compression type (-j and -z) */
    int32_t output_compr;  /*!< the output compression type (-J and -Z) */
    key_t shm_key;  /*!< the shared memory key (-k) */

    // map1/map2/map3 options, but specific to each
    int32_t min_seq_len; /*< the minimum sequence length to examine (-u) */
    int32_t max_seq_len; /*< the maximum sequence length to examine (-U) */

    // map1/map3 options
    int32_t seed_length; /*!< the kmer seed length (-l) */
    int32_t seed_length_set; /*!< 1 if the user has set seed length (-l) */
    
    // map1 options
    int32_t seed_max_diff;  /*!< maximum number of edits in the seed (-s) */
    int32_t seed2_length;  /*!< the secondary seed length (-L) */
    int32_t max_diff; /*!< maximum number of edits (-p) */
    double max_diff_fnr; /*!< false-negative probability assuming a maximum error rate (-p) */ 
    int32_t max_diff_table[TMAP_MAP_UTIL_MAX_DIFF_READ_LENGTH+1]; /*!< the maximum number of differences for varying read lengths */
    double max_err_rate; /*!< the maximum error rate (-P) */
    int32_t max_mm;  /*!< maximum number of mismatches (-m) */
    double max_mm_frac;  /*!< maximum (read length) fraction of mismatches (-m) */
    int32_t max_gapo;  /*!< maximum number of indel opens (-o) */
    double max_gapo_frac;  /*!< maximum (read length) fraction of indel opens (-o) */
    int32_t max_gape;  /*!< maximum number of indel extensions (-e) */
    double max_gape_frac;  /*!< maximum fraction of indel extensions (-e) */
    int32_t max_cals_del;  /*!< the maximum number of CALs to extend a deletion (-d) */
    int32_t indel_ends_bound;  /*!< indels are not allowed within INT number of bps from the end of the read (-i) */
    int32_t max_best_cals;  /*!< stop searching when INT optimal CALs have been found (-b) */
    int32_t max_entries;  /*!< maximum number of alignment nodes (-Q) */
    
    // map2 options
    double yita;  /*!< the error recurrence coefficient (-y)  */
    //double mask_level;  /*!< the mask level (-m) */
    double length_coef;  /*!< the coefficient of length-threshold adjustment (-c) */
    int32_t max_seed_intv;  /*!< the maximum seed interval (-S) */
    int32_t z_best;  /*!< the number of top scoring hits to keep (-b) */
    int32_t seeds_rev;  /*!< the maximum number of seeds for which reverse alignment is triggered (-N) */

    // map3 options
    int32_t max_seed_hits; /*!< the maximum number of hits returned by a seed (-S) */
    int32_t hp_diff; /*!< single homopolymer error difference for enumeration (-H) */

    // mapvsw options
    // None

    // mapall options
    uint32_t algos[2];  /*!< the algorithms that should be run in stage 1 and stage 2, bit-packed */
    int32_t aln_output_mode_ind; /*!< apply the output filter for each algorithm separately (-I) */
    int32_t num_stages;  /*!< the number of stages */ 
    // stage 1/2 mapping algorithm specific options
    struct __tmap_map_opt_t *opt_map1[2]; /*!< map 1 options */
    struct __tmap_map_opt_t *opt_map2[2]; /*!< map 2 options */
    struct __tmap_map_opt_t *opt_map3[2]; /*!< map 3 options */
    struct __tmap_map_opt_t *opt_map_vsw[2]; /*!< map vsw options */
} tmap_map_opt_t;

/*!
  Gets the initialized options
  @return  pointer to the initialized options
  */
tmap_map_opt_t *
tmap_map_opt_init();

/*!
  Destroys the memory associated with these options
  @param  opt  pointer to the options
  */
void
tmap_map_opt_destroy(tmap_map_opt_t *opt);

/*!
  Prints the usage of the map algorithms
  @param  opt  the current options
  @return      always 1
  */
int
tmap_map_opt_usage(tmap_map_opt_t *opt);

/*!
  Parses the command line options and stores them in the options structure
  @param  argc  the number of arguments
  @param  argv  the argument list
  @param  opt   pointer to the options
  @return       1 if successful, 0 otherwise
  */
int32_t
tmap_map_opt_parse(int argc, char *argv[], tmap_map_opt_t *opt);

/*!
  Checks that all options are within range
  @param  opt   pointer to the options
  */
void
tmap_map_opt_check(tmap_map_opt_t *opt);

/*!
  @param  opt   pointer to the options
  */
void
tmap_map_opt_print(tmap_map_opt_t *opt);

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
    uint16_t XE:14;  /*!< the number of supporting seeds */
    int32_t XI;  /*!< the suffix interval size */
} tmap_map_map2_aux_t;

/*! 
  Auxiliary data for map3
  */
typedef struct {
    void *ptr; // NULL
} tmap_map_map3_aux_t;

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
    uint16_t strand:1; /*!< the strand */
    uint32_t seqid;  /*!< the sequence index (0-based) */
    uint32_t pos; /*!< the position (0-based) */
    int16_t mapq; /*!< the mapping quality */
    int32_t score; /*!< the alignment score */
    int32_t ascore;  /*!< the base alignment score (SFF only) */
    int32_t score_subo; /*!< the alignment score of the sub-optimal hit */
    int32_t n_cigar; /*!< the number of cigar operators */
    uint32_t *cigar; /*!< the cigar operator array */
    uint16_t target_len; /*!< internal variable, the target length estimated by the seeding step */ 
    uint16_t n_seeds; /*!< the number seeds in this hit */
    tmap_vsw_result_t *result; // TODO
    union {
        tmap_map_map1_aux_t *map1_aux; /*!< auxiliary data for map1 */
        tmap_map_map2_aux_t *map2_aux; /*!< auxiliary data for map2 */
        tmap_map_map3_aux_t *map3_aux; /*!< auxiliary data for map3 */
        tmap_map_map_vsw_aux_t *map_vsw_aux; /*!< auxiliary data for map_vsw */
    } aux;
} tmap_map_sam_t;

/*!
  Stores multiple mappings for a given read
  */
typedef struct {
    int32_t n; /*!< the number of hits */
    tmap_map_sam_t *sams; /*!< array of hits */
} tmap_map_sams_t;

/*!
  allocates memory for the auxiliary data specific to the algorithm specified by algo_id
  @param  s        the mapping structurem
  @param  algo_id  the algorithm identifier
  */
void
tmap_map_sam_malloc_aux(tmap_map_sam_t *s, int32_t algo_id);

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
  @return  a pointer to the initialized memory
  */
tmap_map_sams_t *
tmap_map_sams_init();

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
  @param  sam_sff_tags  1 if SFF specific SAM tags are to be outputted, 0 otherwise
  */
void
tmap_map_sams_print(tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_map_sams_t *sams, int32_t sam_sff_tags);

/*!
  filters mappings based on the output mode
  @param  sams             the mappings to filter
  @param  aln_output_mode  the output mode
  */
void
tmap_map_sams_filter(tmap_map_sams_t *sams, int32_t aln_output_mode);

/*!
  filters mappings based on the output mode
  @param  sams             the mappings to filter
  @param  aln_output_mode  the output mode
  @param  algo_id          the algorithm identifier
  only filters mappings based on the algorithm id (none process all)
  */
void
tmap_map_sams_filter1(tmap_map_sams_t *sams, int32_t aln_output_mode, int32_t algo_id);

/*!
  removes duplicate alignments that fall within a given window
  @param  sams        the mappings to adjust 
  @param  dup_window  the window size to cluster mappings
  */
void
tmap_map_util_remove_duplicates(tmap_map_sams_t *sams, int32_t dup_window);

/*!
 Computes the mapping quality from the mappings of multiple algorithms
 @param  sams     the sams to update
 @param  seq_len  the sequence length
 @param  opt      the program parameters
 */
inline void
tmap_map_util_mapq(tmap_map_sams_t *sams, int32_t seq_len, tmap_map_opt_t *opt);

/*!
  perform local alignment
  @details              only fills in the score, start and end of the alignments
  @param  refseq        the reference sequence
  @param  sams          the seeded sams
  @param  seq           the query sequence
  @param  opt           the program parameters
  @return               the locally aligned sams
  */
tmap_map_sams_t *
tmap_map_util_sw_gen_score(tmap_refseq_t *refseq,
                 tmap_map_sams_t *sams,
                 tmap_seq_t *seq,
                 tmap_map_opt_t *opt);

/*!
  perform local alignment
  @details              generates the cigar after tmap_map_util_sw_gen_score has been called
  @param  refseq        the reference sequence
  @param  sams          the seeded sams
  @param  seq           the query sequence
  @param  opt           the program parameters
  @return               the locally aligned sams
  */
tmap_map_sams_t *
tmap_map_util_sw_gen_cigar(tmap_refseq_t *refseq,
                 tmap_map_sams_t *sams, 
                 tmap_seq_t *seq,
                 tmap_map_opt_t *opt);

/*!
  @param  sam            the sam record in which to store the results
  @param  target         the target sequence (forward strand)
  @param  target_length  the target sequence length
  @param  query          the query sequence (forward strand)
  @param  query_length   the query sequence length
  @param  seqid          the target sequence id (zero-based)
  @param  pos            the target position (zero-based)
  @param  par            the alignment parameters
  @param  path           the alignment path (must be pre-allocated)
  @param  path_len       the alignment path length
  @param  score_thr      the scoring threshold
  @param  softclip_type  the soft-clipping type
  @param  strand         the strand of the hit
  */
int32_t
tmap_map_util_sw_aux(tmap_map_sam_t *sam,
                 uint8_t *target, int32_t target_length,
                 uint8_t *query, int32_t query_length,
                 uint32_t seqid, uint32_t pos,
                 tmap_sw_param_t *par, tmap_sw_path_t *path, int32_t *path_len,
                 int32_t score_thr, int32_t softclip_type, int32_t strand);

/*!
  re-aligns mappings in flow space
  @param  seq            the seq read sequence
  @param  flow_order      the flow order
  @param  flow_order_len  the flow order length
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
  */
void
tmap_map_util_fsw(tmap_seq_t *seq, 
                  uint8_t *flow_order, int32_t flow_order_len,
                  tmap_map_sams_t *sams, tmap_refseq_t *refseq,
                  int32_t bw, int32_t softclip_type, int32_t score_thr,
                  int32_t score_match, int32_t pen_mm, int32_t pen_gapo, 
                  int32_t pen_gape, int32_t fscore);
#endif
