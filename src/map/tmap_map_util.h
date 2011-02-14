/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_UTIL_H
#define TMAP_MAP_UTIL_H

#include <sys/types.h>

#define TMAP_MAP_UTIL_FSW_OFFSET 2
#define TMAP_MAP_UTIL_SCORE_MATCH 5
#define TMAP_MAP_UTIL_PEN_MM 3
#define TMAP_MAP_UTIL_PEN_GAPO 3
#define TMAP_MAP_UTIL_PEN_GAPE 1
#define TMAP_MAP_UTIL_FSCORE 7 

#define TMAP_MAP_UTIL_FLOW_ORDER "TACG"

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

/*!
  The various algorithm types (flags)
  */
enum {
    TMAP_MAP_ALGO_NONE = 0x0,  /*!< dummy algorithm */
    TMAP_MAP_ALGO_MAP1 = 0x1,  /*!< the map1 algorithm */
    TMAP_MAP_ALGO_MAP2 = 0x2,  /*!< the map2 algorithm */
    TMAP_MAP_ALGO_MAP3 = 0x4,  /*!< the map3 algorithm */
    TMAP_MAP_ALGO_MAPALL = 0x8000 /*!< the mapall algorithm */
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
    char *flow; /*!< the flow order (-x) */
    int32_t bw; /*!< the extra bases to add before and after the target during Smith-Waterman (-w) */
    int32_t aln_global; /*!< align the full read (-g) */
    int32_t dup_window; /*!< remove duplicate alignments from different algorithms within this bp window (-W) */
    int32_t score_thr;  /*!< the score threshold (match-score-scaled) (-T) */
    int32_t reads_queue_size;  /*!< the reads queue size (-q) */
    int32_t num_threads;  /*!< the number of threads (-n) */
    int32_t aln_output_mode;  /*!< specifies how to choose alignments (-a)  */
    char *sam_rg;  /*!< specifies the RG line in the SAM header (-R) */
    int32_t sam_sff_tags;  /*!< specifies to output SFF specific SAM tags (-Y) */
    int32_t input_compr;  /*!< the input compression type (-j and -z) */
    int32_t output_compr;  /*!< the output compression type (-J and -Z) */
    key_t shm_key;  /*!< the shared memory key (-k) */

    // map1/map3 options
    int32_t seed_length; /*!< the kmer seed length (-l) */
    int32_t seed_length_set; /*!< 1 if the user has set seed length (-l) */
    
    // map1 options
    int32_t seed_max_mm;  /*!< maximum number of mismatches in the seed (-s) */
    int32_t max_mm;  /*!< maximum number of mismatches (-m) */
    int32_t seed2_length;  /*!< the secondary seed length (-L) */
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
    int32_t max_seed_band; /*!< the band to group seeds (-b)*/
    int32_t hp_diff; /*!< single homopolymer error difference for enumeration (-H) */

    // mapall options
    uint32_t algos[2];  /*!< the algorithms that should be run in stage 1 and stage 2, bit-packed */
    int32_t aln_output_mode_ind; /*!< apply the output filter for each algorithm separately (-I) */
    int32_t num_stages;  /*!< the number of stages */ 
    // stage 1/2 mapping algorithm specific options
    struct __tmap_map_opt_t *opt_map1[2]; /*!< map 1 options */
    struct __tmap_map_opt_t *opt_map2[2]; /*!< map 2 options */
    struct __tmap_map_opt_t *opt_map3[2]; /*!< map 3 options */

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
    uint16_t n_seeds:15; /*!< the number seeds in this hit */
} tmap_map_map3_aux_t;

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
    union {
        tmap_map_map1_aux_t *map1_aux; /*!< auxiliary data for map1 */
        tmap_map_map2_aux_t *map2_aux; /*!< auxiliary data for map2 */
        tmap_map_map3_aux_t *map3_aux; /*!< auxiliary data for map3 */
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
  re-aligns mappings in flow space
  @param  sff          the sff read sequence
  @param  sams         the mappings to adjust 
  @param  refseq       the reference sequence
  @param  bw           the band width
  @param  aln_global   1 to align the full read, 0 otherwise
  @param  score_thr    the alignment score threshold
  @param  score_match  the match score
  @param  pen_mm       the mismatch penalty
  @param  pen_gapo     the gap open penalty
  @param  pen_gape     the gap extension penalty
  @param  fscore       the flow penalty
  */
void
tmap_map_util_fsw(tmap_sff_t *sff, 
                  tmap_map_sams_t *sams, tmap_refseq_t *refseq,
                  int32_t bw, int32_t aln_global, int32_t score_thr,
                  int32_t score_match, int32_t pen_mm, int32_t pen_gapo, 
                  int32_t pen_gape, int32_t fscore);
#endif
