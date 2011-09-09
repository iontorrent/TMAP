/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_OPT_H
#define TMAP_MAP_OPT_H

#include <sys/types.h>
#include "../sw/tmap_vsw.h"

/*!
  The default offset for homopolymer errors.
  */
#define TMAP_MAP_OPT_FSW_OFFSET 2
/*!
  The default match score.
  */
#define TMAP_MAP_OPT_SCORE_MATCH 1
/*!
  The default mismatch penalty.
  */
#define TMAP_MAP_OPT_PEN_MM 3
/*!
  The default gap open penalty.
  */
#define TMAP_MAP_OPT_PEN_GAPO 5
/*!
  The default gap extension penalty.
  */
#define TMAP_MAP_OPT_PEN_GAPE 2
/*!
  The default offset for homopolymer errors.
  */
#define TMAP_MAP_OPT_FSCORE 2

/*!
  The default flow order.
  */
#define TMAP_MAP_OPT_FLOW_ORDER "TACG"

/*!
  The maximum read length to consider for mapping differences in map1.
  */
#define TMAP_MAP_OPT_MAX_DIFF_READ_LENGTH 250

/*!
  Prints the compression for the input/output.
  @param  _type  the compress type (integer).
  @return  the compression string.
 */
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

/*!
  The various algorithm types (flags)
  */
enum {
    TMAP_MAP_ALGO_NONE = 0x0,  /*!< dummy algorithm */
    TMAP_MAP_ALGO_MAP1 = 0x1,  /*!< the map1 algorithm */
    TMAP_MAP_ALGO_MAP2 = 0x2,  /*!< the map2 algorithm */
    TMAP_MAP_ALGO_MAP3 = 0x4,  /*!< the map3 algorithm */
    TMAP_MAP_ALGO_MAPVSW = 0x4000,  /*!< the mapvsw algorithm */
    TMAP_MAP_ALGO_MAPALL = 0x8000, /*!< the mapall algorithm */
};

/*!
  The various soft-clipping types
  */
enum {
    TMAP_MAP_OPT_SOFT_CLIP_ALL = 0,  /*!< allow soft-clipping on the right and left portion of the read */
    TMAP_MAP_OPT_SOFT_CLIP_LEFT = 1,  /*!< allow soft-clipping on the left portion of the read */
    TMAP_MAP_OPT_SOFT_CLIP_RIGHT = 2,  /*!< allow soft-clipping on the right portion of the read */
    TMAP_MAP_OPT_SOFT_CLIP_NONE = 3,  /*!< do not soft-clip the read */
};

/*!
  The various modes to modify the alignment score
  */
enum {
    TMAP_MAP_OPT_ALN_MODE_UNIQ_BEST      = 0,  /*!< output an alignment only if its score is better than all other alignments */
    TMAP_MAP_OPT_ALN_MODE_RAND_BEST      = 1,  /*!< output a random best scoring alignment */
    TMAP_MAP_OPT_ALN_MODE_ALL_BEST       = 2,  /*!< output all the alignments with the best score */
    TMAP_MAP_OPT_ALN_MODE_ALL            = 3   /*!< Output all alignments */
};
    
/** 
 * A list of command line flags take or available.
 * Taken:
 * abcdefghijklmnopqrstuvwxyz
 * ABCDEFHGIJKLMNOPQRSTUWXYZ
 * 012
 *
 * Available:
 * 3456789
 * 
 * 
*/

typedef struct __tmap_map_opt_t {
    int32_t algo_id;
    int32_t algo_stage;

    // global options
    char **argv;  /*!< the command line argv structure */
    int argc;  /*!< the number of command line arguments passed */
    char *fn_fasta;  /*!< the fasta reference file name (-f) */
    char **fn_reads;  /*!< the reads file name (-r) */
    int32_t fn_reads_num; /*!< the number of read files (-r) */
    int32_t reads_format;  /*!< the reads file format (-F)  */
    char *fn_sam; /*!< the output file name (-0) */
    int32_t score_match;  /*!< the match score (-A) */
    int32_t pen_mm;  /*!< the mismatch penalty (-M) */
    int32_t pen_gapo;  /*!< the indel open penalty (-O) */
    int32_t pen_gape;  /*!< the indel extension penalty (-E) */
    int32_t fscore;  /*!< the flow score penalty (-X) */
    char *flow_order; /*!< the flow order (-x) */
    int32_t flow_order_use_sff; /*!< the flow order should be from the sff (-x) */
    char *key_seq; /*!< the key sequence (-t) */
    int32_t key_seq_use_sff; /*!< the key sequence should be from the sff (-t) */
    int32_t bw; /*!< the extra bases to add before and after the target during Smith-Waterman (-w) */
    int32_t softclip_type; /*!< soft clip type (-g) */
    int32_t softclip_key; /*!< soft clip only the last base of the key (-y) */
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
    int32_t max_diff_table[TMAP_MAP_OPT_MAX_DIFF_READ_LENGTH+1]; /*!< the maximum number of differences for varying read lengths */
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
    double yita;  /*!< the error recurrence coefficient (PERMANENTLY SET) */
    //double mask_level;  /*!< the mask level (-m) */
    double length_coef;  /*!< the coefficient of length-threshold adjustment (-c) */
    int32_t max_seed_intv;  /*!< the maximum seed interval (-S) */
    int32_t z_best;  /*!< the number of top scoring hits to keep (-b) */
    int32_t seeds_rev;  /*!< the maximum number of seeds for which reverse alignment is triggered (-N) */

    // map3 options
    int32_t max_seed_hits; /*!< the maximum number of hits returned by a seed (-S) */
    int32_t hp_diff; /*!< single homopolymer error difference for enumeration (-H) */
    double hit_frac; /*!<   the fraction of seed positions that are under the maximum (-S) (-V) */

    // mapvsw options
    // None

    // mapall options
    uint32_t algos[2];  /*!< the algorithms that should be run in stage 1 and stage 2, bit-packed */
    int32_t aln_output_mode_ind; /*!< apply the output filter for each algorithm separately (-I) */
    int32_t num_stages;  /*!< the number of stages */ 
    int32_t mapall_score_thr;  /*!< the stage one scoring threshold (match-score-scaled) (-C) */
    int32_t mapall_mapq_thr;  /*!< the stage one mapping quality threshold (-D) */
    int32_t mapall_keep_all;  /*!< keep mappings that do not pass the first stage threshold for the next stage (-K) */
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

#endif
