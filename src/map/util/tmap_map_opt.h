/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_OPT_H
#define TMAP_MAP_OPT_H

#include <sys/types.h>
#include <getopt.h>
#include "../../sw/tmap_vsw.h"

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
    TMAP_MAP_ALGO_MAPVSW = 0x400,  /*!< the mapvsw algorithm */
    TMAP_MAP_ALGO_STAGE = 0x800, /*!< the stage options */
    TMAP_MAP_ALGO_MAPALL = 0x1000, /*!< the mapall algorithm */
    TMAP_MAP_ALGO_PAIRING = 0x2000, /*!< flowspace options when printing parameters */
    TMAP_MAP_ALGO_FLOWSPACE = 0x4000, /*!< flowspace options when printing parameters */
    TMAP_MAP_ALGO_GLOBAL = 0x8000, /*!< global options when printing parameters */
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

/*!
  The various option types
  */
enum {
    TMAP_MAP_OPT_TYPE_INT = 0,
    TMAP_MAP_OPT_TYPE_FLOAT,
    TMAP_MAP_OPT_TYPE_NUM,
    TMAP_MAP_OPT_TYPE_FILE,
    TMAP_MAP_OPT_TYPE_STRING,
    TMAP_MAP_OPT_TYPE_NONE
};

/*!
  The print function for the options
  */
typedef void (*tmap_map_opt_option_print_t)(void *opt);

/*!
  The option structure
 */
typedef struct {
    struct option option;
    char *name;
    int32_t type;
    char *description;
    char **multi_options;
    uint32_t algos;
    tmap_map_opt_option_print_t print_func;
} tmap_map_opt_option_t;

/*!
  The options structure
  */
typedef struct {
    tmap_map_opt_option_t *options;
    int32_t n, mem;
    int32_t max_flag_length;
    int32_t max_type_length;
} tmap_map_opt_options_t;

/** 
 * A list of global command line flags take or available.
 *
 * Taken:
 * ABEFGJKLMORTUWXYZ
 * afghijklnqrsvwyz
 *
 * Available:
 * CDHIUV
 * moptux
 * 
 * NB: Lets reserve single character flags for global options. 
*/

/*!
  The command line options structure
 */
typedef struct __tmap_map_opt_t {
    tmap_map_opt_options_t *options;
    int32_t algo_id;
    int32_t algo_stage; /*!< one-based algorithm stage */

    // global options
    char **argv;  /*!< the command line argv structure */
    int argc;  /*!< the number of command line arguments passed */
    char *fn_fasta;  /*!< the fasta reference file name (-f,--fn-fasta) */
    char **fn_reads;  /*!< the reads file name (-r,--fn-fasta) */
    int32_t fn_reads_num; /*!< the number of read files */
    int32_t reads_format;  /*!< the reads file format (-i,--reads-format) */
    char *fn_sam; /*!< the output file name (-s,--fn-sam) */
    int32_t score_match;  /*!< the match score (-A,--score-match) */
    int32_t pen_mm;  /*!< the mismatch penalty (-M,--pen-mismatch) */
    int32_t pen_gapo;  /*!< the indel open penalty (-O,--pen-gap-open) */
    int32_t pen_gape;  /*!< the indel extension penalty (-E,--pen-gap-extension) */
    int32_t bw; /*!< the extra bases to add before and after the target during Smith-Waterman (-w,--band-width) */
    int32_t softclip_type; /*!< soft clip type (-g,--softclip-type) */
    int32_t dup_window; /*!< remove duplicate alignments from different algorithms within this bp window (-W,--duplicate-window) */
    int32_t max_seed_band; /*!< the band to group seeds (-B,--max-seed-band) */
    int32_t no_unroll_banding; /*!< do not unroll the grouped seeds from banding if multiple alignments are found (-U,--no-unroll-banding) */
    int32_t score_thr;  /*!< the score threshold (match-score-scaled) (-T,--score-thres) */
    int32_t reads_queue_size;  /*!< the reads queue size (-q,--reads-queue-size) */
    int32_t num_threads;  /*!< the number of threads (-n,--num-threads) */
    int32_t aln_output_mode;  /*!< specifies how to choose alignments (-a,--aln-output-mode) */
    char *sam_rg;  /*!< specifies the RG line in the SAM header (-R,--sam-read-group) */
    int32_t input_compr;  /*!< the input compression type (-j,--input-bz2 and -z,--input-gz) */
    int32_t output_compr;  /*!< the output compression type (-J,--output-bz2 and -Z,--output-gz) */
    key_t shm_key;  /*!< the shared memory key (-k,--shared-memory-key) */

    // flowspace tags
    int32_t fscore;  /*!< the flow score penalty (-X,--pen-flow-error) */
    char *flow_order; /*!< the flow order (-F,--flow-order) */
    int32_t flow_order_use_file; /*!< the flow order should be from the sff */
    char *key_seq; /*!< the key sequence (-K,--key-sequence) */
    int32_t key_seq_use_file; /*!< the key sequence should be from the sff */
    int32_t softclip_key; /*!< soft clip only the last base of the key (-y,--softclip-key) */
    int32_t sam_sff_tags;  /*!< specifies to output SFF specific SAM tags (-Y,--sam-sff-tags) */
    int32_t ignore_flowgram;  /*!< specifies to ignore the flowgram if available (-S,--ignore-flowgram) */
    int32_t remove_sff_clipping; /*!< removes SFF clipping (-G,--remove-sff-clipping) */

    // pairing options
    int32_t pairing; // TODO
    int32_t strandedness; // TODO
    int32_t positioning; // TODO
    double ins_size_mean; // TODO
    double ins_size_std; // TODO
    double ins_size_std_max_num; // TODO
    int32_t read_rescue; // TODO
    double read_rescue_std_num; // TODO

    // map1/map2/map3 options, but specific to each
    int32_t min_seq_len; /*< the minimum sequence length to examine (--min-seq-length) */
    int32_t max_seq_len; /*< the maximum sequence length to examine (--max-seq-length) */

    // map1/map3 options
    int32_t seed_length; /*!< the kmer seed length (-l) */
    int32_t seed_length_set; /*!< 1 if the user has set seed length (--seed-length) */
    
    // map1 options
    int32_t seed_max_diff;  /*!< maximum number of edits in the seed (--seed-max-diff) */
    int32_t seed2_length;  /*!< the secondary seed length (--seed2-length) */
    int32_t max_diff; /*!< maximum number of edits (--max-diff) */
    double max_diff_fnr; /*!< false-negative probability assuming a maximum error rate (--max-diff-fnr) */ 
    int32_t max_diff_table[TMAP_MAP_OPT_MAX_DIFF_READ_LENGTH+1]; /*!< the maximum number of differences for varying read lengths */
    double max_err_rate; /*!< the maximum error rate (--max-error-rate) */
    int32_t max_mm;  /*!< maximum number of mismatches (--max-mismatches) */
    double max_mm_frac;  /*!< maximum (read length) fraction of mismatches (--max-mismatch-frac) */
    int32_t max_gapo;  /*!< maximum number of indel opens (--max-gap-opens) */
    double max_gapo_frac;  /*!< maximum (read length) fraction of indel opens (--max-gap-open-frac) */
    int32_t max_gape;  /*!< maximum number of indel extensions (--max-gap-extensions) */
    double max_gape_frac;  /*!< maximum fraction of indel extensions (--max-gap-extension-frac) */
    int32_t max_cals_del;  /*!< the maximum number of CALs to extend a deletion (--max-cals-deletion) */
    int32_t indel_ends_bound;  /*!< indels are not allowed within INT number of bps from the end of the read (--indel-ends-bound) */
    int32_t max_best_cals;  /*!< stop searching when INT optimal CALs have been found (--max-best-cals) */
    int32_t max_entries;  /*!< maximum number of alignment nodes (--max-nodes) */
    
    // map2 options
    double yita;  /*!< the error recurrence coefficient (PERMANENTLY SET) */
    //double mask_level;  /*!< the mask level (-m) */
    double length_coef;  /*!< the coefficient of length-threshold adjustment (--length-coef) */
    int32_t max_seed_intv;  /*!< the maximum seed interval (--max-seed-intv) */
    int32_t z_best;  /*!< the number of top scoring hits to keep (--z-best) */
    int32_t seeds_rev;  /*!< the maximum number of seeds for which reverse alignment is triggered (--seeds-rev) */

    // map3 options
    int32_t max_seed_hits; /*!< the maximum number of hits returned by a seed (--max-seed-hits) */
    int32_t hp_diff; /*!< single homopolymer error difference for enumeration (--hp-diff) */
    double hit_frac; /*!< the fraction of seed positions that are under the maximum (--hit-frac) */
    int32_t seed_step; /*!< the number of bases to increase the seed for each seed increase iteration (--seed-step) */ 
    int32_t fwd_search; /*!< perform a forward search instead of a reverse search (--fwd-search) */
    double skip_seed_frac; /*!< the fraction of a seed to skip when a lookup succeeds (--skip-seed-frac) */ 

    // mapvsw options
    // None

    // stage options
    int32_t stage_score_thr;  /*!< the stage one scoring threshold (match-score-scaled) (--stage-score-thres) */
    int32_t stage_mapq_thr;  /*!< the stage one mapping quality threshold (--stage-mapq-thres) */
    int32_t stage_keep_all;  /*!< keep mappings that do not pass the first stage threshold for the next stage (--stage-keep-all) */
    double  stage_seed_freqc; /*!< the minimum frequency a seed must occur in order to be considered for mapping (--stage-seed-freq-cutoff) */

    // sub-options
   struct __tmap_map_opt_t **sub_opts; /*!< sub-options, for multi-stage and multi-mapping */
   int32_t num_sub_opts; /*!< the number of sub-options */
} tmap_map_opt_t;

/*!
  Gets the initialized options
  @return  pointer to the initialized options
  */
tmap_map_opt_t *
tmap_map_opt_init();

/*!
  Add a sub-option to this option list.
  @param  opt  the main option
  @param  algo_id  the algorithm id of the option to create
  @return  a pointer to the sub-option
 */
tmap_map_opt_t*
tmap_map_opt_add_sub_opt(tmap_map_opt_t *opt, int32_t algo_id);

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
  Checks that the global and flowspace options are the same.
  @param  opt_a  the first option
  @param  opt_b  the second option
  */
void
tmap_map_opt_check_global(tmap_map_opt_t *opt_a, tmap_map_opt_t *opt_b); 

/*!
  Checks that the stage options are the same.
  @param  opt_a  the first option
  @param  opt_b  the second option
  */
void
tmap_map_opt_check_stage(tmap_map_opt_t *opt_a, tmap_map_opt_t *opt_b); 

/*!
  Copies only the global parameters from src into opt
  @param  opt_dest  the destination
  @param  opt_src   the soure
  */
void
tmap_map_opt_copy_global(tmap_map_opt_t *opt_dest, tmap_map_opt_t *opt_src);

/*!
  Copies only the stage parameters from src into opt
  @param  opt_dest  the destination
  @param  opt_src   the soure
  */
void
tmap_map_opt_copy_stage(tmap_map_opt_t *opt_dest, tmap_map_opt_t *opt_src);

#endif
