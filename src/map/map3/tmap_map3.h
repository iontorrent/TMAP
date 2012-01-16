/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP3_H
#define TMAP_MAP3_H

#include <config.h>
#include <sys/types.h>
#include "../util/tmap_map_stats.h"

/*! 
  SSAHA2-like Mapping Algorithm
  */


/*!
  Structure for a final hit
  */ 
typedef struct {
    uint16_t strand:1; /*!< the strand */
    uint32_t seqid;  /*!< the sequence index (0-based) */
    uint32_t pos; /*!< the position (0-based) */
    int32_t score; /*!< the alignment score */
    int32_t score_subo; /*!< the alignment score of the sub-optimal hit */
    uint8_t mapq; /*!< the mapping quality */
    uint16_t n_seeds:15; /*!< the number seeds in this hit */
    int32_t n_cigar; /*!< the number of cigar operators */
    uint32_t *cigar; /*!< the cigar operator array */
} tmap_map3_hit_t;

/*!
  Returns the inferred seed length given the reference length
  @param  ref_len  the reference length
  @return          the estimated seed length
  */
int32_t
tmap_map3_get_seed_length(uint64_t ref_len);

/*!
 initializes the mapping routine
 @param  data    pointer to the mapping data pointer
 @param  refseq  the reference sequence
 @param  opt     the program options
 @return         0 if successful, non-zero otherwise
 */
int32_t
tmap_map3_init(void **data, tmap_refseq_t *refseq, tmap_map_opt_t *opt);

/*!
 initializes the mapping routine for a given thread
 @param  data  pointer to the mapping data pointer
 @param  opt   the program options
 @return       0 if successful, non-zero otherwise
 */
int32_t 
tmap_map3_thread_init(void **data, tmap_map_opt_t *opt);

/*!
 runs the mapping routine for a given thread
 @param  data     pointer to the mapping data pointer
 @param  seqs     the sequence to map (forward and reverse compliment)
 @param  seq_len  the sequence length
 @param  index    the reference index
 @param  hash     the occurence hash
 @param  opt      the program options
 @return          the mappings, NULL otherwise
 */
tmap_map_sams_t*
tmap_map3_thread_map_core(void **data, tmap_seq_t *seqs[2], int32_t seq_len, tmap_index_t *index, tmap_bwt_match_hash_t *hash[2], tmap_map_opt_t *opt);

/*!
 runs the mapping routine for a given thread
 @param  data     pointer to the mapping data pointer
 @param  seqs     the sequence to map (forward/reverse compliment/compliment/reverse in integer format)
 @param  index    the reference index
 @param  stat     the alignment statistics
 @param  rand     the random number generator
 @param  hash     the occurence hash (forward/reverse)
 @param  opt      the program options
 @return          the mappings, NULL otherwise
 */
tmap_map_sams_t*
tmap_map3_thread_map(void **data, tmap_seq_t **seqs, 
                     tmap_index_t *index, tmap_map_stats_t *stat, tmap_rand_t *rand, 
                     tmap_bwt_match_hash_t *hash[2],
                     tmap_map_opt_t *opt);

/*!
 cleans up the mapping routine for a given thread
 @param  data  pointer to the mapping data pointer
 @param  opt   the program options
 @return       0 if successful, non-zero otherwise
 */
int32_t
tmap_map3_thread_cleanup(void **data, tmap_map_opt_t *opt);

/*! 
  main-like function for 'tmap map3'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
tmap_map3_main(int argc, char *argv[]);

#endif
