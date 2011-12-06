/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP2_H
#define TMAP_MAP2_H

#include <stdlib.h>
#include <sys/types.h>
#include "../util/tmap_map_stats.h"

/*! 
  The BWA-like (long-read) Mapping Algorithm
  */

/*!
 initializes the mapping routine
 @param  data    pointer to the mapping data pointer
 @param  refseq  the reference sequence
 @param  opt     the program options
 @return         0 if successful, non-zero otherwise
 */
int32_t
tmap_map2_init(void **data, tmap_refseq_t *refseq, tmap_map_opt_t *opt);

/*!
 initializes the mapping routine for a given thread
 @param  data  pointer to the mapping data pointer
 @param  opt   the program options
 @return       0 if successful, non-zero otherwise
 */
int32_t 
tmap_map2_thread_init(void **data, tmap_map_opt_t *opt);

/*!
 runs the mapping routine for a given thread
 @param  data     pointer to the mapping data pointer
 @param  seqs     the sequence to map (forward, reverse compliment, reverse, compliment)
 @param  seq_len  the sequence length
 @param  index    the reference index
 @param  rand     the random number generator
 @param  hash     the occurence hash
 @param  opt      the program options
 @return          the mappings, NULL otherwise
 */
tmap_map_sams_t*
tmap_map2_thread_map_core(void **data, tmap_seq_t *seqs[4], int32_t seq_len, tmap_index_t *index, tmap_bwt_match_hash_t *hash[2], tmap_rand_t *rand, tmap_map_opt_t *opt);

// TODO
tmap_map_sams_t*
tmap_map2_thread_map(void **data, tmap_seq_t **seqs, 
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
tmap_map2_thread_cleanup(void **data, tmap_map_opt_t *opt);

/*! 
  main-like function for 'tmap map2'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
tmap_map2_main(int argc, char *argv[]);

#endif
