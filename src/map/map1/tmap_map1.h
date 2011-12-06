/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP1_H
#define TMAP_MAP1_H

#include <config.h>
#include <sys/types.h>
#include "../util/tmap_map_stats.h"

/*! 
  BWA-like (short-read) Mapping Algorithm
  */

/*!
 Calculates the maximum number of differences allowed given an error rate
 and false negative rate threshold
 @param  l      the read length
 @param  err    the maximum error rate to tolerate
 @param  thres  the false negative (mapping) rate threshold
 @return        the maximum number of differences
 */
int32_t
tmap_map1_cal_maxdiff(int32_t l, double err, double thres);

/*!
 Prints the number of differences for various read lengths
 @param  opt     the program options
 @param  stage   the mapping stage (-1 for no stage)
 */
void
tmap_map1_print_max_diff(tmap_map_opt_t *opt, int32_t stage);

/*!
 initializes the mapping routine
 @param  data    pointer to the mapping data pointer
 @param  refseq  the reference sequence
 @param  opt     the program options
 @return         0 if successful, non-zero otherwise
 */
int32_t
tmap_map1_init(void **data, tmap_refseq_t *refseq, tmap_map_opt_t *opt);

/*!
 initializes the mapping routine for a given thread
 @param  data  pointer to the mapping data pointer
 @param  opt   the program options
 @return       0 if successful, non-zero otherwise
 */
int32_t 
tmap_map1_thread_init(void **data, tmap_map_opt_t *opt);

/*!
 runs the mapping routine for a given thread
 @param  data     pointer to the mapping data pointer
 @param  seqs     the sequence to map (forward, reverse compliment, reverse, compliment)
 @param  index    the reference index
 @param  rand     the random number generator to use
 @param  hash     the occurrence hash
 @param  opt      the program options
 @return          the mappings, NULL otherwise
 */
tmap_map_sams_t*
tmap_map1_thread_map(void **data, tmap_seq_t **seqs, 
                     tmap_index_t *index, 
                     tmap_bwt_match_hash_t *hash,
                     tmap_rand_t *rand, 
                     tmap_map_opt_t *opt);

/*!
 cleans up the mapping routine for a given thread
 @param  data  pointer to the mapping data pointer
 @param  opt   the program options
 @return       0 if successful, non-zero otherwise
 */
int32_t
tmap_map1_thread_cleanup(void **data, tmap_map_opt_t *opt);

/*! 
  main-like function for 'tmap map1'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
tmap_map1_main(int argc, char *argv[]);

#endif
