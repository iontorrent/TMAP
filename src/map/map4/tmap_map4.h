/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP4_H
#define TMAP_MAP4_H

#include <config.h>
#include <sys/types.h>
#include "../util/tmap_map_stats.h"

/*! 
  BWA Fastmap algorithm
  */

/*!
  the query iterator for the bi-directional occurrence search
 */
typedef struct {
    const uint8_t *query; /*!< the query to be searched */
    int32_t start; /*!< the current index into the query */
    int32_t len; /*!< the query length */
    tmap_bwt_smem_intv_vec_t *tmpvec[2]; /*!< temporary memory for occurrence intervals */
    tmap_bwt_smem_intv_vec_t *matches; /*!< the occurence interval matches found by this search */
} tmap_map4_aux_smem_iter_t;

/*!
  Returns the inferred seed length given the reference length
  @param  ref_len  the reference length
  @return          the estimated seed length
  */
int32_t
tmap_map4_get_seed_length(uint64_t ref_len);

/*!
 initializes the mapping routine
 @param  data    pointer to the mapping data pointer
 @param  refseq  the reference sequence
 @param  opt     the program options
 @return         0 if successful, non-zero otherwise
 */
int32_t
tmap_map4_init(void **data, tmap_refseq_t *refseq, tmap_map_opt_t *opt);

/*!
 initializes the mapping routine for a given thread
 @param  data  pointer to the mapping data pointer
 @param  flow_order the flow order
 @param  flow_order_len the flow order length
 @param  key_seq the flow order
 @param  key_seq_len the flow order length
 @param  opt   the program options
 @return       0 if successful, non-zero otherwise
 */
int32_t 
tmap_map4_thread_init(void **data, 
                      uint8_t *flow_order, int32_t flow_order_len,
                      uint8_t *key_seq, int32_t key_seq_len,
                      tmap_map_opt_t *opt);

/*!
 runs the mapping routine for a given thread
 @param  data     pointer to the mapping data pointer
 @param  seqs     the sequence to map (forward, reverse compliment, reverse, compliment)
 @param  index    the reference index
 @param  hash     the occurrence hash
 @param  rand     the random number generator to use
 @param  opt      the program options
 @return          the mappings, NULL otherwise
 */
tmap_map_sams_t*
tmap_map4_thread_map(void **data, tmap_seq_t **seqs, 
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
tmap_map4_thread_cleanup(void **data, tmap_map_opt_t *opt);

/*! 
  main-like function for 'tmap map4'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
tmap_map4_main(int argc, char *argv[]);

#endif
