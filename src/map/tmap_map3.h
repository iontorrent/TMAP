/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP3_H
#define TMAP_MAP3_H

#include <config.h>
#include <sys/types.h>

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
 @param  refseq  the reference sequence
 @param  opt     the program options
 @return         0 if successful, non-zero otherwise
 */
int32_t
tmap_map3_init(tmap_refseq_t *refseq, tmap_map_opt_t *opt);

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
 @param  data    pointer to the mapping data pointer
 @param  seq     the sequence to map
 @param  refseq  the reference sequence
 @param  bwt     the bwt structure
 @param  sa      the sa structure
 @param  opt     the program options
 @return         the mappings, NULL otherwise
 */
tmap_map_sams_t*
tmap_map3_thread_map(void **data, tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2], tmap_map_opt_t *opt);

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
