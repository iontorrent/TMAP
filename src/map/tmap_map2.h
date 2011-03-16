/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP2_H
#define TMAP_MAP2_H

#include <stdlib.h>
#include <sys/types.h>

/*! 
  The BWA-like (long-read) Mapping Algorithm
  */

/*!
 initializes the mapping routine
 @param  refseq  the reference sequence
 @param  opt     the program options
 @return         0 if successful, non-zero otherwise
 */
int32_t
tmap_map2_init(tmap_refseq_t *refseq, tmap_map_opt_t *opt);

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
 @param  data    pointer to the mapping data pointer
 @param  seq     the sequence to map
 @param  refseq  the reference sequence
 @param  bwt     the bwt structure
 @param  sa      the sa structure
 @param  opt     the program options
 @return         the mappings, NULL otherwise
 */
tmap_map_sams_t*
tmap_map2_thread_map(void **data, tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2], tmap_map_opt_t *opt);

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
