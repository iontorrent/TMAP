/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_INDEX_H
#define TMAP_INDEX_H

#define TMAP_INDEX_LARGE_GENOME 0x1000000
// (2^32) - 1
#define TMAP_INDEX_TOO_BIG_GENOME 0xFFFFFFFF

#include "../server/tmap_shm.h"

/*! 
  @details  Constructs the packed reference sequence, BWT string, and Suffix Array.
  */

/*!
  The reference sequence data structure.
 */
typedef struct {
    tmap_refseq_t *refseq; /*!< the packed reference sequence */
    tmap_bwt_t *bwt[2]; /*!< the forward and reverse FM-indexes */
    tmap_sa_t *sa[2]; /*!< the forward and reverse suffix arrays */
    tmap_shm_t *shm; /*!< the shared memory location if loaded from shared memory */
    key_t shm_key; /*!< the shared memory key, zero if not loaded from shared memory */
} tmap_index_t;

/*!
  Initializes the full reference data from file or shared memory.
  @param  fn_fasta  the FASTA file name
  @param  shm_key   the shared memory key, or zero if we are to read in from file
  @return           the full reference index
 */
tmap_index_t*
tmap_index_init(const char *fn_fasta, key_t shm_key);

/*!
  Destroys the index data.
  @param  index  the reference index to destroy
 */
void
tmap_index_destroy(tmap_index_t *index);

/*! 
  structure to store the command line options for 'tmap index'
  */
typedef struct {
    char *fn_fasta;  /*!< the fasta file name (-f) */
    int32_t occ_interval;  /*!< the occurrence array interval (-o) */
    int32_t hash_width;  /*!< the occurrence hash width (-w) */
    int32_t sa_interval;  /*!< the suffix array interval (-i) */
    int32_t is_large;  /*!< 0 to use the short BWT construction algorith, 1 otherwise (large BWT construction algorithm) */
} tmap_index_opt_t;

/*! 
  main-like function for 'tmap index'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int tmap_index(int argc, char *argv[]);

#endif
