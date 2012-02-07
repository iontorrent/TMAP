/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP4_AUX_H
#define TMAP_MAP4_AUX_H

#include "tmap_map4.h"

/*!
  initializes the shared memory search iterator
  @return  the iterator
 */
tmap_map4_aux_smem_iter_t *
tmap_map4_aux_smem_iter_init();

/*!
  destroys the shared memory search iterator
  @param  iter  the iterator
 */
void 
tmap_map4_aux_smem_iter_destroy(tmap_map4_aux_smem_iter_t *iter);


/*!
  Core mapping routine
  @param  seq            the sequence to align (forward)
  @param  refseq         the reference sequence structure (forward)
  @param  bwt            the BWT structure 
  @param  sa             the SA structure 
  @param  hash           the occurrence hash
  @param  iter           the shared memory iterator
  @param  opt            the program options
  @return                the alignments
  the sequences should be in 2-bit format
  */
tmap_map_sams_t *
tmap_map4_aux_core(tmap_seq_t *seq,
                   tmap_refseq_t *refseq,
                   tmap_bwt_t *bwt,
                   tmap_sa_t *sa,
                   tmap_bwt_match_hash_t *hash,
                   tmap_map4_aux_smem_iter_t *iter,
                   tmap_rand_t *rand,
                   tmap_map_opt_t *opt);

#endif
