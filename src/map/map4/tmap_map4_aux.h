/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP4_AUX_H
#define TMAP_MAP4_AUX_H

#include "tmap_map4.h"

// TODO
tmap_map4_aux_smem_iter_t *
tmap_map4_aux_smem_iter_init();

// TODO
void 
tmap_map4_aux_smem_iter_destroy(tmap_map4_aux_smem_iter_t *iter);

tmap_map_sams_t *
tmap_map4_aux_core(tmap_seq_t *seq,
                   tmap_refseq_t *refseq,
                   tmap_bwt_t *bwt,
                   tmap_sa_t *sa,
                   tmap_bwt_match_hash_t *hash,
                   tmap_map4_aux_smem_iter_t *iter,
                   tmap_map_opt_t *opt);

#endif
