#ifndef FMAP_MAP3_AUX_H_
#define FMAP_MAP3_AUX_H_

#include "fmap_map3.h"

// TODO: document

typedef struct {
    uint32_t k;
    uint32_t l;
    uint32_t offset;
} fmap_map3_aux_seed_t;

typedef struct {
    uint32_t seqid; // index
    uint32_t pos; // shift
    uint32_t offset; // offset
} fmap_map3_aux_hit_t;

fmap_map3_aln_t *
fmap_map3_aln_init();

void
fmap_map3_aln_destroy(fmap_map3_aln_t *aln);

fmap_map3_aln_t *
fmap_map3_aux_core(fmap_seq_t *seq[2],
                   fmap_refseq_t *refseq,
                   fmap_bwt_t *bwt,
                   fmap_sa_t *sa,
                   fmap_map3_opt_t *opt);

#endif
