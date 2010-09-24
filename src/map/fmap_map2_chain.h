#ifndef FMAP_MAP2_CHAIN_H_
#define FMAP_MAP2_CHAIN_H_

#include "fmap_map2_aux.h"

// TODO: document
typedef struct {
    uint32_t tbeg, tend;
    int qbeg, qend;
    uint32_t flag:1, idx:31;
    int chain; // also reuse as a counter
} fmap_map2_chain_t;

// TODO: document
void 
fmap_map2_chain_filter(const fmap_map2_opt_t *opt, int len, fmap_map2_aln_t *b[2]);

#endif
