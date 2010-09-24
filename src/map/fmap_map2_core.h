#ifndef FMAP_MAP2_CORE_H_
#define FMAP_MAP2_CORE_H_

#include "../index/fmap_bwtl.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_sa.h"
#include "fmap_map2.h"
#include "fmap_map2_aux.h"
#include "fmap_map2_mempool.h"

#define FMAP_MAP2_MINUS_INF -0x3fffffff

// TODO: document
fmap_map2_aln_t **
fmap_map2_core_aln(const fmap_map2_opt_t *opt, const fmap_bwtl_t *target, 
               const fmap_bwt_t *query_bwt, const fmap_sa_t *query_sa,
               fmap_map2_global_mempool_t *pool);

#endif
