/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP2_CORE_H_
#define TMAP_MAP2_CORE_H_

#include "../index/tmap_bwtl.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_sa.h"
#include "tmap_map2_mempool.h"
#include "tmap_map2_aux.h"
#include "tmap_map2.h"

#define TMAP_MAP2_MINUS_INF -0x3fffffff

/*! 
  */

/*! 
  the core alignment algorithm
  @param  opt         the program options
  @param  target      the target sequence (read)
  @param  query_bwt   the query bwt (reference)
  @param  query_sa    the query sa (reference)
  @param  pool        a global memory pool
  @return             a set of alignments
  */
tmap_map2_aln_t **
tmap_map2_core_aln(const tmap_map_opt_t *opt, const tmap_bwtl_t *target, 
               const tmap_bwt_t *query_bwt, const tmap_sa_t *query_sa,
               tmap_map2_global_mempool_t *pool);

#endif
