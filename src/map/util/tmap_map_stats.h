/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_STATS_H
#define TMAP_MAP_STATS_H

/*!
 * TODO
 */
typedef struct {
    uint64_t num_reads; /*!< the number of reads with at least one mapping */
    uint64_t num_with_mapping; /*!< the number of reads with at least one mapping */
    uint64_t num_after_seeding; /*!< the number of hits after seeding */
    uint64_t num_after_scoring; /*!< the number of hits after scoring */
    uint64_t num_after_rmdup; /*!< the number of hits after duplicate removal */
    uint64_t num_after_filter; /*!< the number of hits after filtering */
} tmap_map_stats_t;

/*!
  @return  a new stats structure
 */
tmap_map_stats_t*
tmap_map_stats_init();

/*!
  @param  s  the mapping driver stats to destroy
 */
void
tmap_map_stats_destroy(tmap_map_stats_t *s);

/*!
  Adds the src stats to the dest stats
  @param  dest  the destination
  @param  src   the source
 */
void
tmap_map_stats_add(tmap_map_stats_t *dest, tmap_map_stats_t *src);

#endif 
