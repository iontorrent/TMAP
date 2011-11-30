/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */

#include <stdlib.h>
#include <stdio.h>
#include <config.h>
#include <unistd.h>
#include "../../util/tmap_error.h"
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_definitions.h"
#include "tmap_map_stats.h"

tmap_map_stats_t*
tmap_map_stats_init()
{
  return tmap_calloc(1, sizeof(tmap_map_stats_t), "return");
}

void
tmap_map_stats_destroy(tmap_map_stats_t *s)
{
  free(s);
}

void
tmap_map_stats_add(tmap_map_stats_t *dest, tmap_map_stats_t *src)
{
  dest->num_reads += src->num_reads;
  dest->num_with_mapping += src->num_with_mapping;
  dest->num_after_seeding += src->num_after_seeding;
  dest->num_after_scoring += src->num_after_scoring;
  dest->num_after_rmdup += src->num_after_rmdup;
  dest->num_after_filter += src->num_after_filter;
}

void
tmap_map_stats_print(tmap_map_stats_t *s)
{
  fprintf(stderr, "num_reads=%lu\n", s->num_reads);
  fprintf(stderr, "num_with_mapping=%lu\n", s->num_with_mapping);
  fprintf(stderr, "num_after_seeding=%lu\n", s->num_after_seeding);
  fprintf(stderr, "num_after_scoring=%lu\n", s->num_after_scoring);
  fprintf(stderr, "num_after_rmdup=%lu\n", s->num_after_rmdup);
  fprintf(stderr, "num_after_filter=%lu\n", s->num_after_filter);
}
