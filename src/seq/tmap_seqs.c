/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "tmap_fq.h"
#include "tmap_sff.h"
#include "tmap_seq.h"
#include "tmap_seqs.h"

tmap_seqs_t *
tmap_seqs_init(int8_t type)
{
  tmap_seqs_t *seqs = NULL;

  seqs = tmap_calloc(1, sizeof(tmap_seqs_t), "seqs");
  seqs->type = type;
  seqs->seqs = NULL;
  seqs->n = seqs->m = 0;

  return seqs;
}

void
tmap_seqs_destroy(tmap_seqs_t *seqs)
{
  int32_t i;
  for(i=0;i<seqs->m;i++) {
      tmap_seq_destroy(seqs->seqs[i]);
  }
  free(seqs->seqs);
  free(seqs);
}

tmap_seqs_t *
tmap_seqs_clone(tmap_seqs_t *seqs)
{
  tmap_seqs_t *ret = NULL;
  int32_t i;

  ret = tmap_calloc(1, sizeof(tmap_seqs_t), "ret");
  ret->type = seqs->type;
  ret->n = seqs->n;
  ret->m = seqs->n; // do not expand memory

  if(0 < seqs->n) {
      ret->seqs = tmap_malloc(seqs->n * sizeof(tmap_seq_t*), "ret->seqs");
      for(i=0;i<ret->n;i++) {
          ret->seqs[i] = tmap_seq_clone(seqs->seqs[i]);
      }
  }

  return ret;
}

void
tmap_seqs_add(tmap_seqs_t *seqs, tmap_seq_t *seq)
{
  // do we need more memory?
  if(seqs->m <= seqs->n) {
      seqs->m++;
      seqs->seqs = tmap_realloc(seqs->seqs, seqs->m * sizeof(tmap_seq_t*), "seqs->seqs");
  }
  seqs->n++;
  seqs->seqs[seqs->n-1] = seq;
}

tmap_seq_t *
tmap_seqs_get(tmap_seqs_t *seqs, int32_t i)
{
  if(seqs->m <= i) { // make room
      seqs->seqs = tmap_realloc(seqs->seqs, (i+1) * sizeof(tmap_seq_t*), "seqs->seqs");
      while(seqs->m <= i) {
          seqs->seqs[seqs->m] = tmap_seq_init(seqs->type);
          seqs->m++;
      }
  }
  return seqs->seqs[i];
}
