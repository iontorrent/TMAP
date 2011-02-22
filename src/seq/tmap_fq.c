/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "tmap_fq.h"

inline tmap_fq_t *
tmap_fq_init()
{
  tmap_fq_t *s = tmap_calloc(1, sizeof(tmap_fq_t), "s");
  s->name = tmap_string_init(0);
  s->comment = tmap_string_init(0);
  s->seq = tmap_string_init(0);
  s->qual= tmap_string_init(0);
  s->is_int = 0;

  return s;
}

inline void
tmap_fq_destroy(tmap_fq_t *fq)
{
  tmap_string_destroy(fq->name);
  tmap_string_destroy(fq->comment);
  tmap_string_destroy(fq->seq);
  tmap_string_destroy(fq->qual);
  free(fq);
}

inline tmap_fq_t*
tmap_fq_clone(tmap_fq_t *fq)
{
  tmap_fq_t *ret = tmap_calloc(1, sizeof(tmap_fq_t), "ret");

  ret->name = tmap_string_clone(fq->name);
  ret->comment = tmap_string_clone(fq->comment);
  ret->seq = tmap_string_clone(fq->seq);
  ret->qual = tmap_string_clone(fq->qual);
  ret->is_int = fq->is_int;

  return ret;
}

void
tmap_fq_reverse(tmap_fq_t *fq)
{
  tmap_string_reverse(fq->seq);
  tmap_string_reverse(fq->qual);
}

void
tmap_fq_reverse_compliment(tmap_fq_t *fq)
{
  tmap_string_reverse_compliment(fq->seq, fq->is_int);
  tmap_string_reverse(fq->qual);
}

void
tmap_fq_compliment(tmap_fq_t *fq)
{
  tmap_string_compliment(fq->seq, fq->is_int);
}

void
tmap_fq_to_int(tmap_fq_t *fq)
{
  int i;
  if(1 == fq->is_int) return;
  for(i=0;i<fq->seq->l;i++) {
      fq->seq->s[i] = tmap_nt_char_to_int[(int)fq->seq->s[i]];
  }
  fq->is_int = 1;
}

void
tmap_fq_to_char(tmap_fq_t *fq)
{
  int i;
  if(0 == fq->is_int) return;
  for(i=0;i<fq->seq->l;i++) {
      fq->seq->s[i] = "ACGTN"[(int)fq->seq->s[i]];
  }
  fq->is_int = 0;
}

inline tmap_string_t *
tmap_fq_get_bases(tmap_fq_t *fq)
{
  return fq->seq;
}

inline tmap_string_t *
tmap_fq_get_qualities(tmap_fq_t *fq)
{
  return fq->qual;
}
