#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "fmap_fq.h"

inline fmap_fq_t *
fmap_fq_init()
{
  fmap_fq_t *s = fmap_calloc(1, sizeof(fmap_fq_t), "s");
  s->name = fmap_string_init(0);
  s->comment = fmap_string_init(0);
  s->seq = fmap_string_init(0);
  s->qual= fmap_string_init(0);
  s->is_int = 0;

  return s;
}

inline void
fmap_fq_destroy(fmap_fq_t *fq)
{
  fmap_string_destroy(fq->name);
  fmap_string_destroy(fq->comment);
  fmap_string_destroy(fq->seq);
  fmap_string_destroy(fq->qual);
  free(fq);
}

inline fmap_fq_t*
fmap_fq_clone(fmap_fq_t *fq)
{
  fmap_fq_t *ret = fmap_calloc(1, sizeof(fmap_fq_t), "ret");

  ret->name = fmap_string_clone(fq->name);
  ret->comment = fmap_string_clone(fq->comment);
  ret->seq = fmap_string_clone(fq->seq);
  ret->qual = fmap_string_clone(fq->qual);
  ret->is_int = fq->is_int;

  return ret;
}

void
fmap_fq_reverse_compliment(fmap_fq_t *fq)
{
  fmap_string_reverse_compliment(fq->seq, fq->is_int);
  fmap_string_reverse(fq->qual);
}

void
fmap_fq_to_int(fmap_fq_t *fq)
{
  int i;
  if(1 == fq->is_int) return;
  for(i=0;i<fq->seq->l;i++) {
      fq->seq->s[i] = fmap_nt_char_to_int[(int)fq->seq->s[i]];
  }
  fq->is_int = 1;
}

void
fmap_fq_to_char(fmap_fq_t *fq)
{
  int i;
  if(0 == fq->is_int) return;
  for(i=0;i<fq->seq->l;i++) {
      fq->seq->s[i] = "ACGTN"[(int)fq->seq->s[i]];
  }
  fq->is_int = 0;
}

inline fmap_string_t *
fmap_fq_get_bases(fmap_fq_t *fq)
{
  return fq->seq;
}

inline fmap_string_t *
fmap_fq_get_qualities(fmap_fq_t *fq)
{
  return fq->qual;
}
