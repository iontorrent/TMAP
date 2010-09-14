#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_definitions.h"
#include "fmap_fq.h"

inline fmap_fq_t *
fmap_fq_init()
{
  fmap_fq_t *s = fmap_calloc(1, sizeof(fmap_fq_t), "s");
  s->name.s = s->comment.s = s->seq.s = s->qual.s = NULL;
  s->is_int = 0;
  return s;
}

inline void
fmap_fq_destroy(fmap_fq_t *seq)
{
  free(seq->name.s); free(seq->comment.s); free(seq->seq.s); free(seq->qual.s);
  free(seq);
}

inline fmap_fq_t*
fmap_fq_clone(fmap_fq_t *seq)
{
  fmap_fq_t *ret = fmap_calloc(1, sizeof(fmap_fq_t), "ret");

  fmap_string_copy(&ret->name, &seq->name);
  fmap_string_copy(&ret->comment, &seq->comment);
  fmap_string_copy(&ret->seq, &seq->seq);
  fmap_string_copy(&ret->qual, &seq->qual);
  ret->is_int = seq->is_int;

  return ret;
}

void 
fmap_fq_reverse(fmap_fq_t *seq, int32_t rev_comp)
{
  int i;

  fmap_string_reverse(&seq->seq);
  fmap_string_reverse(&seq->qual);

  if(1 == rev_comp) {
      if(1 == seq->is_int) { // bases are integer values
          for(i=0;i<seq->seq.l;i++) {
              seq->seq.s[i] = (4 <= seq->seq.s[i]) ? seq->seq.s[i] : 3 - seq->seq.s[i];
          }
      }
      else { // bases are ASCII values
          for(i=0;i<seq->seq.l;i++) {
              seq->seq.s[i] = nt_char_to_rc_char[(int)seq->seq.s[i]];
          }
      }
  }
}

void
fmap_fq_to_int(fmap_fq_t *seq)
{
  int i;
  if(1 == seq->is_int) return;
  for(i=0;i<seq->seq.l;i++) {
      seq->seq.s[i] = nt_char_to_int[(int)seq->seq.s[i]];
  }
  seq->is_int = 1;
}

void
fmap_fq_to_char(fmap_fq_t *seq)
{
  int i;
  if(0 == seq->is_int) return;
  for(i=0;i<seq->seq.l;i++) {
      seq->seq.s[i] = "ACGTN"[(int)seq->seq.s[i]];
  }
  seq->is_int = 0;
}
