#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_definitions.h"
#include "fmap_seq.h"

inline fmap_seq_t *
fmap_seq_init()
{
  fmap_seq_t *s = fmap_calloc(1, sizeof(fmap_seq_t), "s");
  s->name.s = s->comment.s = s->seq.s = s->qual.s = NULL;
  s->is_int = 0;
  return s;
}

inline void
fmap_seq_destroy(fmap_seq_t *seq)
{
  free(seq->name.s); free(seq->comment.s); free(seq->seq.s); free(seq->qual.s);
  free(seq);
}

static inline void
fmap_string_copy(fmap_string_t *dest, fmap_string_t *src)
{
  int32_t i;
  dest->s = fmap_malloc(sizeof(char)*src->m, "dest->s");
  for(i=0;i<src->m;i++) { // copy over all memory
      dest->s[i] = src->s[i];
  }
  dest->l = src->l;
  dest->m = src->m;
}

inline fmap_seq_t*
fmap_seq_clone(fmap_seq_t *seq)
{
  fmap_seq_t *ret = fmap_calloc(1, sizeof(fmap_seq_t), "ret");

  fmap_string_copy(&ret->name, &seq->name);
  fmap_string_copy(&ret->comment, &seq->comment);
  fmap_string_copy(&ret->seq, &seq->seq);
  fmap_string_copy(&ret->qual, &seq->qual);
  ret->is_int = seq->is_int;

  return ret;
}

static inline void
fmap_string_reverse(fmap_string_t *str)
{
  int i;
  for(i = 0; i < (str->l >> 1); ++i) {
      char tmp = str->s[str->l-1-i];
      str->s[str->l-1-i] = str->s[i]; str->s[i] = tmp;
  }
}

void 
fmap_seq_reverse(fmap_seq_t *seq, int32_t rev_comp)
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
fmap_seq_to_int(fmap_seq_t *seq)
{
  int i;
  if(1 == seq->is_int) return;
  for(i=0;i<seq->seq.l;i++) {
      seq->seq.s[i] = nt_char_to_int[(int)seq->seq.s[i]];
  }
  seq->is_int = 1;
}

void
fmap_seq_to_char(fmap_seq_t *seq)
{
  int i;
  if(0 == seq->is_int) return;
  for(i=0;i<seq->seq.l;i++) {
      seq->seq.s[i] = "ACGTN"[(int)seq->seq.s[i]];
  }
  seq->is_int = 0;
}
