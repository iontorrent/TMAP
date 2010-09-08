#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_seq.h"

inline fmap_seq_t *
fmap_seq_init()
{
  fmap_seq_t *s = fmap_calloc(1, sizeof(fmap_seq_t), "s");
  s->name.s = s->comment.s = s->seq.s = s->qual.s = NULL;
  return s;
}

inline void
fmap_seq_destroy(fmap_seq_t *seq)
{
  free(seq->name.s); free(seq->comment.s); free(seq->seq.s); free(seq->qual.s);
  free(seq);
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

  if(rev_comp) {
      for(i = 0;i < seq->seq.l; ++i) {
          seq->seq.s[i] = (4 <= seq->seq.s[i]) ? seq->seq.s[i] : 3 - seq->seq.s[i];
      }
  }
}
