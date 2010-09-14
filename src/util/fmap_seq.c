#include <stdlib.h>

#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_fq.h"
#include "fmap_sff.h"
#include "fmap_seq.h"

fmap_seq_t *
fmap_seq_init(int8_t type)
{
  fmap_seq_t *seq = NULL;

  seq = fmap_calloc(1, sizeof(fmap_seq_t), "seq");
  seq->type = type;

  switch(seq->type) {
    case FMAP_SEQ_TYPE_FQ:
      seq->data.fq = fmap_fq_init();
      break;
    case FMAP_SEQ_TYPE_SFF:
      seq->data.sff = fmap_sff_init();
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }

  return seq;
}

void
fmap_seq_destroy(fmap_seq_t *seq)
{
  switch(seq->type) {
    case FMAP_SEQ_TYPE_FQ:
      fmap_fq_destroy(seq->data.fq);
      break;
    case FMAP_SEQ_TYPE_SFF:
      fmap_sff_destroy(seq->data.sff);
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  free(seq);
}
