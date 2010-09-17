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

fmap_seq_t *
fmap_seq_clone(fmap_seq_t *seq)
{
  fmap_seq_t *ret = NULL;

  ret = fmap_calloc(1, sizeof(fmap_seq_t), "ret");
  ret->type = seq->type;

  switch(seq->type) {
    case FMAP_SEQ_TYPE_FQ:
      ret->data.fq = fmap_fq_clone(seq->data.fq);
      break;
    case FMAP_SEQ_TYPE_SFF:
      ret->data.sff = fmap_sff_clone(seq->data.sff);
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }

  return ret;
}

void
fmap_seq_reverse_compliment(fmap_seq_t *seq)
{
  switch(seq->type) {
    case FMAP_SEQ_TYPE_FQ:
      fmap_fq_reverse_compliment(seq->data.fq);
      break;
    case FMAP_SEQ_TYPE_SFF:
      fmap_sff_reverse_compliment(seq->data.sff);
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
}

void
fmap_seq_to_int(fmap_seq_t *seq)
{
  switch(seq->type) {
    case FMAP_SEQ_TYPE_FQ:
      fmap_fq_to_int(seq->data.fq);
      break;
    case FMAP_SEQ_TYPE_SFF:
      fmap_sff_to_int(seq->data.sff);
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
}

fmap_string_t *
fmap_seq_get_name(fmap_seq_t *seq)
{
  switch(seq->type) {
    case FMAP_SEQ_TYPE_FQ:
      return seq->data.fq->name;
      break;
    case FMAP_SEQ_TYPE_SFF:
      return seq->data.sff->rheader->name;
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return NULL;
}

fmap_string_t *
fmap_seq_get_bases(fmap_seq_t *seq)
{
  switch(seq->type) {
    case FMAP_SEQ_TYPE_FQ:
      return seq->data.fq->seq;
      break;
    case FMAP_SEQ_TYPE_SFF:
      return seq->data.sff->read->bases;
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return NULL;
}

fmap_string_t *
fmap_seq_get_qualities(fmap_seq_t *seq)
{
  switch(seq->type) {
    case FMAP_SEQ_TYPE_FQ:
      return seq->data.fq->qual;
      break;
    case FMAP_SEQ_TYPE_SFF:
      return seq->data.sff->read->quality;
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return NULL;
}
