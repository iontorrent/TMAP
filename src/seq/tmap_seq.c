/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "tmap_fq.h"
#include "tmap_sff.h"
#include "tmap_seq.h"

tmap_seq_t *
tmap_seq_init(int8_t type)
{
  tmap_seq_t *seq = NULL;

  seq = tmap_calloc(1, sizeof(tmap_seq_t), "seq");
  seq->type = type;

  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      seq->data.fq = tmap_fq_init();
      break;
    case TMAP_SEQ_TYPE_SFF:
      seq->data.sff = tmap_sff_init();
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }

  return seq;
}

void
tmap_seq_destroy(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      tmap_fq_destroy(seq->data.fq);
      break;
    case TMAP_SEQ_TYPE_SFF:
      tmap_sff_destroy(seq->data.sff);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  free(seq);
}

tmap_seq_t *
tmap_seq_clone(tmap_seq_t *seq)
{
  tmap_seq_t *ret = NULL;

  ret = tmap_calloc(1, sizeof(tmap_seq_t), "ret");
  ret->type = seq->type;

  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      ret->data.fq = tmap_fq_clone(seq->data.fq);
      break;
    case TMAP_SEQ_TYPE_SFF:
      ret->data.sff = tmap_sff_clone(seq->data.sff);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }

  return ret;
}

void
tmap_seq_reverse(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      tmap_fq_reverse(seq->data.fq);
      break;
    case TMAP_SEQ_TYPE_SFF:
      tmap_sff_reverse(seq->data.sff);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
}

void
tmap_seq_reverse_compliment(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      tmap_fq_reverse_compliment(seq->data.fq);
      break;
    case TMAP_SEQ_TYPE_SFF:
      tmap_sff_reverse_compliment(seq->data.sff);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
}

void
tmap_seq_compliment(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      tmap_fq_compliment(seq->data.fq);
      break;
    case TMAP_SEQ_TYPE_SFF:
      tmap_sff_compliment(seq->data.sff);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
}

void
tmap_seq_to_int(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      tmap_fq_to_int(seq->data.fq);
      break;
    case TMAP_SEQ_TYPE_SFF:
      tmap_sff_to_int(seq->data.sff);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
}

tmap_string_t *
tmap_seq_get_name(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      return seq->data.fq->name;
      break;
    case TMAP_SEQ_TYPE_SFF:
      return seq->data.sff->rheader->name;
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return NULL;
}

inline tmap_string_t *
tmap_seq_get_bases(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      return tmap_fq_get_bases(seq->data.fq);
      break;
    case TMAP_SEQ_TYPE_SFF:
      return tmap_sff_get_bases(seq->data.sff);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return NULL;
}

inline tmap_string_t *
tmap_seq_get_qualities(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      return tmap_fq_get_qualities(seq->data.fq);
      break;
    case TMAP_SEQ_TYPE_SFF:
      return tmap_sff_get_qualities(seq->data.sff);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return NULL;
}

inline void
tmap_seq_remove_key_sequence(tmap_seq_t *seq)
{
  if(TMAP_SEQ_TYPE_SFF != seq->type) return; // ignore
  tmap_sff_remove_key_sequence(seq->data.sff);
}

tmap_seq_t *
tmap_seq_sff2fq(tmap_seq_t *seq)
{
  int32_t i;
  tmap_seq_t *ret= NULL;
  
  if(seq->type == TMAP_SEQ_TYPE_FQ) return tmap_seq_clone(seq);

  //Note:  ignore the comment field
  ret = tmap_seq_init(TMAP_SEQ_TYPE_FQ);
  tmap_string_copy(ret->data.fq->name, seq->data.sff->rheader->name); // name
  tmap_string_copy(ret->data.fq->seq, seq->data.sff->read->bases); // seq
  tmap_string_copy(ret->data.fq->qual, seq->data.sff->read->quality); // qual
  ret->data.fq->is_int = seq->data.sff->is_int; // is in integer format

  // remove key sequence
  for(i=0;i<ret->data.fq->seq->l - seq->data.sff->gheader->key_length;i++) {
      ret->data.fq->seq->s[i] = ret->data.fq->seq->s[i + seq->data.sff->gheader->key_length];
      ret->data.fq->qual->s[i] = ret->data.fq->qual->s[i + seq->data.sff->gheader->key_length];
  }
  ret->data.fq->seq->l -= seq->data.sff->gheader->key_length;
  ret->data.fq->qual->l -= seq->data.sff->gheader->key_length;

  return ret;
}
