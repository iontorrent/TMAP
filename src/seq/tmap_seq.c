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
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      seq->data.sam = tmap_sam_init();
      break;
#endif
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
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_destroy(seq->data.sam);
      break;
#endif
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
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      ret->data.sam = tmap_sam_clone(seq->data.sam);
      break;
#endif
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
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_reverse(seq->data.sam);
      break;
#endif
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
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_reverse_compliment(seq->data.sam);
      break;
#endif
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
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_compliment(seq->data.sam);
      break;
#endif
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
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_to_int(seq->data.sam);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
}

void
tmap_seq_to_char(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      tmap_fq_to_char(seq->data.fq);
      break;
    case TMAP_SEQ_TYPE_SFF:
      tmap_sff_to_char(seq->data.sff);
      break;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_to_char(seq->data.sam);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
}

int32_t
tmap_seq_is_int(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      return seq->data.fq->is_int;
      break;
    case TMAP_SEQ_TYPE_SFF:
      return seq->data.sff->is_int;
      break;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return seq->data.sam->is_int;
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return -1;
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
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return seq->data.sam->name;
      break;
#endif
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
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return tmap_sam_get_bases(seq->data.sam);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return NULL;
}

inline int32_t
tmap_seq_get_bases_length(tmap_seq_t *seq)
{
  return tmap_seq_get_bases(seq)->l;
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
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return tmap_sam_get_qualities(seq->data.sam);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return NULL;
}

inline int32_t
tmap_seq_remove_key_sequence(tmap_seq_t *seq, int32_t remove_clipping, uint8_t *key_seq, int32_t key_seq_len)
{
  if(TMAP_SEQ_TYPE_SFF != seq->type) return 1; // ignore
  return tmap_sff_remove_key_sequence(seq->data.sff, remove_clipping, key_seq, key_seq_len);
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
  for(i=0;i<(int32_t)ret->data.fq->seq->l - seq->data.sff->gheader->key_length;i++) {
      ret->data.fq->seq->s[i] = ret->data.fq->seq->s[i + seq->data.sff->gheader->key_length];
      ret->data.fq->qual->s[i] = ret->data.fq->qual->s[i + seq->data.sff->gheader->key_length];
  }
  ret->data.fq->seq->l -= seq->data.sff->gheader->key_length;
  ret->data.fq->qual->l -= seq->data.sff->gheader->key_length;

  return ret;
}

int32_t
tmap_seq_get_flow_order_int(tmap_seq_t *seq, uint8_t **flow_order)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      break;
    case TMAP_SEQ_TYPE_SFF:
      return tmap_sff_get_flow_order_int(seq->data.sff, flow_order);
      break;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return tmap_sam_get_flow_order_int(seq->data.sam, flow_order);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return 0;
}

int32_t
tmap_seq_get_key_seq_int(tmap_seq_t *seq, uint8_t **key_seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      break;
    case TMAP_SEQ_TYPE_SFF:
      return tmap_sff_get_key_seq_int(seq->data.sff, key_seq);
      break;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return tmap_sam_get_key_seq_int(seq->data.sam, key_seq);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return 0;
}

// NB: includes key bases if present
int32_t
tmap_seq_get_flowgram(tmap_seq_t *seq, uint16_t **flowgram, int32_t mem)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      break;
    case TMAP_SEQ_TYPE_SFF:
      return tmap_sff_get_flowgram(seq->data.sff, flowgram, mem);
      break;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return tmap_sam_get_flowgram(seq->data.sam, flowgram, mem);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return 0;
}

int32_t
tmap_seq_get_flow_start_index(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      return -1;
    case TMAP_SEQ_TYPE_SFF:
      return tmap_sff_get_flow_start_index(seq->data.sff);
      break;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return tmap_sam_get_flow_start_index(seq->data.sam);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return -1;
}

char*
tmap_seq_get_rg_id(tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
    case TMAP_SEQ_TYPE_SFF:
      return NULL;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return tmap_sam_get_rg_id(seq->data.sam);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return NULL;
}
