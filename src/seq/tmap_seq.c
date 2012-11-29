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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      seq->data.sam = tmap_sam_init();
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_destroy(seq->data.sam);
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      ret->data.sam = tmap_sam_clone(seq->data.sam);
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_reverse(seq->data.sam);
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_reverse_compliment(seq->data.sam);
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_compliment(seq->data.sam);
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_to_int(seq->data.sam);
      break;
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_to_char(seq->data.sam);
      break;
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return seq->data.sam->is_int;
      break;
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return seq->data.sam->name;
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return tmap_sam_get_bases(seq->data.sam);
      break;
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
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return tmap_sam_get_qualities(seq->data.sam);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return NULL;
}

inline int32_t
tmap_seq_remove_key_sequence(tmap_seq_t *seq, int32_t remove_clipping)
{
  if(TMAP_SEQ_TYPE_SFF != seq->type) return 1; // ignore
  return tmap_sff_remove_key_sequence(seq->data.sff, remove_clipping);
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

// NB: includes key bases if present
static int32_t
tmap_seq_get_flowgram(tmap_seq_t *seq, uint16_t **flowgram)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      break;
    case TMAP_SEQ_TYPE_SFF:
      return tmap_sff_get_flowgram(seq->data.sff, flowgram);
      break;
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      return tmap_sam_get_flowgram(seq->data.sam, flowgram);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return -1;
}

void
tmap_seq_update(tmap_seq_t *seq, int32_t idx, sam_header_t *header)
{
  char *rg_id = NULL;
  sam_header_records_t *records = NULL;
  sam_header_record_t **record_list = NULL;
  int32_t n = 0;

  // Read Group
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
    case TMAP_SEQ_TYPE_SFF:
      break;
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      rg_id = tmap_sam_get_rg_id(seq->data.sam);
      break;
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  if(NULL == rg_id) { // did not find in SAM/BAM
      // NB: assume that it is from the ith record in the header 
      records = sam_header_get_records(header, "RG");
      if(NULL != records) { // it exists
          if(idx < 0 || records->n <= idx) {
              tmap_error("RG records index was out of bounds", Exit, OutOfRange);
          }
          seq->rg_record = records->records[idx]; // copy over
          if(NULL == seq->rg_record) tmap_bug();
      }
  }
  else { // found in SAM/BAM
      n = 0;
      record_list = sam_header_get_record(header, "RG", "ID", rg_id, &n);
      if(0 == n) {
          fprintf(stderr, "Read Group Identifier: [%s]\n", rg_id);
          tmap_error("Did not find the @RG.ID in the SAM/BAM Header", Exit, OutOfRange);
      }
      else if(1 < n) {
          fprintf(stderr, "Read Group Identifier: [%s]\n", rg_id);
          tmap_error("Found more than one @RG.ID in the SAM/BAM Header", Exit, OutOfRange);
      }
      seq->rg_record = record_list[0];
      free(record_list); // NB: shallow copied
  }

  // Program Group
  // NB: assumes the last item in the header
  records = sam_header_get_records(header, "PG");
  if(NULL != records && 0 < records->n) { // it exists
      seq->pg_record = records->records[records->n-1]; // copy over
  }
  else {
      seq->pg_record = NULL;
  }

  // key sequence and flow order
  seq->fo_start_idx = -1;
  if(NULL != seq->rg_record) { // It should exist in the SAM/BAM Header
      seq->ks = sam_header_record_get(seq->rg_record, "KS");
      seq->fo = sam_header_record_get(seq->rg_record, "FO");
        
      // flow order index start
      if(NULL != seq->ks && NULL != seq->fo && TMAP_SEQ_TYPE_SFF == seq->type) { // only if it is an SFF
          // in addition, remove key sequence and trimming
          seq->fo_start_idx = tmap_seq_remove_key_sequence(seq, 1); 
      }
      else if(TMAP_SEQ_TYPE_SAM == seq->type || TMAP_SEQ_TYPE_BAM == seq->type) { // Try the ZF tag...
          seq->fo_start_idx = tmap_sam_get_fo_start_idx(seq->data.sam);
      }

      // flowgram information...
      seq->flowgram_len = tmap_seq_get_flowgram(seq, &seq->flowgram);

      // check if all flowspace information is available
      /*if((NULL == seq->ks || NULL == seq->fo || -1 == seq->fo_start_idx || NULL == seq->flowgram)// anything missing
         && (NULL != seq->ks || NULL != seq->fo || -1 != seq->fo_start_idx || NULL != seq->flowgram)) { // anything exists
          fprintf(stderr, "@RG.KS %s present.\n", (NULL == seq->ks) ? "is not" : "is");
          fprintf(stderr, "@RG.FO %s present.\n", (NULL == seq->fo) ? "is not" : "is");
          fprintf(stderr, "@SAM.FZ %s present.\n", (NULL == seq->flowgram) ? "is not" : "is");
          fprintf(stderr, "@SAM.ZF %s present.\n", (-1 == seq->fo_start_idx) ? "is not" : "is");
          tmap_error("Not all flowspace information available (@RG.KS and @RG.FO, and @SAM.FZ and @SAM.ZF)", Exit, OutOfRange);
      }*/
  }
}
