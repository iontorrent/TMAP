/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <config.h>

#include "../samtools/kstring.h"
#include "../samtools/sam.h"
#include "../samtools/bam.h"
#include "../samtools/sam.h"
#include "../samtools/sam_header.h"

#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_string.h"
#include "tmap_sam.h"
#include "../io/tmap_sam_io.h"

tmap_sam_t *
tmap_sam_init()
{
  return tmap_calloc(1, sizeof(tmap_sam_t), "sam");
}

void
tmap_sam_destroy(tmap_sam_t *sam)
{
  if(NULL == sam) return;
  if(NULL != sam->name) tmap_string_destroy(sam->name);
  if(NULL != sam->seq) tmap_string_destroy(sam->seq);
  if(NULL != sam->qual) tmap_string_destroy(sam->qual);
  if(NULL != sam->b) bam_destroy1(sam->b);
  if(NULL != sam->flowgram) free(sam->flowgram);
  free(sam);
}

inline tmap_sam_t*
tmap_sam_clone(tmap_sam_t *sam)
{
  tmap_sam_t *ret = tmap_calloc(1, sizeof(tmap_sam_t), "ret");

  ret->name = tmap_string_clone(sam->name);
  ret->seq = tmap_string_clone(sam->seq);
  ret->qual = tmap_string_clone(sam->qual);
  ret->is_int = sam->is_int;

  // do not clone flow space info

  return ret;
}

void
tmap_sam_reverse(tmap_sam_t *sam)
{
  tmap_string_reverse(sam->seq);
  tmap_string_reverse(sam->qual);
}

void
tmap_sam_reverse_compliment(tmap_sam_t *sam)
{
  tmap_string_reverse_compliment(sam->seq, sam->is_int);
  tmap_string_reverse(sam->qual);
}

void
tmap_sam_compliment(tmap_sam_t *sam)
{
  tmap_string_compliment(sam->seq, sam->is_int);
}

void
tmap_sam_to_int(tmap_sam_t *sam)
{
  int i;
  if(1 == sam->is_int) return;
  for(i=0;i<sam->seq->l;i++) {
      sam->seq->s[i] = tmap_nt_char_to_int[(int)sam->seq->s[i]];
  }
  sam->is_int = 1;
}

void
tmap_sam_to_char(tmap_sam_t *sam)
{
  int i;
  if(0 == sam->is_int) return;
  for(i=0;i<sam->seq->l;i++) {
      sam->seq->s[i] = "ACGTN"[(int)sam->seq->s[i]];
  }
  sam->is_int = 0;
}

inline tmap_string_t *
tmap_sam_get_bases(tmap_sam_t *sam)
{
    return sam->seq;
}

inline tmap_string_t *
tmap_sam_get_qualities(tmap_sam_t *sam)
{
    return sam->qual;
}

void
tmap_sam_update_flow_info(tmap_sam_t *sam, tmap_sam_io_t *samio)
{
  uint8_t *tag;
  // init
  sam->fo = NULL;
  sam->ks = NULL;
  sam->rg_id = NULL;
  if(NULL != sam->flowgram) free(sam->flowgram);
  sam->flowgram = NULL;
  sam->flowgram_len = 0;
  // flowgram
  tag = bam_aux_get(sam->b, "FZ");
  if(NULL != tag) {
      sam->flowgram = bam_auxB2S(tag, &sam->flowgram_len);
  }
  // ZF
  tag = bam_aux_get(sam->b, "ZF");
  if(NULL != tag) sam->flow_start_index = bam_aux2i(tag);
  else sam->flow_start_index = -1;
  // RG
  tag = bam_aux_get(sam->b, "RG");
  if(NULL == tag) return;
  sam->rg_id = bam_aux2Z(tag);
  if(NULL == sam->rg_id) return;
  //sam->fo = sam_tbl_get(samio->rg_tbls[TMAP_SAM_RG_FO], (const char*)sam->rg_id);
  //sam->ks = sam_tbl_get(samio->rg_tbls[TMAP_SAM_RG_KS], (const char*)sam->rg_id);
}

int32_t
tmap_sam_get_flow_order_int(tmap_sam_t *sam, uint8_t **flow_order)
{
  int32_t i, flow_order_len;
  const char *fo = sam->fo;
  if(NULL == fo) return 0;
  flow_order_len = strlen(fo); // TODO: cache this
  (*flow_order) = tmap_malloc(sizeof(uint8_t) * flow_order_len, "flow_order");
  for(i=0;i<flow_order_len;i++) {
      (*flow_order)[i] = tmap_nt_char_to_int[(int)fo[i]];
  }
  return flow_order_len;
}

int32_t
tmap_sam_get_key_seq_int(tmap_sam_t *sam, uint8_t **key_seq)
{
  int32_t i, key_seq_len;
  const char *ks = sam->ks;
  if(NULL == ks) return 0;
  key_seq_len = strlen(ks); // TODO: cache this
  (*key_seq) = tmap_malloc(sizeof(uint8_t) * key_seq_len, "key_seq");
  for(i=0;i<key_seq_len;i++) {
      (*key_seq)[i] = tmap_nt_char_to_int[(int)ks[i]];
  }
  return key_seq_len;
}

int32_t
tmap_sam_get_flowgram(tmap_sam_t *sam, uint16_t **flowgram, int32_t mem)
{
  int32_t flowgram_len = sam->flowgram_len;
  if(NULL == sam->flowgram) return 0;
  if(mem < flowgram_len) (*flowgram) = tmap_realloc((*flowgram), sizeof(uint16_t) * flowgram_len, "(*flowgram)");
  memcpy((*flowgram), sam->flowgram, flowgram_len * sizeof(uint16_t));
  return flowgram_len;
}

int32_t
tmap_sam_get_flow_start_index(tmap_sam_t *sam)
{
  return sam->flow_start_index;
}

char*
tmap_sam_get_rg_id(tmap_sam_t *sam)
{
  return (char*)sam->rg_id;
}
