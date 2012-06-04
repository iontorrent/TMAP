/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <config.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_definitions.h"
#include "../seq/tmap_sam.h"
#include "tmap_file.h"
#include "tmap_sam_io.h"
#include "tmap_seq_io.h"

#include "../samtools/sam_header.h"
// from sam.c
#define TYPE_BAM  1
#define TYPE_READ 2

extern char *bam_nt16_rev_table;

static inline tmap_sam_io_t *
tmap_sam_io_init_helper(const char *fn, int32_t is_bam)
{
  tmap_sam_io_t *samio = NULL;
  int32_t i;

  // initialize memory
  samio = tmap_calloc(1, sizeof(tmap_sam_io_t), "samio");
  if(0 == is_bam) {
      samio->fp = samopen(fn, "r", NULL);
  }
  else {
      samio->fp = samopen(fn, "rb", NULL);
  }
  if(NULL == samio->fp) {
      tmap_error(fn, Exit, OpenFileError);
  }

  // check if there are sequences in the header
  if(samio->fp->header->n_targets == 0) {
      tmap_error("Found no @SQ lines in the SAM header", Exit, OutOfRange);
  }

  /* HERE
  // create the rg tables
  samio->rg_tbls = tmap_calloc(TMAP_SAM_RG_NUM, sizeof(void*), "samio->rg_tabls");
  if(NULL == samio->fp->header->dict) samio->fp->header->dict = sam_header_parse2(samio->fp->header->text); // parse the header dictionary
  for(i=0;i<TMAP_SAM_RG_NUM;i++) {
      // table that goes from RG id to tag value
      samio->rg_tbls[i] = sam_header2tbl(samio->fp->header->dict, "RG", "ID", (char*)TMAP_SAM_RG_TAGS[i]);
  }
  // table of rg ids
  samio->rg_ids = sam_header2list(samio->fp->header->dict, "RG", (char*)TMAP_SAM_RG_TAGS[TMAP_SAM_RG_ID], &samio->rg_ids_num);
  */

  return samio;
}

inline tmap_sam_io_t *
tmap_sam_io_init(const char *fn)
{
  return tmap_sam_io_init_helper(fn, 0);
}

inline tmap_sam_io_t *
tmap_bam_io_init(const char *fn)
{
  return tmap_sam_io_init_helper(fn, 1);
}

void
tmap_sam_io_destroy(tmap_sam_io_t *samio)
{
  int32_t i;
  /* HERE
  for(i=0;i<TMAP_SAM_RG_NUM;i++) {
      sam_tbl_destroy(samio->rg_tbls[i]);
  }
  free(samio->rg_tbls);
  free(samio->rg_ids);
  */
  samclose(samio->fp);
  free(samio);
}

static void
tmap_sam_io_update_string(tmap_string_t **dst, char *src, int32_t len)
{
  if(NULL == (*dst)) (*dst) = tmap_string_init(len+1);
  else if((*dst)->m < len + 1) {
      tmap_string_destroy((*dst));
      (*dst) = tmap_string_init(len+1);
  }
  if(NULL != src) memcpy((*dst)->s, src, len);
  (*dst)->l = len;
}

int32_t
tmap_sam_io_read(tmap_sam_io_t *samio, tmap_sam_t *sam)
{
  if(NULL != sam->b) {
      bam_destroy1(sam->b);
  }
  sam->b = bam_init1();

  if(0 < samread(samio->fp, sam->b)) {
      char *str;
      int32_t i, len;

      // name
      str = bam1_qname(sam->b);
      len = strlen(str);
      tmap_sam_io_update_string(&sam->name, str, len);
      sam->name->s[len] = '\0';
      // seq and qual
      len = sam->b->core.l_qseq;
      tmap_sam_io_update_string(&sam->seq, NULL, len);
      tmap_sam_io_update_string(&sam->qual, (char*)bam1_qual(sam->b), len);
      for(i=0;i<len;i++) {
          sam->seq->s[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(sam->b), i)];
          sam->qual->s[i] = QUAL2CHAR(sam->qual->s[i]);
      }
      sam->seq->s[len] = sam->qual->s[len] = '\0';
      // reverse compliment if necessary
      if((sam->b->core.flag & BAM_FREVERSE)) {
          tmap_sam_reverse_compliment(sam);
      }
      // save for later
      tmap_sam_update_flow_info(sam, samio);
      return 1;
  }
  
  return -1;
}

int32_t
tmap_sam_io_read_buffer(tmap_sam_io_t *samio, tmap_sam_t **sam_buffer, int32_t buffer_length)
{
  int32_t n = 0;

  if(buffer_length <= 0) return 0;

  while(n < buffer_length && 0 <= tmap_sam_io_read(samio, sam_buffer[n])) {
      n++;
  }

  return n;
}

char***
tmap_sam_io_get_rg_header(tmap_sam_io_t *samio, int32_t *n)
{
  int32_t i, j;
  char ***header = NULL;
  char *rg_id = NULL;

  /* HERE
  if(0 == samio->rg_ids_num) {
      *n = 0;
      return NULL;
  }

  *n = samio->rg_ids_num;
  header = tmap_calloc(samio->rg_ids_num, sizeof(char**), "header");
  for(i=0;i<samio->rg_ids_num;i++) { // for each rg id
      header[i] = tmap_calloc(TMAP_SAM_RG_NUM, sizeof(char*), "header[i]");
      rg_id = samio->rg_ids[i];
      for(j=0;j<TMAP_SAM_RG_NUM;j++) { // for each tag
          header[i][j] = (char*)sam_tbl_get(samio->rg_tbls[j], (const char*)rg_id);
      }
  }
  */

  return header;
}
      
sam_header_t*
tmap_sam_io_get_sam_header(tmap_sam_io_t *samio)
{
  return samio->fp->header->header;
}
