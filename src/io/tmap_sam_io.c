/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <config.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../seq/tmap_sam.h"
#include "tmap_file.h"
#include "tmap_sam_io.h"
#include "tmap_seq_io.h"

#ifdef HAVE_SAMTOOLS
#include <sam_header.h>
// from sam.c
#define TYPE_BAM  1
#define TYPE_READ 2

extern char *bam_nt16_rev_table;

static inline tmap_sam_io_t *
tmap_sam_io_init_helper(const char *fn, int32_t is_bam)
{
  static char *rg_tags[] = {"ID", "CN", "DS", "DT", "FO", "KS", "LB", "PG", "PI", "PL", "PU", "SM"};
  tmap_sam_io_t *samio = NULL;
  char **p = NULL;
  int32_t i, n;

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

  for(i=0;i<12;i++) {
      p = sam_header2list(samio->fp->header->dict, "RG", rg_tags[i], &n);
      if(1 == n) samio->rg_tags[i] = p[0];
      else samio->rg_tags[i] = NULL;
      free(p); p = NULL;
  }

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
  samclose(samio->fp);
  free(samio);
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
      sam->name = tmap_string_init(len+1);
      memcpy(sam->name->s, str, len+1);
      sam->name->l = len;
      // seq and qual
      len = sam->b->core.l_qseq;
      sam->seq = tmap_string_init(len+1);
      sam->qual = tmap_string_init(len+1);
      memcpy(sam->qual->s, (char*)bam1_qual(sam->b), len);
      for(i=0;i<len;i++) {
          sam->seq->s[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(sam->b), i)];
          sam->qual->s[i] = QUAL2CHAR(sam->qual->s[i]);
      }
      sam->seq->l = sam->qual->l = len;
      sam->seq->s[len] = sam->qual->s[len] = '\0';
      // reverse compliment if necessary
      if((sam->b->core.flag & BAM_FREVERSE)) {
          tmap_sam_reverse_compliment(sam);
      }
      // save for later
      sam->io = samio;
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
#endif
