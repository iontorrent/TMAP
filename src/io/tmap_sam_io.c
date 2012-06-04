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

static void
tmap_sam_io_init2_fs_and_add(tmap_seqs_io_t *io_in,
                         sam_header_records_t *records,
                         sam_header_record_t *record,
                         int32_t sam_flowspace_tags)
{
  char tag[2];
  // add @RG.KS and @RG.FO
  if(io_in->type == TMAP_SEQ_TYPE_SFF && sam_flowspace_tags) {
      if(io_in->n <= records->n) tmap_error("Too many read groups specified", Exit, OutOfRange);
      // @RG.KS
      tag[0]='K';tag[1]='S';
      if(0 == sam_header_record_add(record, tag, tmap_sff_io_get_rg_ks(io_in->seqios[records->n]->io.sff))) {
          tmap_error("Could not add the KS tag", Exit, OutOfRange);
      }
      // @RG.FO
      tag[0]='F';tag[1]='O';
      if(0 == sam_header_record_add(record, tag, tmap_sff_io_get_rg_fo(io_in->seqios[records->n]->io.sff))) {
          tmap_error("Could not add the FO tag", Exit, OutOfRange);
      }
  }
  // check for the @RG.ID and @RG.SM tags
  if(NULL == sam_header_record_get(record, "ID")) tmap_bug(); // should not happen
  if(NULL == sam_header_record_get(record, "SM")) tmap_error("SM tag missing from the read group", Exit, OutOfRange);
  // add the read group
  if(0 == sam_header_records_add(records, record)) tmap_bug(); 
}

inline tmap_sam_io_t *
tmap_sam_io_init2(const char *fn, const char *mod,
                  tmap_refseq_t *refseq,
                  tmap_seqs_io_t *io_in,
                  const char **rg_sam, int32_t rg_sam_num,
                  int32_t sam_flowspace_tags,
                  int32_t argc, const char *argv[])
{
  tmap_sam_io_t *io = NULL;
  sam_header_t *header = NULL; // the output header
  sam_header_records_t *records = NULL;
  sam_header_record_t *record = NULL;
  char tag[2];
  char *command_line= NULL;
  int32_t i;

  // @HD
  if(io_in->type == TMAP_SEQ_TYPE_SAM || io_in->type == TMAP_SEQ_TYPE_BAM) {
      // should be only one input file
      if(1 != io_in->n) {
          tmap_bug();
      }
      // get the current header
      header = io_in->seqios[0]->io.samio->fp->header->header; // wow, that's a lot of pointers
      if(NULL == header) {
          tmap_bug();
      }
  }
  else {
      // empty header
      header = sam_header_init();
      // @HD - header line
      records = sam_header_get_records(header, "HD"); // get the header line
      if(NULL == records) tmap_bug();
      record = sam_header_record_init("HD"); // new header line
      if(0 == sam_header_record_add(record, "VN", "1.4")) tmap_bug(); // version number
      if(0 == sam_header_records_add(records, record)) tmap_bug(); // add the header line
      // nullify
      record = NULL;
      records = NULL;
  }

  // @SQ
  if(NULL != refseq) {
      records = sam_header_get_records(header, "SQ"); // get the sequence dictionary
      if(NULL == records) tmap_bug();
      for(i=0;i<refseq->num_annos;i++) { // for each reference sequence
          char num[32];
          record = sam_header_record_init("SQ"); // new reference sequence record
          if(0 == sam_header_record_add(record, "SN", refseq->annos[i].name->s)) tmap_bug(); // reference sequence name
          if(sprintf(num, "%d", refseq->annos[i].len) < 0) tmap_bug(); // integer to string
          if(0 == sam_header_record_add(record, "SN", num)) tmap_bug(); // reference sequence name
          if(0 == sam_header_records_add(records, record)) tmap_bug(); // add the reference sequence record
      }
      record = NULL;
      records = NULL;
  }

  // @RG - read group
  if(0 < rg_sam_num) { // @RG specified on the command line
      // Check for SAM/BAM
      if(io_in->type == TMAP_SEQ_TYPE_SAM || io_in->type == TMAP_SEQ_TYPE_BAM) {
          tmap_error("Cannot specify the read groups on the command line when using SAM/BAM as input."
                     "  Please embed in the SAM/BAM header instead.", Exit, OutOfRange);
      }
      record = NULL;
      records = sam_header_get_records(header, "RG"); // get the read group
      if(NULL == records) tmap_bug();
      // go through the command line arguments
      for(i=0;i<sam_rg_num;i++) {
          if(strlen(sam_rg[i]) < 4) tmap_error("Read group too small", Exit, OutOfRange);
          if(':' != sam_rg[i][2]) tmap_error("Read group improperly formatted (no colon)", Exit, OutOfRange);

          // check for id
          if('I' == sam_rg[i][0] && 'D' == sam_rg[i][1]) { // new read group
              if(NULL != record) { // add the record
                  tmap_sam_io_init2_fs_and_add(io_in, records, record, sam_flowspace_tags); // add @RG.KS and @RG.FO
              }
              record = sam_header_record_init("RG"); // new read group
          }
          // add the tag/value to the record
          tag[0]=sam_rg[i][0]; tag[1]=sam_rg[i][1]; // setup the tag
          if(0 == sam_header_record_add(record, tag, sam_rg[i]+3)) tmap_bug(); // add the tag/value
      }
      if(NULL != record) { // add the record
          tmap_sam_io_init2_fs_and_add(io_in, records, record, sam_flowspace_tags); // add @RG.KS and @RG.FO
      }
      // check that the # of read groups added was the same as the # of input files...
      if(records->n != io_in->n) tmap_error("The number of read groups did not match the number of input files", Exit, OutOfRange);
  }
  else if(io_in->type != TMAP_SEQ_TYPE_SAM && io_in->type != TMAP_SEQ_TYPE_BAM) { // dummy...
      records = sam_header_get_records(header, "RG"); // get the read group
      if(NULL == records) tmap_bug();
      for(i=0;i<io_in->n;i++) { // for each input file
          char buf[32];
          record = sam_header_record_init("RG"); // new read group
          if(1 == io_in->n) strcpy(buf, "NOID");
          else if(sprintf(buf, "NOID.%d", i+1) < 0) tmap_bug();
          if(0 == sam_header_record_add(record, "ID", buf)) tmap_bug(); // dummy ID
          if(0 == sam_header_record_add(record, "SM", "NOSM")) tmap_bug(); // dummy SM, for Picard validation
          tmap_sam_io_init2_fs_and_add(io_in, records, record, sam_flowspace_tags); // add @RG.KS and @RG.FO
      }
  }

  // @PG - program group
  // TODO: check for previous program group ID and set @PG.PP
  records = sam_header_get_records(header, "PG"); // get the program group records
  if(NULL == records) tmap_bug();
  record = sam_header_record_init("PG"); // new program group
  if(0 == sam_header_record_add(record, "ID", PACKAGE_NAME)) tmap_bug(); // @PG.ID
  if(0 == sam_header_record_add(record, "VN", PACKAGE_VERSION)) tmap_bug(); // @PG.VN
  // @PG.CL
  command_line = NULL;
  j = 1; // for the EOL
  for(i=0;i<argc;i++) {
      if(0 < i) j++;
      j += strlen(argv[i]);
      command_line = tmap_realloc(command_line, sizeof(char) * j, "command_line");
      if(0 < i) strcat(command_line, " ");
      strcat(command_line, argv[i]);
  }
  if(0 == sam_header_record_add(record, "CL", command_line)) tmap_bug(); // @PG.CL
  free(command_line);

  // TODO: create bam_header_t

  // 2. Open the BAM file...

  return io;
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
