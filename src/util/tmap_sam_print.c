/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <config.h>

#ifdef HAVE_SAMTOOLS
#include <kstring.h>
#include <sam.h>
#include <bam.h>
#endif

#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_string.h"
#include "../io/tmap_file.h"
#include "../io/tmap_seq_io.h"
#include "../sw/tmap_sw.h"
#include "tmap_sam_print.h"

static char tmap_sam_rg_id[1024]="NOID";
static int32_t tmap_sam_rg_id_use = 0;

#define TMAP_SAM_NO_RG_SM "NOSM"

/*
static char **
tmap_sam_parse_rg(char *rg)
{
  char **header = NULL;
  int32_t i, j, len, tag_i;
  int32_t tags_found[TMAP_SAM_RG_NUM];
  int32_t num_found = 0; 

  if(NULL == rg) return NULL;

  header = tmap_calloc(TMAP_SAM_RG_NUM, sizeof(char*), "header");

  len = strlen(rg);

  // init
  for(i=0;i<TMAP_SAM_RG_NUM;i++) {
      tags_found[i] = 0;
  }

  // convert strings of "\t" to tab characters '\t'
  for(i=0;i<len-1;i++) {
      if(rg[i] == '\\' && rg[i+1] == 't') {
          rg[i] = '\t';
          // shift down
          for(j=i+1;j<len-1;j++) {
              rg[j] = rg[j+1];
          }
          len--;
          rg[len]='\0';
      }
  }

  // must have at least "@RG\t"
  if(len < 4
     || 0 != strncmp(rg, "@RG\t", 4)) {
      tmap_error("Malformed RG line", Exit, OutOfRange);
  }
  i = 3;

  while(i<len) {
      if('\t' == rg[i]) {
          i++; // move past the tab
          if(len <= i+2) { // must have "XX:" 
              tmap_error("Improper tag in the RG line", Exit, OutOfRange);
          }
          for(tag_i=0;tag_i<TMAP_SAM_RG_NUM;tag_i++) {
              if(TMAP_SAM_RG_TAGS[tag_i][0] == rg[i] && TMAP_SAM_RG_TAGS[tag_i][1] == rg[i+1]) { // found!
                  tags_found[tag_i]++;
                  num_found++;
                  // copy over
                  if(1 < tags_found[tag_i]) {
                      tmap_file_fprintf(tmap_file_stderr, "\nFound multiple %s tags for the RG SAM header.\n", TMAP_SAM_RG_TAGS[tag_i]);
                      tmap_error(NULL, Exit, OutOfRange);
                  }
                  else { // 1 == tags_found[tag_i]
                      int32_t l = 0;
                      // get length
                      for(j=i+3;j<len && '\t' != rg[j];j++) {
                          l++;
                      }
                      if(l <= 0) {
                          tmap_file_fprintf(tmap_file_stderr, "\nFound an empty tag in the RG SAM header: %s.\n", TMAP_SAM_RG_TAGS[tag_i]);
                          tmap_error(NULL, Exit, OutOfRange);
                      }
                      header[tag_i] = tmap_malloc(sizeof(char) * (1 + l), "header[tag_i]");
                      for(j=0;j<l;j++) {
                          header[tag_i][j] = rg[j + i + 3];
                      }
                      header[tag_i][l] = '\0';
                      i += l + 3 - 1; // subtract one since below we add one to i
                  }
                  break;
              }
          }
          if(TMAP_SAM_RG_NUM == tag_i) {
              tmap_error("Improper tag in the RG line", Exit, OutOfRange);
          }
      }
      i++;
  }
  if(0 == num_found) {
      free(header);
      return NULL;
  }
  return header;
}
*/

void
tmap_sam_print_header(tmap_file_t *fp, tmap_refseq_t *refseq, tmap_seq_io_t *seqio, char *sam_rg, 
                      int32_t sam_flowspace_tags, int32_t ignore_rg_sam_tags, 
                      int argc, char *argv[])
{
  int32_t i, j, header_n = 0;
  char **header_a = NULL;
  char ***header_b = NULL;

  // SAM header
  tmap_file_fprintf(fp, "@HD\tVN:%s\tSO:unsorted\n",
                    TMAP_SAM_PRINT_VERSION);
  if(NULL != refseq) {
      for(i=0;i<refseq->num_annos;i++) {
          tmap_file_fprintf(fp, "@SQ\tSN:%s\tLN:%d\n",
                            refseq->annos[i].name->s, (int)refseq->annos[i].len);
      }
  }
  // RG
  /* HERE
  header_a = tmap_sam_parse_rg(sam_rg); // parse the input read group line
  if(1 == ignore_rg_sam_tags) { // do not get the header from the input file
      if(1 == sam_flowspace_tags) { // ... except for the RG.FS/RG.KO
          // get the RG header from the input file
          header_b = tmap_seq_io_get_rg_header(seqio, &header_n);
          if(1 < header_n) { 
              // TODO: we could check to see that FO/KS are the same across all
              // input read groups
              tmap_error("Command line read group found with multiple read groups from the input file", Exit, OutOfRange);
          }
          else if(1 == header_n) {
              if(NULL == header_a) {
                  header_a = tmap_calloc(TMAP_SAM_RG_NUM, sizeof(char*), "header_a");
                  // copy over default RG.ID
                  header_a[TMAP_SAM_RG_ID] = tmap_malloc(sizeof(char) * (strlen(tmap_sam_rg_id) + 1), "header_a[TMAP_SAM_RG_ID]");
                  strcpy(header_a[TMAP_SAM_RG_ID], tmap_sam_rg_id);
              }
              for(i=0;i<TMAP_SAM_RG_NUM;i++) { // for each RG.TAG
                  switch(i) {
                    case TMAP_SAM_RG_FO:
                    case TMAP_SAM_RG_KS:
                      if(NULL != header_a[i] && NULL != header_b[0][i]) {
                          tmap_error("Command line and input read groups share tags", Exit, OutOfRange);
                      }
                      else if(NULL != header_b[0][i]) { // copy over
                          header_a[i] = tmap_malloc(sizeof(char) * (strlen(header_b[0][i]) + 1), "header_a[i]");
                          strcpy(header_a[i], header_b[0][i]);
                      }
                    default:
                      break;
                  }
              }
          }
          // free header_b, it is no longer in use
          for(i=0;i<header_n;i++) {
              free(header_b[i]);
          }
          free(header_b);
          header_b = NULL;
          header_n = 0;
      }
  }
  else { 
      // get the RG header from the input file
      header_b = tmap_seq_io_get_rg_header(seqio, &header_n);
  }

  // reconcile the RG headers
  if(NULL != header_a) { // header a exists
      if(NULL != header_b && 1 == header_n) { // header b exists, and only one line...
          // check to see if they are mutually exclusive
          for(i=0;i<TMAP_SAM_RG_NUM;i++) {
              if(NULL != header_a[i] && NULL != header_b[0][i]) {
                  tmap_file_fprintf(tmap_file_stderr, "\nFound both command line and input file read group information for the same tag: %s.\n", TMAP_SAM_RG_TAGS[i]);
                  tmap_error(NULL, Exit, OutOfRange);
              }
              else if(NULL == header_a[i] && NULL != header_b[0][i]) { // copy over
                  header_a[i] = tmap_calloc(1+strlen(header_b[0][i]), sizeof(char), "header_a[i]");
                  strcpy(header_a[i], header_b[0][i]);
              }
          }
          // free
          free(header_b[0]);
          free(header_b);
          header_b = NULL;
          header_n = 0;
      }
      if(0 == header_n) { // no header b
          if(NULL != header_a[TMAP_SAM_RG_ID]) {
              strcpy(tmap_sam_rg_id, header_a[TMAP_SAM_RG_ID]);
          }
          else {
              header_a[TMAP_SAM_RG_ID] = tmap_malloc(sizeof(char) * (strlen(tmap_sam_rg_id) + 1), "header_a[i]");
              strcpy(header_a[TMAP_SAM_RG_ID], tmap_sam_rg_id);
          }
          if(NULL == header_a[TMAP_SAM_RG_SM]) { // for Picard
              header_a[TMAP_SAM_RG_SM] = tmap_malloc(sizeof(char) * (strlen(TMAP_SAM_NO_RG_SM) + 1), "header_a[i]");
              strcpy(header_a[TMAP_SAM_RG_SM], TMAP_SAM_NO_RG_SM);
          }
          tmap_sam_rg_id_use = 1;
          tmap_file_fprintf(fp, "@RG");
          for(i=0;i<TMAP_SAM_RG_NUM;i++) {
              if(NULL != header_a[i]) {
                  tmap_file_fprintf(fp, "\t%s:%s", TMAP_SAM_RG_TAGS[i], header_a[i]);
              }
          }
          tmap_file_fprintf(fp, "\n");
      }
      else { // both header_a and header_b exist 
          tmap_error("Found both command line and input file read group information", Exit, OutOfRange);
      }
  }
  else { // no header_a exists
      if(NULL != header_b) { // no header_b exists
          tmap_sam_rg_id_use = 0;
          for(i=0;i<header_n;i++) { // for each RG.ID
              if(NULL == header_b[i][TMAP_SAM_RG_ID]) {
                  if(1 == header_n && TMAP_SEQ_TYPE_SFF == seqio->type) { // make an exception for SFF files
                      header_b[i][TMAP_SAM_RG_ID] = tmap_sam_rg_id;
                      tmap_sam_rg_id_use = 1;
                  }
                  else {
                      tmap_error("missing RG.ID found in the RG SAM Header", Exit, OutOfRange);
                  }
              }
              // RG.SM for picard
              if(NULL == header_b[i][TMAP_SAM_RG_SM]) {
                  header_b[i][TMAP_SAM_RG_SM] = TMAP_SAM_NO_RG_SM;
              }
              tmap_file_fprintf(fp, "@RG");
              for(j=0;j<TMAP_SAM_RG_NUM;j++) { // for each RG.TAG
                  if(NULL != header_b[i][j]) {
                      tmap_file_fprintf(fp, "\t%s:%s", TMAP_SAM_RG_TAGS[j], header_b[i][j]);
                  }
              }
              tmap_file_fprintf(fp, "\n");
          }
      }
      else {
          header_a = tmap_calloc(TMAP_SAM_RG_NUM, sizeof(char*), "header_a");
          // RG.ID
          header_a[TMAP_SAM_RG_ID] = tmap_malloc(sizeof(char) * (strlen(tmap_sam_rg_id) + 1), "header_a[i]");
          strcpy(header_a[TMAP_SAM_RG_ID], tmap_sam_rg_id);
          // RG.SM for Picard
          header_a[TMAP_SAM_RG_SM] = tmap_malloc(sizeof(char) * (strlen(TMAP_SAM_NO_RG_SM) + 1), "header_a[i]");
          strcpy(header_a[TMAP_SAM_RG_SM], TMAP_SAM_NO_RG_SM);
          tmap_sam_rg_id_use = 1;
          tmap_file_fprintf(fp, "@RG");
          for(i=0;i<TMAP_SAM_RG_NUM;i++) {
              if(NULL != header_a[i]) {
                  tmap_file_fprintf(fp, "\t%s:%s", TMAP_SAM_RG_TAGS[i], header_a[i]);
              }
          }
          tmap_file_fprintf(fp, "\n");
      }
  }
  */

  // PG
  tmap_file_fprintf(fp, "@PG\tID:%s\tVN:%s\tCL:",
                    PACKAGE_NAME, PACKAGE_VERSION);
  for(i=0;i<argc;i++) {
      if(0 < i) tmap_file_fprintf(fp, " ");
      tmap_file_fprintf(fp, "%s", argv[i]);
  }
  tmap_file_fprintf(fp, "\n");

  // free
  for(i=0;i<header_n;i++) {
      free(header_b[i]);
  }
  free(header_b);
  if(NULL != header_a) {
      for(i=0;i<TMAP_SAM_RG_NUM;i++) {
          free(header_a[i]);
      }
  }
  free(header_a);
}

static inline void
tmap_sam_print_flowgram(tmap_file_t *fp, uint16_t *flowgram, int32_t length)
{
  int32_t i;
  tmap_file_fprintf(fp, "\tFZ:B:S");
  for(i=0;i<length;i++) {
      tmap_file_fprintf(fp, ",%u", flowgram[i]);
  }
}

static inline void
tmap_sam_print_rg(tmap_file_t *fp, tmap_seq_t *seq)
{
  // RG 
  if(1 == tmap_sam_rg_id_use) {
      tmap_file_fprintf(fp, "\tRG:Z:%s", tmap_sam_rg_id);
  }
  else if(0 == tmap_sam_rg_id_use) {
      char *id = tmap_seq_get_rg_id(seq);
      if(NULL == id) {
          tmap_error("Missing Record RG.ID in the input file", Exit, OutOfRange);
      }
      tmap_file_fprintf(fp, "\tRG:Z:%s", id);
  }
}

static inline void 
tmap_sam_print_fz_and_zf(tmap_file_t *fp, tmap_seq_t *seq)
{
  uint16_t *flowgram = NULL;
  int32_t flow_start_index;
  int32_t flowgram_len;
  flowgram_len = tmap_seq_get_flowgram(seq, &flowgram, 0);
  if(NULL != flowgram) {
      tmap_sam_print_flowgram(fp, flowgram, flowgram_len);
      free(flowgram);
  }
  flow_start_index = tmap_seq_get_flow_start_index(seq);
  if(0 <= flow_start_index) {
      tmap_file_fprintf(fp, "\tZF:i:%d", flow_start_index);
  }
}

inline void
tmap_sam_print_unmapped(tmap_file_t *fp, tmap_seq_t *seq, int32_t sam_flowspace_tags, int32_t bidirectional, tmap_refseq_t *refseq,
                      uint32_t end_num, uint32_t m_unmapped, uint32_t m_prop, 
                      uint32_t m_strand, uint32_t m_seqid, uint32_t m_pos,
                      const char *format, ...)
{
  uint32_t flag = 0;
  tmap_string_t *name=NULL, *bases=NULL, *qualities=NULL;
  va_list ap;

  name = tmap_seq_get_name(seq);
  bases = tmap_seq_get_bases(seq);
  qualities = tmap_seq_get_qualities(seq);
  
  // set the flag
  flag = 0x4;
  if(0 < end_num) { // mate info
      flag |= 0x1;
      if(1 == m_prop) flag |= 0x2; // properly aligned
      if(1 == m_unmapped) flag |= 0x8; // unmapped
      else if(1 == m_strand) flag |= 0x20; // strand 
      flag |= (1 == end_num) ? 0x40 : 0x80; // first/second end
  }

  // name, flag, seqid, pos, mapq, cigar
  tmap_file_fprintf(fp, "%s\t%u\t*\t%u\t%u\t*", name->s, flag, 0, 0);

  // NB: hard clipped portions of the read is not reported
  // mate info
  if(0 == end_num) { // no mate
      tmap_file_fprintf(fp, "\t*\t0\t0");
  }
  else if(1 == m_unmapped) { // unmapped mate
      tmap_file_fprintf(fp, "\t*\t0\t0");
  }
  else if(NULL != refseq) { // mapped mate
      tmap_file_fprintf(fp, "\t%s\t%u\t%u",
                        refseq->annos[m_seqid].name->s,
                        m_pos+1,
                        0);
  }

  // bases and qual
  tmap_file_fprintf(fp, "\t%s\t%s",
                    (0 == bases->l) ? "*" : bases->s, (0 == qualities->l) ? "*" : qualities->s);

  // RG 
  tmap_sam_print_rg(fp, seq);

  // PG
  tmap_file_fprintf(fp, "\tPG:Z:%s", PACKAGE_NAME);

  // FZ and ZF
  if(1 == sam_flowspace_tags) {
      tmap_sam_print_fz_and_zf(fp, seq);
  }
  if(1 == bidirectional) {
      tmap_file_fprintf(fp, "\tXB:i:1");
  }
  
  // optional tags
  if(NULL != format) {
      va_start(ap, format);
      tmap_file_vfprintf(fp, format, ap);
      va_end(ap);
  }
  tmap_file_fprintf(fp, "\n");
}

static inline tmap_string_t *
tmap_sam_md(tmap_refseq_t *refseq, char *read_bases, // read bases are characters
            uint32_t seqid, uint32_t pos, // seqid and pos are 0-based
            uint32_t *cigar, int32_t n_cigar, int32_t *nm, char *read_bases_eq)
{
  int32_t i, j;
  uint32_t ref_i, read_i, ref_start, ref_end;
  int32_t l = 0; // the length of the last md op
  uint8_t read_base, ref_base;
  tmap_string_t *md=NULL;
  uint8_t *target = NULL;;

  md = tmap_string_init(32);
  (*nm) = 0;

  ref_start = ref_end = pos + 1; // make one-based
  for(i=0;i<n_cigar;i++) { // go through each cigar operator
      int32_t op_len;
      op_len = cigar[i] >> 4;
      switch(cigar[i]&0xf) {
        case BAM_CMATCH:
        case BAM_CDEL:
        case BAM_CREF_SKIP:
          ref_end += op_len; break;
        default:
          break;
      }
  }
  ref_end--;
      
  target = tmap_refseq_subseq2(refseq, seqid+1, ref_start, ref_end, NULL, 0, NULL);
  if(NULL == target) {
      tmap_bug();
  }

  if(0 == n_cigar) {
      tmap_bug();
  }

  read_i = ref_i = 0;
  for(i=0;i<n_cigar;i++) { // go through each cigar operator
      int32_t op_len, op;

      op_len = cigar[i] >> 4;
      op = cigar[i] & 0xf;

      if(BAM_CMATCH == op) {
          for(j=0;j<op_len;j++) {
              if(refseq->len <= refseq->annos[seqid].offset + pos + ref_i) break; // out of boundary

              read_base = tmap_nt_char_to_int[(int)read_bases[read_i]]; 
              ref_base = target[ref_i];

              if(read_base == ref_base) { // a match
                  if(NULL != read_bases_eq) read_bases_eq[read_i] = '=';
                  l++;
              }
              else {
                  if(NULL != read_bases_eq) read_bases_eq[read_i] = read_bases[read_i];
                  tmap_string_lsprintf(md, md->l, "%d%c", l, tmap_iupac_int_to_char[ref_base]);
                  l = 0;
                  (*nm)++;
              }
              read_i++;
              ref_i++; 
          }
          if(j < op_len) break;
      }
      else if(BAM_CINS == op) {
          if(NULL != read_bases_eq) {
              for(j=0;j<op_len;j++) {
                  read_bases_eq[read_i+j] = read_bases[read_i+j];
              }
          }
          read_i += op_len;
          (*nm) += op_len;
      }
      else if(BAM_CDEL == op) {
          tmap_string_lsprintf(md, md->l, "%d^", l);
          for(j=0;j<op_len;j++) {
              if(refseq->len <= refseq->annos[seqid].offset + pos + ref_i) break; // out of boundary
              ref_base = target[ref_i];
              tmap_string_lsprintf(md, md->l, "%c", tmap_iupac_int_to_char[ref_base]);
              ref_i++;
          }
          if(j < op_len) break;
          (*nm) += op_len;
          l=0;
      }
      else if(BAM_CREF_SKIP == op) {
          ref_i += op_len;
      }
      else if(BAM_CSOFT_CLIP == op) {
          if(NULL != read_bases_eq) {
              for(j=0;j<op_len;j++) {
                  read_bases_eq[read_i+j] = read_bases[read_i+j];
              }
          }
          read_i += op_len;
      }
      else if(BAM_CHARD_CLIP == op) {
          // ignore
      }
      else if(BAM_CPAD == op) {
          // ignore
      }
      else {
          tmap_error("could not understand the cigar operator", Exit, OutOfRange);
      }
  }
  tmap_string_lsprintf(md, md->l, "%d", l);
  if(NULL != read_bases_eq) read_bases_eq[read_i] = '\0';

  free(target);

  return md;
}

inline void
tmap_sam_print_mapped(tmap_file_t *fp, tmap_seq_t *seq, int32_t sam_flowspace_tags, int32_t bidirectional, int32_t seq_eq, tmap_refseq_t *refseq,
                      uint8_t strand, uint32_t seqid, uint32_t pos, int32_t aln_num,
                      uint32_t end_num, uint32_t m_unmapped, uint32_t m_prop, double m_num_std, uint32_t m_strand,
                      uint32_t m_seqid, uint32_t m_pos, uint32_t m_tlen,
                      uint8_t mapq, uint32_t *cigar, int32_t n_cigar,
                      int32_t score, int32_t ascore, int32_t pscore, int32_t nh, int32_t algo_id, int32_t algo_stage,
                      const char *format, ...)
{
  va_list ap;
  int32_t i;
  tmap_string_t *name=NULL, *bases=NULL, *qualities=NULL;
  char *bases_eq=NULL;
  uint32_t flag;
  tmap_string_t *md;
  int32_t nm;

  /*
  fprintf(stderr, "end_num=%d m_unmapped=%d m_prop=%d m_strand=%d m_seqid=%d m_pos=%d m_tlen=%d\n",
          end_num, m_unmapped, m_prop, m_strand, m_seqid, m_pos, m_tlen);
  */

  name = tmap_seq_get_name(seq);
  bases = tmap_seq_get_bases(seq);
  qualities = tmap_seq_get_qualities(seq);

  if(1 == strand) { // reverse for the output
      tmap_string_reverse_compliment(bases, 0);
      tmap_string_reverse(qualities);
  }

  if(0 == pos + 1) {
      tmap_error("position is out of range", Exit, OutOfRange);
  }

  // compute the MD/NM
  if(1 == seq_eq && 0 < bases->l) {
      bases_eq = tmap_calloc((1 + bases->l), sizeof(char), "bases_eq");
  }
  else {
      bases_eq = NULL;
  }
  md = tmap_sam_md(refseq, bases->s, seqid, pos, cigar, n_cigar, &nm, bases_eq);

  // flag
  flag = 0;
  if(1 == strand) flag |= 0x10; // strand
  if(0 < aln_num) flag |= 0x100; // secondary alignment
  if(0 < end_num) { // mate info
      flag |= 0x1;
      if(0 == m_unmapped && 1 == m_prop) flag |= 0x2; // properly aligned
      if(1 == m_unmapped) flag |= 0x8; // unmapped
      else if(1 == m_strand) flag |= 0x20; // strand 
      flag |= (1 == end_num) ? 0x40 : 0x80; // first/second end
  }

  tmap_file_fprintf(fp, "%s\t%u\t%s\t%u\t%u\t",
                    name->s, flag, refseq->annos[seqid].name->s,
                    pos + 1,
                    mapq);

  // print out the cigar
  if(TMAP_SEQ_TYPE_SFF == seq->type) {
      if(0 == strand && 0 < seq->data.sff->rheader->clip_left) {
          tmap_file_fprintf(fp, "%dH", seq->data.sff->rheader->clip_left);
      }
      else if(1 == strand && 0 < seq->data.sff->rheader->clip_right) {
          tmap_file_fprintf(fp, "%dH", seq->data.sff->rheader->clip_right);
      }
  }
  for(i=0;i<n_cigar;i++) {
      tmap_file_fprintf(fp, "%d%c",
                        cigar[i]>>4, "MIDNSHP"[cigar[i]&0xf]);
  }
  if(TMAP_SEQ_TYPE_SFF == seq->type) {
      if(1 == strand && 0 < seq->data.sff->rheader->clip_left) {
          tmap_file_fprintf(fp, "%dH", seq->data.sff->rheader->clip_left);
      }
      else if(0 == strand && 0 < seq->data.sff->rheader->clip_right) {
          tmap_file_fprintf(fp, "%dH", seq->data.sff->rheader->clip_right);
      }
  }
  
  // mate info
  if(0 == end_num) { // no mate
      tmap_file_fprintf(fp, "\t*\t0\t0");
  }
  else if(1 == m_unmapped) { // unmapped mate
      tmap_file_fprintf(fp, "\t%s\t%u\t%u",
                        "=",
                        pos + 1,
                        0);
  }
  else { // mapped mate
      tmap_file_fprintf(fp, "\t%s\t%u\t%d",
                        refseq->annos[m_seqid].name->s,
                        m_pos+1,
                        m_tlen);
  }

  // bases and qualities
  if(1 == seq_eq && NULL != bases_eq) {
      tmap_file_fprintf(fp, "\t%s\t%s", bases_eq, (0 == qualities->l) ? "*" : qualities->s);
  }
  else {
      tmap_file_fprintf(fp, "\t%s\t%s", bases->s, (0 == qualities->l) ? "*" : qualities->s);
  }
  
  // RG 
  tmap_sam_print_rg(fp, seq);

  // PG
  tmap_file_fprintf(fp, "\tPG:Z:%s", PACKAGE_NAME);

  // MD and NM
  tmap_file_fprintf(fp, "\tMD:Z:%s\tNM:i:%d", md->s, nm);

  // AS
  tmap_file_fprintf(fp, "\tAS:i:%d", score);

  // NH
  if(1 < nh) tmap_file_fprintf(fp, "\tNH:i:%d", nh);
  
  // FZ and ZF
  if(1 == sam_flowspace_tags) {
      tmap_sam_print_fz_and_zf(fp, seq);
  }

  // XA
  if(0 < algo_stage) {
      tmap_file_fprintf(fp, "\tXA:Z:%s-%d", tmap_algo_id_to_name(algo_id), algo_stage);
  }
  
  // XZ
  if(TMAP_SEQ_TYPE_SFF == seq->type && INT32_MIN != ascore) {
      tmap_file_fprintf(fp, "\tXZ:i:%d", ascore);
  }
  
  if(0 < end_num) { // mate info
      tmap_file_fprintf(fp, "\tYP:i:%d", pscore);
      if(0 == m_unmapped) {
          tmap_file_fprintf(fp, "\tYS:f:%f", m_num_std);
      }
  }
  if(1 == bidirectional) {
      tmap_file_fprintf(fp, "\tXB:i:1");
  }

  // optional tags
  if(NULL != format) {
      va_start(ap, format);
      tmap_file_vfprintf(fp, format, ap);
      va_end(ap);
  }
  // new line
  tmap_file_fprintf(fp, "\n");
  if(1 == strand) { // reverse back
      tmap_string_reverse_compliment(bases, 0);
      tmap_string_reverse(qualities);
  }

  // free
  tmap_string_destroy(md);
  free(bases_eq);
}

#ifdef HAVE_SAMTOOLS
// from bam_md.c in SAMtools
// modified not fill in the NM tag, and not to start the reference a c->pos
static void 
tmap_sam_md1_core(bam1_t *b, char *ref)
{
  uint8_t *seq = bam1_seq(b);
  uint32_t *cigar = bam1_cigar(b);
  bam1_core_t *c = &b->core;
  int i, x, y, u = 0;
  kstring_t *str;
  uint8_t *old_md, *old_nm;
  int32_t old_nm_i=-1, nm=0;

  str = (kstring_t*)calloc(1, sizeof(kstring_t));
  for (i = y = x = 0; i < c->n_cigar; ++i) {
      int j, l = cigar[i]>>4, op = cigar[i]&0xf;
      if (op == BAM_CMATCH) {
          for (j = 0; j < l; ++j) {
              int z = y + j;
              int c1 = bam1_seqi(seq, z), c2 = bam_nt16_table[(int)ref[x+j]];
              if (ref[x+j] == 0) break; // out of boundary
              if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
                  ++u;
              } else {
                  ksprintf(str, "%d", u);
                  kputc(ref[x+j], str);
                  u = 0; 
                  nm++;
              }
          }
          if (j < l) break;
          x += l; y += l;
      } else if (op == BAM_CDEL) {
          ksprintf(str, "%d", u);
          kputc('^', str);
          for (j = 0; j < l; ++j) {
              if (ref[x+j] == 0) break;
              kputc(ref[x+j], str);
          }
          u = 0;
          if (j < l) break;
          x += l; 
          nm += l;
      } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
          y += l;
          if (op == BAM_CINS) nm += l;
      } else if (op == BAM_CREF_SKIP) {
          x += l;
      }
  }
  ksprintf(str, "%d", u);

  // update MD
  old_md = bam_aux_get(b, "MD");
  if(NULL == old_md) {
      bam_aux_append(b, "MD", 'Z', str->l + 1, (uint8_t*)str->s);
  }
  else {
      int is_diff = 0;
      if(strlen((char*)old_md+1) == str->l) {
          for(i = 0; i < str->l; ++i) {
            if(toupper(old_md[i+1]) != toupper(str->s[i])) {
              break;
            }
          }
          if(i < str->l) {
              is_diff = 1;
          }
      } 
      else {
          is_diff = 1;
      }
      if(1 == is_diff) {
          bam_aux_del(b, old_md);
          bam_aux_append(b, "MD", 'Z', str->l + 1, (uint8_t*)str->s);
      }
  }

  // update NM
  old_nm = bam_aux_get(b, "NM");
  if(NULL != old_nm) {
      old_nm_i = bam_aux2i(old_nm);
      if(old_nm_i != nm) {
          bam_aux_del(b, old_nm);
          bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm);
      }
  }

  free(str->s); free(str);
}

// from bam_md.c in SAMtools
void 
tmap_sam_md1(bam1_t *b, char *ref, int32_t len)
{
  int32_t i, j;
  char *ref_tmp = NULL;
  ref_tmp = tmap_malloc(sizeof(char) * (1 + len), "ref_tmp");
  for(i=j=0;i<len;i++) {
      if('-' != ref[i] && 'H' != ref[i]) {
          ref_tmp[j] = ref[i];
          j++;
      }
  }
  ref_tmp[j]='\0';
  tmap_sam_md1_core(b, ref_tmp);
  free(ref_tmp);
}

// soft-clipping is not supported
static inline int
tmap_sam_get_type(char ref, char read)
{
  if('-' == ref) { // insertion
      return TMAP_SW_FROM_I;
  }
  else if('-' == read) { // deletion
      return TMAP_SW_FROM_D;
  }
  else { // match/mismatch
      return TMAP_SW_FROM_M;
  }
}

void
tmap_sam_update_cigar_and_md(bam1_t *b, char *ref, char *read, int32_t len)
{
  int32_t i, n_cigar, last_type;
  uint32_t *cigar;
  int32_t diff;
  int32_t soft_clip_start_i, soft_clip_end_i;

  if(b->data_len - b->l_aux != bam1_aux(b) - b->data) {
      tmap_error("b->data_len - b->l_aux != bam1_aux(b) - b->data", Exit, OutOfRange);
  }

  // keep track of soft clipping
  n_cigar = soft_clip_start_i = soft_clip_end_i = 0;
  cigar = bam1_cigar(b);
  if(BAM_CSOFT_CLIP == TMAP_SW_CIGAR_OP(cigar[0])) {
      soft_clip_start_i = 1;
      n_cigar++;
  }
  if(1 < b->core.n_cigar && BAM_CSOFT_CLIP == TMAP_SW_CIGAR_OP(cigar[b->core.n_cigar-1])) {
      soft_clip_end_i = 1;
      n_cigar++;
  }
  cigar = NULL;

  // get the # of cigar operators
  last_type = tmap_sam_get_type(ref[0], read[0]);
  n_cigar++;
  for(i=1;i<len;i++) {
      int32_t cur_type = tmap_sam_get_type(ref[i], read[i]);
      if(cur_type != last_type) {
          n_cigar++;
      }
      last_type = cur_type;
  }

  // resize the data field if necessary
  if(n_cigar < b->core.n_cigar) {
      diff = sizeof(uint32_t) * (b->core.n_cigar - n_cigar);
      // shift down
      for(i=b->core.l_qname;i<b->data_len - diff;i++) {
          b->data[i] = b->data[i + diff];
      }
      b->data_len -= diff;
      b->core.n_cigar = n_cigar;
  }
  else if(b->core.n_cigar < n_cigar) {
      diff = sizeof(uint32_t) * (n_cigar - b->core.n_cigar);
      // realloc
      if(b->m_data <= (b->data_len + diff)) {
          b->m_data = b->data_len + diff + 1;
          tmap_roundup32(b->m_data);
          b->data = tmap_realloc(b->data, sizeof(uint8_t) * b->m_data, "b->data");
      }
      // shift up
      for(i=b->data_len-1;b->core.l_qname<=i;i--) {
          b->data[i + diff] = b->data[i];
      }
      b->data_len += diff;
      b->core.n_cigar = n_cigar;
  }
  if(b->data_len - b->l_aux != bam1_aux(b) - b->data) {
      tmap_error("b->data_len - b->l_aux != bam1_aux(b) - b->data", Exit, OutOfRange);
  }

  // create the cigar
  cigar = bam1_cigar(b);
  for(i=soft_clip_start_i;i<n_cigar-soft_clip_end_i;i++) {
      cigar[i] = 0;
  }
  n_cigar = soft_clip_start_i; // skip over soft clipping etc.
  last_type = tmap_sam_get_type(ref[0], read[0]);
  TMAP_SW_CIGAR_STORE(cigar[n_cigar], last_type, 1);
  for(i=1;i<len;i++) {
      int32_t cur_type = tmap_sam_get_type(ref[i], read[i]);
      if(cur_type == last_type) {
          // add to the cigar length
          TMAP_SW_CIGAR_ADD_LENGTH(cigar[n_cigar], 1);
      }
      else {
          // add to the cigar
          n_cigar++;
          TMAP_SW_CIGAR_STORE(cigar[n_cigar], cur_type, 1);
      }
      last_type = cur_type;
  }

  // Note: the md tag must be updated
  tmap_sam_md1(b, ref, len);
}
#endif
