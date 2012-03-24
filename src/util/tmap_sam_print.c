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

static char tmap_sam_rg_id[1024]="ID";

#define TMAP_SAM_PRINT_RG_HEADER_TAGS 12

// Notes: we could add all the tags as input
static void
tmap_sam_parse_rg(char *rg, const char *fo, const char *ks, const char *pg)
{
  int32_t i, j, len;
  // ID, CN, DS, DT, LB, PG, PI, PL, PU, SM
  int32_t tags_found[TMAP_SAM_PRINT_RG_HEADER_TAGS] = {0,0,0,0,0,0,0,0,0,0,0,0};
  char *tags_name[TMAP_SAM_PRINT_RG_HEADER_TAGS] = {"ID","CN","DS","DT","FO","KS","LB","PG","PI","PL","PU","SM"};
  char *tags_value[TMAP_SAM_PRINT_RG_HEADER_TAGS] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

  len = strlen(rg);

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
          int32_t tag_i = -1;
          i++; // move past the tab
          if(len <= i+2) { // must have "XX:" 
              tmap_error("Improper tag in the RG line", Exit, OutOfRange);
          }
          else if('I' == rg[i] && 'D' == rg[i+1]) {
              tag_i = 0;
              // copy over the id
              for(j=i+3;j<len;j++) {
                  if('\t' == rg[j]) break;
                  tmap_sam_rg_id[j-i-3] = rg[j];
              }
              if(j == i) tmap_error("Malformed RG line", Exit, OutOfRange);
              tmap_sam_rg_id[j-i]='\0'; // null terminator
          }
          else if('C' == rg[i] && 'N' == rg[i+1]) tag_i=1;
          else if('D' == rg[i] && 'S' == rg[i+1]) tag_i=2;
          else if('D' == rg[i] && 'T' == rg[i+1]) tag_i=3;
          else if('F' == rg[i] && 'O' == rg[i+1]) {
              tag_i=4;
              if(NULL != fo) {
                  tmap_error("FO tag not allowed in the RG line", Exit, OutOfRange);
              }
          }
          else if('K' != rg[i] && 'S' == rg[i+1]) {
              tag_i=5;
              if(NULL == ks) {
                  tmap_error("KS tag not allowed in the RG line", Exit, OutOfRange);
              }
          }
          else if('L' == rg[i] && 'B' == rg[i+1]) tag_i=6;
          else if('P' == rg[i] && 'G' == rg[i+1]) {
              tag_i=7;
              if(NULL != pg) {
                  tmap_error("PG tag not allowed in the RG line", Exit, OutOfRange);
              }
          }
          else if('P' == rg[i] && 'I' == rg[i+1]) tag_i=8;
          else if('P' == rg[i] && 'L' == rg[i+1]) tag_i=9;
          else if('P' == rg[i] && 'U' == rg[i+1]) tag_i=10;
          else if('S' == rg[i] && 'M' == rg[i+1]) tag_i=11;
          else {
              tmap_error("Improper tag in the RG line", Exit, OutOfRange);
          }
          tags_found[tag_i]++;
          if(1 == tags_found[tag_i]) {
              tags_value[tag_i] = tmap_malloc(sizeof(char) * (len + 1), "tags_value[tag_i]");
              for(j=i;j<len && '\t' != rg[j];j++) {
                  tags_value[tag_i][j-i] = rg[j];
              }
              if(j - i <= 3) {
                  tmap_file_fprintf(tmap_file_stderr, "\nFound an empty tag in the RG SAM header: %s\n", tags_name[tag_i]);
                  tmap_error(NULL, Exit, OutOfRange);
              }
              tags_value[tag_i][j-i] = '\0';
          }
      }
      i++;
  }

  strcpy(rg, "@RG");
  for(i=0;i<TMAP_SAM_PRINT_RG_HEADER_TAGS;i++) {
      if(1 < tags_found[i]) {
          tmap_file_fprintf(tmap_file_stderr, "\nFound multiple %s tags for the RG SAM header\n", tags_name[i]);
          tmap_error(NULL, Exit, OutOfRange);
      }
      else if(1 == tags_found[i]) {
          strcat(rg, "\t");
          strcat(rg, tags_value[i]);
          // copy over the tag
          free(tags_value[i]);
      }
  }
}

void
tmap_sam_print_header(tmap_file_t *fp, tmap_refseq_t *refseq, tmap_seq_io_t *seqio, char *sam_rg, 
                      char *flow_order, char *key_seq, int32_t sam_sff_tags, int argc, char *argv[])
{
  int32_t i;
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
  if(NULL != seqio && 1 == sam_sff_tags) {
      if(NULL != flow_order) { // this should not happen, since it should be checked upstream
          tmap_error("flow order was specified when using sam sff tags", Exit, OutOfRange);
      }
      if(NULL != key_seq) { // this should not happen, since it should be checked upstream
          tmap_error("key sequence was specified when using sam sff tags", Exit, OutOfRange);
      }
      if(NULL == (flow_order = tmap_seq_io_get_rg_fo(seqio))) {
          tmap_error("flow order could not be retrieved from the input", Exit, OutOfRange);
      }
      if(NULL == (key_seq = tmap_seq_io_get_rg_ks(seqio))) {
          tmap_error("key sequence could not be retrieved from the input", Exit, OutOfRange);
      }
      if(NULL != sam_rg) { // SAM RG is user-specified
          tmap_sam_parse_rg(sam_rg, 
                            flow_order,
                            key_seq,
                            PACKAGE_NAME);
          tmap_file_fprintf(fp, "%s\tFO:%s\tKS:%s\tPG:%s\n",
                            sam_rg,
                            flow_order,
                            key_seq,
                            PACKAGE_NAME);
      }
      else {
          tmap_file_fprintf(fp, "@RG\tID:%s\tFO:%s\tKS:%s\tPG:%s\n",
                            tmap_sam_rg_id,
                            flow_order,
                            key_seq,
                            PACKAGE_NAME);
      }
  }
  else {
      if(NULL != sam_rg) {
          tmap_sam_parse_rg(sam_rg, NULL, flow_order, PACKAGE_NAME);
          tmap_file_fprintf(fp, "%s\n", sam_rg);
      }
      else if(NULL == flow_order && NULL == key_seq) {
          tmap_file_fprintf(fp, "@RG\tID:%s\tPG:%s\n",
                            tmap_sam_rg_id,
                            PACKAGE_NAME);
      }
      else {
          tmap_file_fprintf(fp, "@RG\tID:%s", tmap_sam_rg_id);
          if(NULL != flow_order) {
              tmap_file_fprintf(fp, "\tFO:%s", flow_order);
          }
          if(NULL != key_seq) {
              tmap_file_fprintf(fp, "\tKS:%s", key_seq);
          }
          tmap_file_fprintf(fp, "\tPG:%s\n", PACKAGE_NAME);
      }
  }
  tmap_file_fprintf(fp, "@PG\tID:%s\tVN:%s\tCL:",
                    PACKAGE_NAME, PACKAGE_VERSION);
  for(i=0;i<argc;i++) {
      if(0 < i) tmap_file_fprintf(fp, " ");
      tmap_file_fprintf(fp, "%s", argv[i]);
  }
  tmap_file_fprintf(fp, "\n");
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
tmap_sam_print_unmapped(tmap_file_t *fp, tmap_seq_t *seq, int32_t sam_sff_tags, int32_t bidirectional, tmap_refseq_t *refseq,
                      uint32_t end_num, uint32_t m_unmapped, uint32_t m_prop, 
                      uint32_t m_strand, uint32_t m_seqid, uint32_t m_pos)
{
  uint32_t flag = 0;
  tmap_string_t *name=NULL, *bases=NULL, *qualities=NULL;

  name = tmap_seq_get_name(seq);
  bases = tmap_seq_get_bases(seq);
  qualities = tmap_seq_get_qualities(seq);
  
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
  // optional tags
  tmap_file_fprintf(fp, "\tRG:Z:%s\tPG:Z:%s",
                    tmap_sam_rg_id,
                    PACKAGE_NAME);
  // FZ and ZF
  if(1 == sam_sff_tags) {
      tmap_sam_print_fz_and_zf(fp, seq);
  }
  if(1 == bidirectional) {
      tmap_file_fprintf(fp, "\tXB:i:1");
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
tmap_sam_print_mapped(tmap_file_t *fp, tmap_seq_t *seq, int32_t sam_sff_tags, int32_t bidirectional, int32_t seq_eq, tmap_refseq_t *refseq,
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
  
  // RG and PG
  tmap_file_fprintf(fp, "\tRG:Z:%s\tPG:Z:%s",
                    tmap_sam_rg_id,
                    PACKAGE_NAME);

  // MD and NM
  tmap_file_fprintf(fp, "\tMD:Z:%s\tNM:i:%d", md->s, nm);

  // AS
  tmap_file_fprintf(fp, "\tAS:i:%d", score);

  // NH
  if(1 < nh) tmap_file_fprintf(fp, "\tNH:i:%d", nh);
  
  // FZ and ZF
  if(1 == sam_sff_tags) {
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
