/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <config.h>
#include <math.h>

#include "../samtools/kstring.h"
#include "../samtools/sam.h"
#include "../samtools/bam.h"

#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_string.h"
#include "../io/tmap_file.h"
#include "../io/tmap_seq_io.h"
#include "../sw/tmap_sw.h"
#include "tmap_sam_convert.h"

/**
 * TODO:
 * - is there a more efficient way to query, add, and replace optional tags?
 */

static inline void
tmap_sam_convert_replace_tagA(bam1_t *b, const char tag[2], char value)
{
  uint8_t *ex_s = NULL;
  char ex_value = 0;
  // check if it exists
  ex_s = bam_aux_get(b, tag);
  if(NULL != ex_s) { // exists
      ex_value = bam_aux2A(ex_s); // get value
      if(value == ex_value) return; // OK
      else bam_aux_del(b, ex_s); // Delete
  }
  // append
  bam_aux_append(b, tag, 'A', 1, (uint8_t*)&value);
}

static inline void
tmap_sam_convert_replace_tagi(bam1_t *b, const char tag[2], int32_t value)
{
  uint8_t *ex_s = NULL;
  int32_t ex_value = 0;
  // check if it exists
  ex_s = bam_aux_get(b, tag);
  if(NULL != ex_s) { // exists
      ex_value = bam_aux2i(ex_s); // get value
      if(value == ex_value) return; // OK
      else bam_aux_del(b, ex_s); // Delete
  }
  // append
  bam_aux_append(b, tag, 'i', sizeof(int32_t), (uint8_t*)&value);
}

static inline void
tmap_sam_convert_replace_tagf(bam1_t *b, const char tag[2], float value)
{
  uint8_t *ex_s = NULL;
  float ex_value = 0;
  // check if it exists
  ex_s = bam_aux_get(b, tag);
  if(NULL != ex_s) { // exists
      ex_value = bam_aux2f(ex_s); // get value
      if(!(value != ex_value)) return; // OK
      else bam_aux_del(b, ex_s); // Delete
  }
  // append
  bam_aux_append(b, tag, 'f', sizeof(float), (uint8_t*)&value);
}

static inline void
tmap_sam_convert_replace_tagZ(bam1_t *b, const char tag[2], const char *value)
{
  uint8_t *ex_s = NULL;
  char *ex_value = NULL;
  // check if it exists
  ex_s = bam_aux_get(b, tag);
  if(NULL != ex_s) { // exists
      ex_value = bam_aux2Z(ex_s); // get value
      if(0 == strcmp(value, ex_value)) return; // OK
      else bam_aux_del(b, ex_s); // Delete
  }
  // append
  bam_aux_append(b, tag, 'Z', 1+strlen(value), (uint8_t*)value);
}

static inline void
tmap_sam_convert_replace_tagB_S(bam1_t *b, const char tag[2], int32_t len, const uint16_t *value)
{
  uint8_t *ex_s = NULL;
  uint16_t *ex_value = NULL;
  int32_t ex_len = 0;
  // check if it exists
  ex_s = bam_aux_get(b, tag);
  if(NULL != ex_s) { // exists
      ex_value = bam_auxB2S(ex_s, &len);
      if(len == ex_len && 0 == memcmp(ex_value, value, len)) return; // OK
      else bam_aux_del(b, ex_s); // Delete
  }
  // append
  bam_aux_appendB(b, tag, 'B', 'S', len, (uint8_t*)value);
}

// adapted from src/samtools/padding.c
static void 
tmap_sam_convert_replace_cigar(bam1_t *b, int n, uint32_t *cigar)
{
  if (n != b->core.n_cigar) {
      int o = b->core.l_qname + b->core.n_cigar * 4;
      if (b->data_len + (n - b->core.n_cigar) * 4 > b->m_data) {
          b->m_data = b->data_len + (n - b->core.n_cigar) * 4;
          tmap_roundup32(b->m_data);
          b->data = (uint8_t*)realloc(b->data, b->m_data);
      }
      memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->data_len - o);
      memcpy(b->data + b->core.l_qname, cigar, n * 4);
      b->data_len += (n - b->core.n_cigar) * 4;
      b->core.n_cigar = n;
  } else memcpy(b->data + b->core.l_qname, cigar, n * 4);
}

static inline void
tmap_sam_convert_rg(tmap_seq_t *seq, bam1_t *b)
{
  char *id = NULL;
  if(NULL == seq->rg_record) return;
  id = sam_header_record_get(seq->rg_record, "ID");
  if(NULL == id) return;
  // append
  tmap_sam_convert_replace_tagZ(b, "RG", id);
}

static inline void
tmap_sam_convert_pg(tmap_seq_t *seq, bam1_t *b)
{
  char *id = NULL;
  if(NULL == seq->pg_record) return;
  id = sam_header_record_get(seq->pg_record, "ID");
  if(NULL == id) return;
  // append
  tmap_sam_convert_replace_tagZ(b, "PG", id);
}

static inline void 
tmap_sam_convert_fz_and_zf(tmap_seq_t *seq, bam1_t *b)
{
  // ZF
  if(0 <= seq->fo_start_idx) { // exists
      tmap_sam_convert_replace_tagi(b, "ZF", seq->fo_start_idx);
  }

  // FZ
  if(NULL != seq->flowgram && 0 < seq->flowgram_len) {
      tmap_sam_convert_replace_tagB_S(b, "FZ", seq->flowgram_len, seq->flowgram);
  }
}

static inline void
tmap_sam_convert_add_xa(bam1_t *b, int32_t algo_id, int32_t algo_stage)
{
  char *str = NULL;
  char *name = NULL;

  name = tmap_algo_id_to_name(algo_id);
  str = tmap_malloc(sizeof(char) * (1 + strlen(name) + (int)(1.0 + log10(algo_stage)) + 1), "str"); // EOL + name + num + dash
  if(0 == sprintf(str, "%s-%d", name, algo_stage)) tmap_bug(); 
  tmap_sam_convert_replace_tagZ(b, "XA", str);
  free(str);
}

static inline void
tmap_sam_convert_add_optional(bam1_t *b, const char *format, va_list ap)
{
  int32_t start, len;

  if(NULL == format) return;

  len = strlen(format);
  start = 0;
  while(start < len) {
      char tag[2];
      char type;

      // tab
      if('\t' != format[start]) {
          start++;
          continue;
      }
      start++;

      // check for space
      if(len - start < 7) tmap_bug();

      // tag
      tag[0] = format[start]; start++;
      tag[1] = format[start]; start++;
      // colon
      if(':' != format[start]) tmap_bug();
      start++;
      // type
      type = format[start]; start++;
      // colon
      if(':' != format[start]) tmap_bug();
      start++;
      // ignore percentage sign and specifier

      // append based on type
      char A;
      int32_t i;
      float f;
      char *Z;
      switch(type) {
        case 'A':
          A = (char)va_arg(ap, int32_t);
          tmap_sam_convert_replace_tagA(b, tag, A);
          break;
        case 'i':
          i = va_arg(ap, int32_t);
          tmap_sam_convert_replace_tagi(b, tag, i);
          break;
        case 'f':
          f = (float)va_arg(ap, double);
          tmap_sam_convert_replace_tagf(b, tag, f);
          break;
        case 'Z':
          Z = va_arg(ap, char*);
          tmap_sam_convert_replace_tagZ(b, tag, Z);
          break;
        default:
          // NB: not supported
          tmap_bug();
      }
  }
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

static void
tmap_sam_convert_init_aux(bam1_t *b, tmap_string_t *qname, 
                          uint32_t *cigar, int32_t n_cigar,
                          tmap_string_t *seq, tmap_string_t *qual) 
{
  int32_t i;

  b->l_aux = 0;

  // NB: use bam1_aux to calculate the pre-aux length. So need the following
  // set:
  // 1. (b)->core.n_cigar
  b->core.n_cigar = (uint32_t)n_cigar;
  // 2. (b)->core.l_qname
  b->core.l_qname = (uint32_t)qname->l + 1; // zero-tailing
  // 3. (b)->core.l_qseq
  b->core.l_qseq = (uint32_t)seq->l;

  // TODO: missing qual string
  if(seq->l != qual->l) tmap_bug(); 

  // from bam1_aux
  b->data_len = b->m_data = b->core.n_cigar*4 + b->core.l_qname + b->core.l_qseq + (b->core.l_qseq + 1)/2;

  // re-allocate
  b->data = tmap_calloc(b->m_data, sizeof(uint8_t), "b->data");

  // qname
  for(i=0;i<qname->l;i++) {
      bam1_qname(b)[i] = qname->s[i];
  }
  bam1_qname(b)[i] = '\0'; // include trailing '\0'

  // cigar
  for(i=0;i<b->core.n_cigar;i++) {
      // TODO: validate that the format is the same
      bam1_cigar(b)[i] = cigar[i];
  }

  // seq
  for(i=0;i<seq->l;i++) {
      bam1_seq_seti(bam1_seq(b), i, bam_nt16_table[(int)seq->s[i]]); 
  }

  // qual
  for(i=0;i<qual->l;i++) {
      bam1_qual(b)[i] = qual->s[i] - 33;
  }
}

inline bam1_t*
tmap_sam_convert_unmapped(tmap_file_t *fp, tmap_seq_t *seq, int32_t sam_flowspace_tags, int32_t bidirectional, tmap_refseq_t *refseq,
                      uint32_t end_num, uint32_t m_unmapped, uint32_t m_prop, 
                      uint32_t m_strand, uint32_t m_seqid, uint32_t m_pos,
                      const char *format, ...)
{
  int32_t i;
  va_list ap;
  bam1_t *b = NULL;

  b = tmap_calloc(1, sizeof(bam1_t), "b"); 

  // From SAM/BAM
  if(TMAP_SEQ_TYPE_SAM == seq->type || TMAP_SEQ_TYPE_BAM == seq->type) {
      // copy from original
      bam1_t *o = seq->data.sam->b; 
      (*b) = (*o); // shallow copy
      // copy auxiliary data
      b->data = tmap_calloc(o->m_data, sizeof(uint8_t), "b->data");
      for(i=0;i<o->data_len;i++) {
          b->data[i] = o->data[i];
      }
      // NB: name/bases/qualities should already be set
      // check the cigar
      tmap_sam_convert_replace_cigar(b, 0, NULL);
  }
  else {
      // init the BAM structure
      tmap_sam_convert_init_aux(b, tmap_seq_get_name(seq),
                                NULL, 0, 
                                tmap_seq_get_bases(seq),
                                tmap_seq_get_qualities(seq));
  }

  // set the flag
  // NB: these are additive
  b->core.flag |= 0x4;
  if(0 < end_num) { // mate info
      b->core.flag |= 0x1;
      if(1 == m_prop) b->core.flag |= 0x2; // properly aligned
      if(1 == m_unmapped) b->core.flag |= 0x8; // unmapped
      else if(1 == m_strand) b->core.flag |= 0x20; // strand 
      b->core.flag |= (1 == end_num) ? 0x40 : 0x80; // first/second end
  }

  // tid and pos
  b->core.tid = -1;
  b->core.pos = -1;

  // mapq
  b->core.qual = 0;

  // NB: hard clipped portions of the read are not reported
  // mate info
  b->core.mtid = b->core.mpos = -1;
  b->core.isize = 0;
  if(0 == end_num) { // no mate
      // do nothing
  }
  else if(1 == m_unmapped) { // unmapped mate
      // do nothing
  }
  else if(NULL != refseq) { // mapped mate
      b->core.mtid = m_seqid;
      b->core.mpos = m_pos;
  }

  // RG 
  tmap_sam_convert_rg(seq, b);

  // PG
  tmap_sam_convert_pg(seq, b);
  
  // FZ and ZF
  if(1 == sam_flowspace_tags) {
      tmap_sam_convert_fz_and_zf(seq, b);
  }
  
  // XB
  if(1 == bidirectional) {
      tmap_sam_convert_replace_tagi(b, "XB", bidirectional);
  }
  
  // optional tags
  va_start(ap, format);
  tmap_sam_convert_add_optional(b, format, ap);
  va_end(ap);

  return b;
}

inline bam1_t*
tmap_sam_convert_mapped(tmap_file_t *fp, tmap_seq_t *seq, int32_t sam_flowspace_tags, int32_t bidirectional, int32_t seq_eq, tmap_refseq_t *refseq,
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
  tmap_string_t *md = NULL;
  int32_t nm;
  bam1_t *b = NULL;
  uint32_t *cur_cigar;

  /*
  fprintf(stderr, "end_num=%d m_unmapped=%d m_prop=%d m_strand=%d m_seqid=%d m_pos=%d m_tlen=%d\n",
          end_num, m_unmapped, m_prop, m_strand, m_seqid, m_pos, m_tlen);
  */

  // Copy the cigar, in case of clipping 
  cur_cigar = cigar;
  cigar = tmap_calloc(n_cigar, sizeof(uint32_t), "cigar");
  memcpy(cigar, cur_cigar, sizeof(uint32_t) * n_cigar);

  // Add hard clips to the cigar...
  if(TMAP_SEQ_TYPE_SFF == seq->type) {
      if(0 == strand && 0 < seq->data.sff->rheader->clip_left) {
          cigar = tmap_realloc(cigar, sizeof(uint32_t) * (n_cigar+1), "cigar");
          for(i=n_cigar;0<i;i--) {
              cigar[i] = cigar[i-1];
          }
          cigar[0] = (seq->data.sff->rheader->clip_left << 4) | BAM_CHARD_CLIP;
          n_cigar++;
      }
      else if(1 == strand && 0 < seq->data.sff->rheader->clip_right) {
          cigar = tmap_realloc(cigar, sizeof(uint32_t) * (n_cigar+1), "cigar");
          for(i=n_cigar;0<i;i--) {
              cigar[i] = cigar[i-1];
          }
          cigar[0] = (seq->data.sff->rheader->clip_right << 4) | BAM_CHARD_CLIP;
          n_cigar++;
      }
      if(1 == strand && 0 < seq->data.sff->rheader->clip_left) {
          cigar = tmap_realloc(cigar, sizeof(uint32_t) * (n_cigar+1), "cigar");
          cigar[n_cigar] = (seq->data.sff->rheader->clip_left << 4) | BAM_CHARD_CLIP;
          n_cigar++;
      }
      else if(0 == strand && 0 < seq->data.sff->rheader->clip_right) {
          cigar = tmap_realloc(cigar, sizeof(uint32_t) * (n_cigar+1), "cigar");
          cigar[n_cigar] = (seq->data.sff->rheader->clip_right << 4) | BAM_CHARD_CLIP;
          n_cigar++;
      }
  }

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

  // BAM structure
  b = tmap_calloc(1, sizeof(bam1_t), "b"); 

  // From SAM/BAM
  if(TMAP_SEQ_TYPE_SAM == seq->type || TMAP_SEQ_TYPE_BAM == seq->type) {
      // copy from original
      bam1_t *o = seq->data.sam->b; 
      // HERE END
      (*b) = (*o); // shallow copy
      // copy auxiliary data
      b->data = tmap_calloc(o->m_data, sizeof(uint8_t), "b->data");
      memcpy(b->data, o->data, sizeof(uint8_t) * o->m_data); // o->data_len
      // NB: name/bases/qualities should already be set
      // check the cigar
      if(0 < b->core.n_cigar) tmap_bug(); // TODO: reset the cigar...
      if(1 == strand) { // reset bases and qualities on the reverse strand (should be the same length)
          // seq
          for(i=0;i<bases->l;i++) {
              bam1_seq_seti(bam1_seq(b), i, bam_nt16_table[(int)bases->s[i]]); 
          }
          // qual
          for(i=0;i<qualities->l;i++) {
              bam1_qual(b)[i] = qualities->s[i] - 33;
          }
      }
      // update the cigar
      if(0 < n_cigar) {
          tmap_sam_convert_replace_cigar(b, n_cigar, cigar);
      }
  }
  else {
      // init the BAM structure
      tmap_sam_convert_init_aux(b, name,
                                cigar, n_cigar,
                                bases, qualities);
  }

  // flag
  // NB: these are additive
  b->core.flag = 0;
  if(1 == strand) b->core.flag |= 0x10; // strand
  if(0 < aln_num) b->core.flag |= 0x100; // secondary alignment
  if(0 < end_num) { // mate info
      b->core.flag |= 0x1;
      if(0 == m_unmapped && 1 == m_prop) b->core.flag |= 0x2; // properly aligned
      if(1 == m_unmapped) b->core.flag |= 0x8; // unmapped
      else if(1 == m_strand) b->core.flag |= 0x20; // strand 
      b->core.flag |= (1 == end_num) ? 0x40 : 0x80; // first/second end
  }

  // tid and pos
  b->core.tid = seqid;
  b->core.pos = pos;

  // mapq
  b->core.qual = mapq;
  
  // mate info
  b->core.mtid = b->core.mpos = -1;
  b->core.isize = 0;
  if(0 == end_num) { // no mate
      // do nothing
  }
  else if(1 == m_unmapped) { // unmapped mate
      b->core.mtid = seqid;
      b->core.mpos = pos;
  }
  else { // mapped mate
      b->core.mtid = m_seqid;
      b->core.mpos = m_pos;
      b->core.isize = m_tlen;
  }

  // TODO support 'seq_eq'
  // TODO support missing qualities
  // TODO: optional tags
  
  // RG 
  tmap_sam_convert_rg(seq, b);
  
  // PG
  tmap_sam_convert_pg(seq, b);

  // MD
  tmap_sam_convert_replace_tagZ(b, "MD", md->s);

  // NM
  tmap_sam_convert_replace_tagi(b, "NM", nm);

  // AS
  tmap_sam_convert_replace_tagi(b, "AS", score);

  // NH
  tmap_sam_convert_replace_tagi(b, "NH", nh);

  // FZ and ZF
  if(1 == sam_flowspace_tags) {
      tmap_sam_convert_fz_and_zf(seq, b);
  }
  
  // XA
  if(0 < algo_stage) {
      tmap_sam_convert_add_xa(b, algo_id, algo_stage);
  }
  
  // XZ
  if(TMAP_SEQ_TYPE_SFF == seq->type && INT32_MIN != ascore) {
      tmap_sam_convert_replace_tagi(b, "XZ", ascore);
  }
  
  // YP/YS
  if(0 < end_num) { // mate info
      tmap_sam_convert_replace_tagi(b, "YP", pscore);
      if(0 == m_unmapped) {
          tmap_sam_convert_replace_tagi(b, "YS", m_num_std);
      }
  }

  // XB
  if(1 == bidirectional) {
      tmap_sam_convert_replace_tagi(b, "XB", bidirectional);
  }

  // optional tags
  va_start(ap, format);
  tmap_sam_convert_add_optional(b, format, ap);
  va_end(ap);

  if(1 == strand) { // reverse back
      tmap_string_reverse_compliment(bases, 0);
      tmap_string_reverse(qualities);
  }

  // free
  tmap_string_destroy(md);
  free(bases_eq);
  free(cigar);

  return b;
}

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
