#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <config.h>
#include "fmap_alloc.h"
#include "../io/fmap_file.h"
#include "fmap_sam.h"

void
fmap_sam_print_header(fmap_file_t *fp, fmap_refseq_t *refseq, int argc, char *argv[])
{
  int32_t i;
  // SAM header
  for(i=0;i<refseq->num_annos;i++) {
      fmap_file_fprintf(fp, "@SQ\tSN:%s\tLN:%d\n",
                        refseq->annos[i].name->s, (int)refseq->annos[i].len);
  }
  fmap_file_fprintf(fp, "@PG\tID:%s\tVN:%s\tCL:",
                    PACKAGE_NAME, PACKAGE_VERSION);
  for(i=0;i<argc;i++) {
      if(0 < i) fmap_file_fprintf(fp, " ");
      fmap_file_fprintf(fp, "%s", argv[i]);
  }
  fmap_file_fprintf(fp, "\n");
}

inline void
fmap_sam_print_unmapped(fmap_file_t *fp, fmap_seq_t *seq)
{
  uint16_t flag = 0x0004;
  fmap_string_t *name=NULL, *bases=NULL, *qualities=NULL;

  name = fmap_seq_get_name(seq);
  bases = fmap_seq_get_bases(seq);
  qualities = fmap_seq_get_qualities(seq);

  fmap_file_fprintf(fp, "%s\t%u\t%s\t%u\t%u\t*\t*\t0\t0\t%s\t%s\n",
                    name->s, flag, "*",
                    0, 0, bases->s, qualities->s);
}

static inline fmap_string_t *
fmap_sam_md(fmap_refseq_t *refseq, char *read_bases, // read bases are characters
            uint32_t seqid, uint32_t pos, // seqid and pos are 0-based
            uint32_t *cigar, int32_t n_cigar)
{
  int32_t i, j;
  int32_t ref_i, read_i;
  int32_t l = 0; // the length of the last md op
  uint8_t read_base, ref_base;
  fmap_string_t *md=NULL;

  md = fmap_string_init(32);

  read_i = 0;
  ref_i = refseq->annos[seqid].offset + pos;
              
  for(i=0;i<n_cigar;i++) { // go through each cigar operator
      int32_t op_len, op;
              
      op_len = cigar[i] >> 4;
      op = cigar[i] & 0xf;

      if(BAM_CMATCH == op) {
          for(j=0;j<op_len;j++) {
              if(refseq->len <= ref_i) break; // out of boundary

              read_base = nt_char_to_int[(int)read_bases[read_i]]; 
              ref_base = fmap_refseq_seq_i(refseq, ref_i);

              if(read_base == ref_base) { // a match
                  l++;
              }
              else {
                  fmap_string_lsprintf(md, md->l, "%d%c", l, "ACGTN"[ref_base]);
                  l = 0;
              }
              read_i++;
              ref_i++; 
          }
          if(j < op_len) break;
      }
      else if(BAM_CINS == op) {
          read_i += op_len;
      }
      else if(BAM_CDEL == op) {
          fmap_string_lsprintf(md, md->l, "%d", l);
          for(j=0;j<op_len;j++) {
              if(refseq->len <= ref_i) break; // out of boundary
              ref_base = fmap_refseq_seq_i(refseq, ref_i);
              fmap_string_lsprintf(md, md->l, "%c", "ACGTN"[ref_base]);
              ref_i++;
          }
          if(j < op_len) break;
          l=0;
      }
      else if(BAM_CREF_SKIP == op) {
          ref_i += op_len;
      }
      else if(BAM_CSOFT_CLIP == op) {
          read_i += op_len;
      }
      else if(BAM_CHARD_CLIP == op) {
          // ignore
      }
      else if(BAM_CPAD == op) {
          // ignore
      }
      else {
          fmap_error("could not understand the cigar operator", Exit, OutOfRange);
      }
  }
  fmap_string_lsprintf(md, md->l, "%d", l);

  return md;
}

inline void
fmap_sam_print_mapped(fmap_file_t *fp, fmap_seq_t *seq, fmap_refseq_t *refseq,
                      uint8_t strand, uint32_t seqid, uint32_t pos, 
                      uint8_t mapq, uint32_t *cigar, int32_t n_cigar,
                      const char *format, ...)
{
  va_list ap;
  int32_t i, sff_soft_clip = 0;
  fmap_string_t *name=NULL, *bases=NULL, *qualities=NULL;
  uint32_t *cigar_tmp = NULL, cigar_tmp_allocated = 0;
  fmap_string_t *md;

  name = fmap_seq_get_name(seq);
  bases = fmap_seq_get_bases(seq);
  qualities = fmap_seq_get_qualities(seq);

  if(1 == strand) { // reverse for the output
      fmap_string_reverse_compliment(bases, 0);
      fmap_string_reverse(qualities);
  }

  fmap_file_fprintf(fp, "%s\t%u\t%s\t%u\t%u\t",
                    name->s, (1 == strand) ? 0x10 : 0, refseq->annos[seqid].name->s,
                    pos + 1,
                    mapq);

  // add the soft clipping of from an SFF
  cigar_tmp = cigar;
  if(FMAP_SEQ_TYPE_SFF == seq->type) {
      sff_soft_clip = seq->data.sff->gheader->key_length; // soft clip the key sequence
      if(0 < sff_soft_clip && 0 < n_cigar) {
          if(0 == strand) {  // forward strand sff soft clip
              if(4 == cigar[0]) { // add to an existing soft clip
                  cigar[0] = (((cigar[0]>>4) + sff_soft_clip) << 4) | 4;
              }
              else { // add a new soft clip to the front
                  cigar_tmp_allocated = 1;
                  cigar_tmp = fmap_calloc(n_cigar+1, sizeof(uint32_t), "cigar_tmp");
                  cigar_tmp[0] = (sff_soft_clip << 4) | 4; 
                  for(i=0;i<n_cigar;i++) {
                      cigar_tmp[i+1] = cigar[i];
                  }
                  n_cigar++;
              }
          }
          else {  // reverse strand sff soft clip
              if(4 == cigar[n_cigar-1]) { // add to an existing soft clip
                  cigar[n_cigar-1] = (((cigar[n_cigar-1]>>4) + sff_soft_clip) << 4) | 4;
              }
              else { // add a new soft clip to the end
                  cigar_tmp_allocated = 1;
                  cigar_tmp = fmap_calloc(n_cigar+1, sizeof(uint32_t), "cigar_tmp");
                  cigar_tmp[n_cigar] = (sff_soft_clip << 4) | 4; 
                  for(i=0;i<n_cigar;i++) {
                      cigar_tmp[i] = cigar[i];
                  }
                  n_cigar++;
              }
          }
      }
  }
  
  // print out the cigar
  for(i=0;i<n_cigar;i++) {
      fmap_file_fprintf(fp, "%d%c",
                        cigar_tmp[i]>>4, "MIDNSHP"[cigar_tmp[i]&0xf]);
  }

  // bases and qualities
  fmap_file_fprintf(fp, "\t*\t0\t0\t%s\t%s",
                    bases->s, qualities->s);

  // MD
  md = fmap_sam_md(refseq, bases->s, seqid, pos, cigar_tmp, n_cigar);
  fmap_file_fprintf(fp, "\tMD:Z:%s\n", md->s);
  fmap_string_destroy(md);

  // optional tags
  if(NULL != format) {
      va_start(ap, format);
      fmap_file_vfprintf(fp, format, ap);
      va_end(ap);
  }
  // new line
  fmap_file_fprintf(fp, "\n");
  if(1 == strand) { // reverse back
      fmap_string_reverse_compliment(bases, 0);
      fmap_string_reverse(qualities);
  }

  if(1 == cigar_tmp_allocated) {
      free(cigar_tmp);
  }
}
