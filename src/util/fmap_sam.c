#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <config.h>

#ifdef HAVE_SAMTOOLS
#include <kstring.h>
#include <sam.h>
#include <bam.h>
#endif

#include "fmap_alloc.h"
#include "../io/fmap_file.h"
#include "../sw/fmap_sw.h"
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
  uint32_t ref_i, read_i;
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

              read_base = fmap_nt_char_to_int[(int)read_bases[read_i]]; 
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
          fmap_string_lsprintf(md, md->l, "%d^", l);
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
  fmap_file_fprintf(fp, "\tMD:Z:%s", md->s);
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

#ifdef HAVE_SAMTOOLS
// from bam_md.c in SAMtools
static void 
bam_fillmd1_core(bam1_t *b, char *ref, int is_equal)
{
  uint8_t *seq = bam1_seq(b);
  uint32_t *cigar = bam1_cigar(b);
  bam1_core_t *c = &b->core;
  int i, x, y, u = 0;
  kstring_t *str;
  uint8_t *old_md;

  str = (kstring_t*)calloc(1, sizeof(kstring_t));
  for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
      int j, l = cigar[i]>>4, op = cigar[i]&0xf;
      if (op == BAM_CMATCH) {
          for (j = 0; j < l; ++j) {
              int z = y + j;
              int c1 = bam1_seqi(seq, z), c2 = bam_nt16_table[(int)ref[x+j]];
              if (ref[x+j] == 0) break; // out of boundary
              if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
                  if (is_equal) seq[z/2] &= (z&1)? 0xf0 : 0x0f;
                  ++u;
              } else {
                  ksprintf(str, "%d", u);
                  kputc(ref[x+j], str);
                  u = 0; 
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
      } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
          y += l;
      } else if (op == BAM_CREF_SKIP) {
          x += l;
      }
  }
  ksprintf(str, "%d", u);

  // update MD
  old_md = bam_aux_get(b, "MD");
  if (!old_md) {
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
  free(str->s); free(str);
}

// from bam_md.c in SAMtools
static void 
bam_fillmd1(bam1_t *b, char *ref, int is_equal)
{
  bam_fillmd1_core(b, ref, is_equal);
}

// soft-clipping is not supported
static inline int
fmap_sam_get_type(char ref, char read)
{
  if('-' == ref) { // insertion
      return FMAP_SW_FROM_I;
  }
  else if('-' == read) { // deletion
      return FMAP_SW_FROM_D;
  }
  else { // match/mismatch
      return FMAP_SW_FROM_M;
  }
}

void
fmap_sam_left_justify(bam1_t *b, char *ref, char *read, int32_t len)
{
  int32_t i, j, n_cigar, last_type;
  char *ref_tmp=NULL;
  uint32_t *cigar;

  // get the # of cigar operators
  last_type = fmap_sam_get_type(ref[0], read[0]);
  for(i=n_cigar=1;i<len;i++) {
      int32_t cur_type = fmap_sam_get_type(ref[i], read[i]);
      if(cur_type != last_type) {
          n_cigar++;
      }
      last_type = cur_type;
  }

  // resize the data field if necessary
  if(n_cigar < b->core.n_cigar) {
      // shift down
      for(i=0;i<b->l_aux - b->core.n_cigar + n_cigar;i++) {
          b->data[i] = b->data[i + b->core.n_cigar - n_cigar];
      }
      b->l_aux -= b->core.n_cigar - n_cigar;
      b->core.n_cigar = n_cigar;
  }
  else if(b->core.n_cigar < n_cigar) {
      // realloc
      b->data = fmap_realloc(b->data, sizeof(uint8_t) * (b->l_aux + n_cigar - b->core.n_cigar), "b->data");
      // shift up
      for(i=b->l_aux-1;0<=i;i--) {
          b->data[i + n_cigar - b->core.n_cigar] = b->data[i];
      }
      b->l_aux += n_cigar - b->core.n_cigar;
      b->core.n_cigar = n_cigar;
  }

  // create the cigar
  cigar = bam1_cigar(b);
  for(i=0;i<n_cigar;i++) {
      cigar[0] = 0;
  }
  cigar[0] = 1u << 4 | fmap_sam_get_type(ref[0], read[0]);
  for(i=1,n_cigar=0;i<len;i++) {
      int32_t cur_type = fmap_sam_get_type(ref[i], read[i]);
      if(cur_type == last_type) {
          // add to the cigar length
          cigar[n_cigar] += 1u << 4; 
      }
      else {
          // add to the cigar
          n_cigar++;
          cigar[n_cigar] = 1u << 4 | cur_type; 
      }
      last_type = cur_type;
  }

  // update MD
  ref_tmp = fmap_malloc(sizeof(char) * (1 + len), "ref_tmp");
  for(i=j=0;i<len;i++) {
      if('-' != ref[i]) {
          ref_tmp[j] = ref[i];
          j++;
      }
  }
  ref_tmp[i]='\0';
  bam_fillmd1(b, ref_tmp, 0);
  free(ref_tmp);
}
#endif
