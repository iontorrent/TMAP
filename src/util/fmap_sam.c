#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <config.h>
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

inline void
fmap_sam_print_mapped(fmap_file_t *fp, fmap_seq_t *seq, fmap_refseq_t *refseq,
                    uint8_t strand, uint32_t seqid, uint32_t pos, 
                    int32_t mapq, uint32_t *cigar, int32_t n_cigar,
                    const char *format, ...)
{
  va_list ap;
  int32_t i, sff_soft_clip = 0, cigar_start, cigar_end;
  fmap_string_t *name=NULL, *bases=NULL, *qualities=NULL;

  name = fmap_seq_get_name(seq);
  bases = fmap_seq_get_bases(seq);
  qualities = fmap_seq_get_qualities(seq);

  if(FMAP_SEQ_TYPE_SFF == seq->type) {
      sff_soft_clip = seq->data.sff->gheader->key_length; // soft clip the key sequence
  }

  if(1 == strand) { // reverse for the output
      fmap_string_reverse_compliment(bases, 0);
      fmap_string_reverse(qualities);
  }

  fmap_file_fprintf(fp, "%s\t%u\t%s\t%u\t%u\t",
                    name->s, (1 == strand) ? 0x10 : 0, refseq->annos[seqid].name->s,
                    pos + 1,
                    mapq);
  // Note: we must check if the cigar starts or ends with a soft clip
  cigar_start = 0;
  cigar_end = n_cigar;
  if(0 < sff_soft_clip) {
      if(0 == strand) {  // forward strand sff soft clip
          if(0 < n_cigar && 4 == cigar[0]) {
              sff_soft_clip += (cigar[0]>>4);
              cigar_start++; // do not print out the first cigar op, this will be printed out later
          }
          fmap_file_fprintf(fp, "%dS", sff_soft_clip);
      }
      else {  // reverse strand sff soft clip
          if(0 < n_cigar && 4 == cigar[cigar_end-1]) {
              sff_soft_clip += (cigar[cigar_end-1]>>4);
              cigar_end--; // do not print out the last cigar op, this will be printed out later
          }
      }
  }
  // print out the cigar
  for(i=cigar_start;i<cigar_end;i++) {
      fmap_file_fprintf(fp, "%d%c",
                        cigar[i]>>4, "MIDNSHP"[cigar[i]&0xf]);
  }
  // add trailing soft clipping if necessary
  if(0 < sff_soft_clip && 1 == strand) {  // reverse strand sff soft clip
      fmap_file_fprintf(fp, "%dS", sff_soft_clip);
  }
  fmap_file_fprintf(fp, "\t*\t0\t0\t%s\t%s",
                    bases->s, qualities->s);
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
}
