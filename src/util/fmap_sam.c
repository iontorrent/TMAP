#include <stdlib.h>
#include <stdio.h>
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
