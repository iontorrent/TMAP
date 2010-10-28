#ifndef FMAP_SAM2FS_H_
#define FMAP_SAM2FS_H_

/*!
  @brief
  @param  bam          the SAM/BAM structure to examine
  @param  flow_order    the flow order of the four DNA bases
  @param  flow_score    the flow score for the flow-space Smith-Waterman re-alignment
  @param  flow_offfset  the maximum homopolymer offset to examine
  @param  aln_global   the global alignment will be used if 1, fitting alignment otherwise
  @param  ref          the returned reference alignment string
  @param  read         the returned read alignment string 
  @param  aln          the returned reference/read alignment string
  */
void 
fmap_sam2fs_aux(bam1_t *bam, char *flow_order, int32_t flow_score, int32_t flow_offset, int32_t aln_global,
                                char **ref, char **read, char **aln);

/*! 
  @brief        main-like function for 'fmap sam2fs'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
fmap_sam2fs_main(int argc, char *argv[]);

#endif
