#ifndef FMAP_SAM2FS_H_
#define FMAP_SAM2FS_H_

#include <config.h>

// TODO: doc
enum {
    FMAP_SAM2FS_OUTPUT_ALN = 0,
    FMAP_SAM2FS_OUTPUT_SAM = 1,
    FMAP_SAM2FS_OUTPUT_BAM = 2
};

typedef struct {
    char *flow_order; /*!< flow order (-f) */
    int32_t flow_score; /*!< flow penalty (-F) */
    int32_t flow_offset; /*!< search for homopolymer errors +- offset during re-alignment (-o) */
    int32_t aln_global; /*!< run global alignment (otherwise read fitting) (-g) */
    int32_t output_type; /*!< the output type: 0-flow space alignment 1-base space alignment 2-SAM (-z) */
} fmap_sam2fs_opt_t;

#ifdef HAVE_SAMTOOLS

/*!
  @brief
  @param  bam          the SAM/BAM structure to examine
  @param  flow_order    the flow order of the four DNA bases
  @param  flow_score    the flow score for the flow-space Smith-Waterman re-alignment
  @param  flow_offfset  the maximum homopolymer offset to examine
  @param  aln_global   the global alignment will be used if 1, fitting alignment otherwise
  @param  output_type  the output type
  @return              the original SAM/BAM, unless the output type is SAM
  @details             if the output type is 2, then the bam will be modified
  */
bam1_t *
fmap_sam2fs_aux(bam1_t *bam, char *flow_order, int32_t flow_score, int32_t flow_offset, 
                int32_t aln_global, int32_t output_type);

#endif

/*! 
  @brief        main-like function for 'fmap sam2fs'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
fmap_sam2fs_main(int argc, char *argv[]);

#endif
