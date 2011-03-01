/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SAM2FS_H
#define TMAP_SAM2FS_H

#include <config.h>

// TODO: doc
enum {
    TMAP_SAM2FS_OUTPUT_ALN = 0,
    TMAP_SAM2FS_OUTPUT_SAM = 1,
    TMAP_SAM2FS_OUTPUT_BAM = 2
};

typedef struct {
    char *flow_order; /*!< flow order (-f) */
    int32_t flow_score; /*!< flow penalty (-F) */
    int32_t flow_offset; /*!< search for homopolymer errors +- offset during re-alignment (-o) */
    int32_t softclip_type; /*!< the soft clip type (-g) */
    int32_t output_type; /*!< the output type: 0-flow space alignment 1-base space alignment 2-SAM (-z) */
    int32_t j_type; /*!< how indels are justified in alignment (-l) */
    int32_t reads_queue_size;  /*!< the reads queue size (-q) */
    int32_t num_threads;  /*!< the number of threads (-n) */
} tmap_sam2fs_opt_t;

#ifdef HAVE_SAMTOOLS

/*!
  @brief
  @param  bam          the SAM/BAM structure to examine
  @param  flow_order    the flow order of the four DNA bases
  @param  flow_score    the flow score for the flow-space Smith-Waterman re-alignment
  @param  flow_offfset  the maximum homopolymer offset to examine
  @param  aln_global   the global alignment will be used if 1, fitting alignment otherwise
  @param  output_type  the output type
  @param  j_type       how indels are justified in alignment 
  @return              the original SAM/BAM, unless the output type is SAM
  @details             if the output type is 2, then the bam will be modified
  */
bam1_t *
tmap_sam2fs_aux(bam1_t *bam, char *flow_order, int32_t flow_score, int32_t flow_offset, 
                int32_t aln_global, int32_t output_type, int32_t j_type);

#endif

/*! 
  @brief        main-like function for 'tmap sam2fs'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
tmap_sam2fs_main(int argc, char *argv[]);

#endif
