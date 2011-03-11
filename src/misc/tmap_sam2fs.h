/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SAM2FS_H
#define TMAP_SAM2FS_H

#include <config.h>

/*!
 The type of output for sam2fs
 */ 
enum {
    TMAP_SAM2FS_OUTPUT_ALN = 0, /*!< pretty-print alignment */
    TMAP_SAM2FS_OUTPUT_SAM = 1, /*!< SAM file */
    TMAP_SAM2FS_OUTPUT_BAM = 2 /*!< BAM file */
};

typedef struct {
    char *flow_order; /*!< flow order (-f) */
    int32_t score_match;  /*!< the match score (-A) */
    int32_t pen_mm;  /*!< the mismatch penalty (-M) */
    int32_t pen_gapo;  /*!< the indel open penalty (-O) */
    int32_t pen_gape;  /*!< the indel extension penalty (-E) */
    int32_t fscore;  /*!< the flow score penalty (-X) */
    int32_t flow_offset; /*!< search for homopolymer errors +- offset during re-alignment (-o) */
    int32_t softclip_type; /*!< the soft clip type (-g) */
    int32_t output_type; /*!< the output type: 0-flow space alignment 1-base space alignment 2-SAM (-z) */
    int32_t j_type; /*!< how indels are justified in alignment (-l) */
    int32_t reads_queue_size;  /*!< the reads queue size (-q) */
    int32_t num_threads;  /*!< the number of threads (-n) */
} tmap_sam2fs_opt_t;

/*! 
  @brief        main-like function for 'tmap sam2fs'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
tmap_sam2fs_main(int argc, char *argv[]);

#endif
