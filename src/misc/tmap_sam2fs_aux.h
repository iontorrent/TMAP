/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SAM2FS_AUX_H
#define TMAP_SAM2FS_AUX_H

#include "../io/tmap_file.h"

typedef struct {
    uint8_t *flow_order; /*!< the flow order in 2-bit format */
    int32_t flow_order_len; /*!< the flow order length */
    uint16_t *jump_fwd; /*!< the forward jump table */
    uint16_t *jump_rev; /*!< the reverse jump table */
} tmap_sam2fs_aux_flow_order_t;

typedef struct {
    // for flow space smith waterman
    char *fsw_qseq; /*< query bases */
    char *fsw_tseq; /*< target bases */
    char *fsw_aln; /*< flow space smith waterman alignment string */
    int32_t score; /*< flow space smith waterman score */
    // for flow space alignment
    char *sam2fs_flow_order; /*< flow order */
    int8_t *sam2fs_qseq; /*< query flow signals */
    int8_t *sam2fs_tseq; /*< target flow signals */
    char *sam2fs_aln; /*< flow space alignment string */
    int32_t sam2fs_len; /*< alignment length */
} tmap_sam2fs_aln_t;

/*!
  @param  flow_order      the flow order bases, in character format
  @param  flow_order_len  the number of bases in the flow order
  @return                the initialized flow order for sam2fs
  */
tmap_sam2fs_aux_flow_order_t*
tmap_sam2fs_aux_flow_order_init(char *flow_order);

/*!
  @param  f  the flow order to destroy, for sam2fs
  */
void
tmap_sam2fs_aux_flow_order_destroy(tmap_sam2fs_aux_flow_order_t *f);

// TODO: document
void
tmap_sam2fs_aux_flow_align(tmap_file_t *fp, uint8_t *qseq, int32_t qseq_len,
                           uint8_t *tseq, int32_t tseq_len, 
                           tmap_sam2fs_aux_flow_order_t *flow_order, 
                           int8_t strand, char sep,
                           tmap_sam2fs_aln_t *aln_ret);

#endif
