#ifndef FMAP_SAM2FS_AUX_H_
#define FMAP_SAM2FS_AUX_H_

#include "../io/fmap_file.h"

// TODO: document
void
fmap_sam2fs_aux_flow_align(fmap_file_t *fp, uint8_t *qseq, int32_t qseq_len, uint8_t *tseq, int32_t tseq_len, uint8_t *flow_order, int8_t strand);

#endif
