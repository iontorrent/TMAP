#ifndef FMAP_MAP2_AUX_H_
#define FMAP_MAP2_AUX_H_

#include "../util/fmap_string.h"
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_sa.h"
#include "fmap_map2_mempool.h"
#include "fmap_map2.h"

#define FMAP_MAP2_MASK_LEVEL 0.90f

// TODO
typedef struct {
    uint32_t k, l, flag:18, n_seeds:14;
    int32_t len, G, G2;
    int32_t beg, end;
} fmap_map2_hit_t;

// TODO
typedef struct {
    int32_t n, max;
    fmap_map2_hit_t *hits;
    int32_t *n_cigar;
    uint32_t **cigar;
} fmap_map2_aln_t;

// TODO
typedef struct {
    uint8_t strand:1; // 1-bit
    uint32_t seqid, pos; // zero-based
    uint8_t mapq;
    int32_t n_cigar;
    uint32_t *cigar;
    int32_t AS;
    int32_t XS;
    uint16_t XF:2, XE:14;
    int32_t XI;
} fmap_map2_sam_entry_t;

// TODO
// everything but the NAME, SEQ, and QUAL fields
typedef struct __fmap_map2_sam_t {
    int32_t num_entries;
    fmap_map2_sam_entry_t *entries;
} fmap_map2_sam_t;

// TODO
void
fmap_map2_aln_destroy(fmap_map2_aln_t *a);

// TODO
int32_t
fmap_map2_aux_resolve_duphits(const fmap_bwt_t *bwt, const fmap_sa_t *sa, fmap_map2_aln_t *b, int32_t IS);

// TODO
fmap_map2_sam_t *
fmap_map2_sam_init(int32_t n);

// TODO
fmap_map2_sam_t *
fmap_map2_sam_realloc(fmap_map2_sam_t *sam, int32_t n);

// TODO
void
fmap_map2_sam_destroy(fmap_map2_sam_t *sam);

// TODO
fmap_map2_sam_t *
fmap_map2_aux_core(fmap_map2_opt_t *_opt,
                   fmap_seq_t *p,
                   fmap_refseq_t *refseq,
                   fmap_bwt_t *bwt[2],
                   fmap_sa_t *sa[2],
                   fmap_map2_global_mempool_t *pool);

#endif
