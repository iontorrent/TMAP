#ifndef FMAP_MAP2_MEMPOOL_H_
#define FMAP_MAP2_MEMPOOL_H_

#include <stdint.h>
#include "../util/fmap_vec.h"

struct __fmap_map2_mempool_t;

// TODO: document
typedef struct {
    uint32_t qk, ql;
    int32_t I, D, G;
    uint32_t pj:2, qlen:30;
    int32_t tlen;
    int32_t ppos, upos;
    int32_t cpos[4];
} fmap_map2_cell_t;

// TODO: document
typedef struct {
    int32_t n, max;
    uint32_t tk, tl;
    fmap_map2_cell_t *array;
} fmap_map2_entry_t, *fmap_map2_entry_p;

// TODO: document
typedef struct {
    int32_t n_pending;
    fmap_vec_t(fmap_map2_entry_p) stack0, pending;
    struct __fmap_map2_mempool_t *pool;
} fmap_map2_stack_t;

// TODO: document
typedef struct __fmap_map2_mempool_t {
    int32_t cnt; // if cnt!=0, then there must be memory leak
    fmap_vec_t(fmap_map2_entry_p) pool;
} fmap_map2_mempool_t;

// TODO: document
typedef struct {
    fmap_map2_stack_t *stack;
    int32_t max_l;
    uint8_t *aln_mem;
} fmap_map2_global_mempool_t;

// TODO: document
void 
fmap_map2_stack_destroy(fmap_map2_stack_t *stack);

// TODO: document
inline void
fmap_map2_stack_push0(fmap_map2_stack_t *s, fmap_map2_entry_p e);

// TODO: document
inline fmap_map2_entry_p 
fmap_map2_stack_pop(fmap_map2_stack_t *s);

// TODO: document
#define fmap_map2_stack_isempty(s) (fmap_vec_size(s->stack0) == 0 && s->n_pending == 0)

// TODO: document
inline fmap_map2_entry_p
fmap_map2_mempool_alloc(fmap_map2_mempool_t *mp);

// TODO: document
void
fmap_map2_mempool_free(fmap_map2_mempool_t *mp, fmap_map2_entry_t *e);

// TODO: document
void
fmap_map2_mempool_destroy(fmap_map2_mempool_t *mp);

// TODO: document
fmap_map2_mempool_t *
fmap_map2_mempool_init();

// TODO: document
fmap_map2_global_mempool_t *
fmap_map2_global_mempool_init();

// TODO: document
void
fmap_map2_global_mempool_destroy(fmap_map2_global_mempool_t *global);

#endif
