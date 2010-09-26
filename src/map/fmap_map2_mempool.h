#ifndef FMAP_MAP2_MEMPOOL_H_
#define FMAP_MAP2_MEMPOOL_H_

#include <stdint.h>
#include "../util/fmap_vec.h"

/*! @header
  @abstract  Memory Pools for Map2
  */

// TODO: document
/*! @typedef
  @abstract
  @field  qk    lower suffix array interval of the query
  @field  ql    upper suffix array interval of the query
  @field  I     insertion score
  @field  D     deletion score
  @field  G     maximum of the match/mismatch/insertion/deletion scores
  @field  pj    the zero-based index of the parent in the prefix trie [0-3]
  @field  qlen  the number of bases used in the query 
  @field  tlen  the number of bases used in the target
  @field  ppos  the parents position in the memory array
  @field  upos  the child position of its parent in the memory array
  @field  cpos  undetermined
 */
typedef struct {
    uint32_t qk, ql;
    int32_t I, D, G;
    uint32_t pj:2, qlen:30;
    int32_t tlen;
    int32_t ppos, upos;
    int32_t cpos[4];
} fmap_map2_cell_t;

/*! @typedef
  @field  n      the index of the next cell
  @field  max    the number of cells allocated
  @field  tk     lower suffix array interval of the target
  @field  tl     upper suffix array interval of the target 
  @field  array  the array of cells
*/
typedef struct {
    int32_t n, max;
    uint32_t tk, tl;
    fmap_map2_cell_t *array;
} fmap_map2_entry_t, *fmap_map2_entry_p;

/*! @typedef
  @abstract     a memory pool
  @field  cnt   the size of the memory pool
  @field  pool  the memory pool vector
*/
typedef struct __fmap_map2_mempool_t {
    int32_t cnt; // if cnt!=0, then there must be memory leak
    fmap_vec_t(fmap_map2_entry_p) pool;
} fmap_map2_mempool_t;

/*! 
  @brief            a two-level memory stack
  @param  n_pending  the number of elements in the pending stack
  @param  stack0     the main stack
  @param  pending    the pending stack
  @param  pool       a memory pool
*/
typedef struct {
    int32_t n_pending;
    fmap_vec_t(fmap_map2_entry_p) stack0, pending;
    struct __fmap_map2_mempool_t *pool;
} fmap_map2_stack_t;

/*! @typedef
  @abstract       a global memory pool
  @field  stack    the main two-level memory stack
  @field  max_l    the working memory length
  @field  aln_mem  working memory
*/
typedef struct {
    fmap_map2_stack_t *stack;
    int32_t max_l;
    uint8_t *aln_mem;
} fmap_map2_global_mempool_t;

/*! @function
  @abstract      destroys the stack
  @param  stack  a pointer to the stack
  */
void 
fmap_map2_stack_destroy(fmap_map2_stack_t *stack);

/*! @function
  @abstract      push an entry onto the main stack
  @param  stack  a pointer to the stack
  @param  e      the element to push
  */
inline void
fmap_map2_stack_push0(fmap_map2_stack_t *stack, fmap_map2_entry_p e);

/*! @function
  @abstract      pop an entry from the main stack
  @param  stack  a pointer to the stack
  @return        the popped element
  @discussion    this will fail (error) if the stack is empty
  */
inline fmap_map2_entry_p 
fmap_map2_stack_pop(fmap_map2_stack_t *stack);

/*! @define
  @abstract      checks if is the stack empty
  @param  stack  a pointer to the stack
  @return        true if the stack is empty, false otherwise 
  */
#define fmap_map2_stack_isempty(stack) (fmap_vec_size(stack->stack0) == 0 && stack->n_pending == 0)

/*! @function
  @abstract   gets another entry from the memory pool
  @param  mp  pointer to a memory pool
  @return     an element from the memory pool
  @discussion  this will allocate more memory if necessary
  */
inline fmap_map2_entry_p
fmap_map2_mempool_alloc(fmap_map2_mempool_t *mp);

/*! @function
  @abstract   destroys the memory pool
  @param  mp  pointer to a memory pool
  */
void
fmap_map2_mempool_free(fmap_map2_mempool_t *mp, fmap_map2_entry_t *e);

/*! @function
  @abstract   destroys the memory pool
  @param  mp  pointer to a memory pool
  */
void
fmap_map2_mempool_destroy(fmap_map2_mempool_t *mp);

/*! @function
  @abstract  initializes the memory pool
  @return  pointer to a memory pool
  */
fmap_map2_mempool_t *
fmap_map2_mempool_init();

/*! @function
  @abstract  initializes the global memory pool
  @return  pointer to a global memory pool
  */
fmap_map2_global_mempool_t *
fmap_map2_global_mempool_init();

/*! @function
  @abstract       destroys a global memory pool
  @param  global  pointer to a global memory pool
  */
void
fmap_map2_global_mempool_destroy(fmap_map2_global_mempool_t *global);

#endif
