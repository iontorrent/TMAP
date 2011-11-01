/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP2_MEMPOOL_H
#define TMAP_MAP2_MEMPOOL_H

#include <stdint.h>
#include "../../util/tmap_vec.h"
#include "../../index/tmap_bwt_match.h"

/*! 
  Memory Pools for Map2
  */

/*! 
 The Dynamic Programming cell for the DAWG vs. BWT alignment
 */
typedef struct {
    tmap_bwt_match_occ_t match_sa; /*!< lower and upper interval for the query */
    int32_t I;  /*!< insertion score */
    int32_t D;  /*!< deletion score */
    int32_t G;  /*!< maximum of the match/mismatch/insertion/deletion scores */
    uint32_t pj:2;  /*!< the zero-based index of the parent in the prefix trie [0-3] */
    uint32_t qlen:30;  /*!< the number of bases used in the query  */
    int32_t tlen;  /*!< the number of bases used in the target */
    int32_t ppos;  /*!< the parents position in the memory array */
    int32_t upos;  /*!< the child position of its parent in the memory array */
    int32_t cpos[4];  /*!< undetermined */
} tmap_map2_cell_t;

/*! 
 The Dynamic Programming row for the DAWG vs. BWT alignment
*/
typedef struct {
  int32_t n;  /*!< the index of the next cell */
  int32_t max;  /*!< the number of cells allocated */
  uint32_t tk;  /*!< lower suffix array interval of the target */
  uint32_t tl;  /*!< upper suffix array interval of the target  */
  tmap_map2_cell_t *array;  /*!< the array of cells */
} tmap_map2_entry_t, *tmap_map2_entry_p;
/*! 
  a memory pool
*/
typedef struct __tmap_map2_mempool_t {
  int32_t cnt;  /*!< the number of entries being used in the memory pool */
  tmap_vec_t(tmap_map2_entry_p) pool;  /*!< the memory pool vector */
} tmap_map2_mempool_t;
/*! 
  @brief            a two-level memory stack
  @param  n_pending  the number of elements in the pending stack
  @param  stack0     the main stack
  @param  pending    the pending stack
  @param  pool       a memory pool
*/
typedef struct {
  int32_t n_pending;
  tmap_vec_t(tmap_map2_entry_p) stack0, pending;
  struct __tmap_map2_mempool_t *pool;
} tmap_map2_stack_t;
/*! 
  a global memory pool
*/
typedef struct {
  tmap_map2_stack_t *stack;  /*!< the main two-level memory stack */
  int32_t max_l;  /*!< the working memory length */
  uint8_t *aln_mem;  /*!< working memory */
} tmap_map2_global_mempool_t;
/*! 
  destroys the stack
  @param  stack  a pointer to the stack
  */
void 
tmap_map2_stack_destroy(tmap_map2_stack_t *stack);

/*! 
  push an entry onto the main stack
  @param  stack  a pointer to the stack
  @param  e      the element to push
  */
inline void
tmap_map2_stack_push0(tmap_map2_stack_t *stack, tmap_map2_entry_p e);

/*! 
  pop an entry from the main stack
  @param  stack  a pointer to the stack
  @return        the popped element
  @details       this will fail (error) if the stack is empty
  */
inline tmap_map2_entry_p 
tmap_map2_stack_pop(tmap_map2_stack_t *stack);

/*! 
  checks if is the stack empty
  @param  stack  a pointer to the stack
  @return        true if the stack is empty, false otherwise 
  */
#define tmap_map2_stack_isempty(stack) (tmap_vec_size(stack->stack0) == 0 && stack->n_pending == 0)

/*! 
  gets another entry from the memory pool
  @param  mp  pointer to a memory pool
  @return     an element from the memory pool
  @details    this will allocate more memory if necessary
  */
inline tmap_map2_entry_p
tmap_map2_mempool_pop(tmap_map2_mempool_t *mp);

/*! 
  relinquishes the entry back to the memory pool
  @param  mp  pointer to a memory pool
  @param  e   pointer to the entry to relinquish
  */
void
tmap_map2_mempool_push(tmap_map2_mempool_t *mp, tmap_map2_entry_t *e);

/*! 
  destroys the memory pool
  @param  mp  pointer to a memory pool
  */
void
tmap_map2_mempool_destroy(tmap_map2_mempool_t *mp);

/*! 
  initializes the memory pool
  @return  pointer to a memory pool
  */
tmap_map2_mempool_t *
tmap_map2_mempool_init();

/*! 
  initializes the global memory pool
  @return  pointer to a global memory pool
  */
tmap_map2_global_mempool_t *
tmap_map2_global_mempool_init();

/*! 
  destroys a global memory pool
  @param  global  pointer to a global memory pool
  */
void
tmap_map2_global_mempool_destroy(tmap_map2_global_mempool_t *global);

#endif
