#include <stdlib.h>
#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_vec.h"
#include "fmap_map2_mempool.h"

void 
fmap_map2_stack_destroy(fmap_map2_stack_t *stack) 
{ 
  fmap_map2_mempool_destroy(stack->pool);
  fmap_vec_destroy(stack->stack0); 
  fmap_vec_destroy(stack->pending); 
  free(stack); 
}

inline void 
fmap_map2_stack_push0(fmap_map2_stack_t *s, fmap_map2_entry_p e) 
{ 
  fmap_vec_push(fmap_map2_entry_p, s->stack0, e); 
}

inline fmap_map2_entry_p 
fmap_map2_stack_pop(fmap_map2_stack_t *s)
{
  if(fmap_vec_size(s->stack0) == 0 && s->n_pending != 0) {
      fmap_error(NULL, Exit, OutOfRange);
  }
  return fmap_vec_pop(s->stack0);
}

inline fmap_map2_entry_p
fmap_map2_mempool_alloc(fmap_map2_mempool_t *mp)
{
  mp->cnt++;
  if(fmap_vec_size(mp->pool) == 0) {
      return (fmap_map2_entry_t*)fmap_calloc(1, sizeof(fmap_map2_entry_t), "entry");
  }
  else {
      return fmap_vec_pop(mp->pool);
  }
}

void 
fmap_map2_mempool_free(fmap_map2_mempool_t *mp, fmap_map2_entry_t *e)
{
  mp->cnt--;
  e->n = 0;
  fmap_vec_push(fmap_map2_entry_p, mp->pool, e);
}

void 
fmap_map2_mempool_destroy(fmap_map2_mempool_t *mp)
{
  int32_t i;

  for(i=0;i<fmap_vec_size(mp->pool);i++) {
      free(fmap_vec_A(mp->pool, i)->array);
      free(fmap_vec_A(mp->pool, i));
  }
  fmap_vec_destroy(mp->pool);
  free(mp);
}

fmap_map2_global_mempool_t *
fmap_map2_global_mempool_init()
{
  fmap_map2_global_mempool_t *pool= NULL;
  fmap_map2_stack_t *stack = NULL;

  pool = fmap_calloc(1, sizeof(fmap_map2_global_mempool_t), "pool");
  stack = fmap_calloc(1, sizeof(fmap_map2_stack_t), "stack");
  stack->pool = fmap_calloc(1, sizeof(fmap_map2_mempool_t), "stack->pool");
  pool->stack = stack;

  return pool;
}

void
fmap_map2_global_mempool_destroy(fmap_map2_global_mempool_t *global)
{
  fmap_map2_stack_destroy((fmap_map2_stack_t*)global->stack);
  free(global->aln_mem);
  free(global);
}
