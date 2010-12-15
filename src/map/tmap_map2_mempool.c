#include <stdlib.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_vec.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_bwtl.h"
#include "../index/tmap_sa.h"
#include "tmap_map_util.h"
#include "tmap_map2.h"
#include "tmap_map2_aux.h"
#include "tmap_map2_core.h"
#include "tmap_map2_mempool.h"

void 
tmap_map2_stack_destroy(tmap_map2_stack_t *stack) 
{ 
  tmap_map2_mempool_destroy(stack->pool);
  tmap_vec_destroy(stack->stack0); 
  tmap_vec_destroy(stack->pending); 
  free(stack); 
}

inline void 
tmap_map2_stack_push0(tmap_map2_stack_t *s, tmap_map2_entry_p e) 
{ 
  tmap_vec_push(tmap_map2_entry_p, s->stack0, e); 
}

inline tmap_map2_entry_p 
tmap_map2_stack_pop(tmap_map2_stack_t *s)
{
  if(tmap_vec_size(s->stack0) == 0 && s->n_pending != 0) {
      tmap_error(NULL, Exit, OutOfRange);
  }
  return tmap_vec_pop(s->stack0);
}

inline tmap_map2_entry_p
tmap_map2_mempool_pop(tmap_map2_mempool_t *mp)
{
  mp->cnt++;
  if(tmap_vec_size(mp->pool) == 0) {
      return (tmap_map2_entry_t*)tmap_calloc(1, sizeof(tmap_map2_entry_t), "entry");
  }
  else {
      return tmap_vec_pop(mp->pool);
  }
}

void 
tmap_map2_mempool_push(tmap_map2_mempool_t *mp, tmap_map2_entry_t *e)
{
  mp->cnt--;
  e->n = 0;
  tmap_vec_push(tmap_map2_entry_p, mp->pool, e);
}

void 
tmap_map2_mempool_destroy(tmap_map2_mempool_t *mp)
{
  int32_t i;

  if(0 != mp->cnt) tmap_error("memory leak detected", Exit, OutOfRange);

  for(i=0;i<tmap_vec_size(mp->pool);i++) {
      free(tmap_vec_A(mp->pool, i)->array);
      free(tmap_vec_A(mp->pool, i));
  }
  tmap_vec_destroy(mp->pool);
  free(mp);
}

tmap_map2_global_mempool_t *
tmap_map2_global_mempool_init()
{
  tmap_map2_global_mempool_t *pool= NULL;

  pool = tmap_calloc(1, sizeof(tmap_map2_global_mempool_t), "pool");
  pool->stack = tmap_calloc(1, sizeof(tmap_map2_stack_t), "stack");
  pool->stack->pool = tmap_calloc(1, sizeof(tmap_map2_mempool_t), "stack->pool");

  return pool;
}

void
tmap_map2_global_mempool_destroy(tmap_map2_global_mempool_t *global)
{
  tmap_map2_stack_destroy((tmap_map2_stack_t*)global->stack);
  free(global->aln_mem);
  free(global);
}
