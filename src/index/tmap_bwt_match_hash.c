#include <stdio.h>
#include <stdlib.h>
#include "../util/tmap_hash.h"
#include "tmap_bwt_match_hash.h"

// TODO: change this if the size of the pacpos
TMAP_HASH_MAP_INIT_INT(tmap_bwt_match_hash, uint32_t) 

tmap_bwt_match_hash_t*
tmap_bwt_match_hash_init()
{
  tmap_bwt_match_hash_t *h = NULL;
  h = tmap_calloc(sizeof(tmap_bwt_match_hash_t), 1, "hash");
  h->hash = (void*)tmap_hash_init(tmap_bwt_match_hash);
  return h;
}

void
tmap_bwt_match_hash_destroy(tmap_bwt_match_hash_t *h)
{
  tmap_hash_t(tmap_bwt_match_hash) *hash = (tmap_hash_t(tmap_bwt_match_hash)*)(h->hash);
  tmap_hash_destroy(tmap_bwt_match_hash, hash); 
  free(h);
}

void
tmap_bwt_match_hash_clear(tmap_bwt_match_hash_t *h)
{
  tmap_hash_t(tmap_bwt_match_hash) *hash = (tmap_hash_t(tmap_bwt_match_hash)*)(h->hash);
  tmap_hash_clear(tmap_bwt_match_hash, hash); 
  free(h);
}

int32_t
tmap_bwt_match_hash_put(tmap_bwt_match_hash_t *h, uint32_t key, uint32_t val)
{
  int32_t ret;
  tmap_hash_int_t iter;
  tmap_hash_t(tmap_bwt_match_hash) *hash = (tmap_hash_t(tmap_bwt_match_hash)*)(h->hash);
  // get the iter
  iter = tmap_hash_put(tmap_bwt_match_hash, hash, key, &ret); 
  // set the value
  tmap_hash_val(hash, iter) = val;
  return ret;
}

uint32_t
tmap_bwt_match_hash_get(tmap_bwt_match_hash_t *h, uint32_t key)
{
  tmap_hash_int_t iter;
  tmap_hash_t(tmap_bwt_match_hash) *hash = (tmap_hash_t(tmap_bwt_match_hash)*)(h->hash);
  // get the iter
  iter = tmap_hash_get(tmap_bwt_match_hash, hash, key); 
  // check it was found
  if(tmap_hash_end(hash) == iter) {
      return UINT32_MAX;
  }
  // return the value
  return (uint32_t)(tmap_hash_val(hash, iter));
}
