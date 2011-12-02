#ifndef TMAP_BWT_MATCH_HASH_H_
#define TMAP_BWT_MATCH_HASH_H_

/*!
  The BWT hash wrapper
 */
typedef struct {
   void *hash; /*! the hash used by bwt match, the type is defined in the source */ 
} tmap_bwt_match_hash_t;

/*!
  @return  the initialized hash structure
 */
tmap_bwt_match_hash_t*
tmap_bwt_match_hash_init();

/*!
  @param  h  the hash to destroy
 */
void
tmap_bwt_match_hash_destroy(tmap_bwt_match_hash_t *h);

/*!
  @param  h  the hash to clear
 */
void
tmap_bwt_match_hash_clear(tmap_bwt_match_hash_t *h);

/*!
  @param  h  the hash 
  @param  key  the hash key
  @param  val  the hash value
  @return  0 if the key is present in the hash table; 1 if the bucket is empty (never used); 2 if the element in the bucket has been deleted 
  @details  this does not check if the value is overwritten
 */
int32_t
tmap_bwt_match_hash_put(tmap_bwt_match_hash_t *h, uint32_t key, uint32_t val);

/*!
  @param  h  the hash 
  @param  key  the hash key
  @return the hash value, or UINT32_MAX if not found
 */
uint32_t
tmap_bwt_match_hash_get(tmap_bwt_match_hash_t *h, uint32_t);

#endif
