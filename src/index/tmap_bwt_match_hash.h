/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_BWT_MATCH_HASH_H_
#define TMAP_BWT_MATCH_HASH_H_

/*!
  The BWT hash wrapper and an API for BWT Index Lookups with further hashing
  @details  This API facilitates a secondary hash into the BWT
  */

/*!
  Hash structure
  */
typedef struct {
   void *hash[4]; /*! the hash used by bwt match for each possible next base, the type is defined in the source */ 
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
  @param  c    the base in integer format
  @param  val  the hash value
  @return  0 if the key is present in the hash table; 1 if the bucket is empty (never used); 2 if the element in the bucket has been deleted 
  @details  this does not check if the value is overwritten
 */
int32_t
tmap_bwt_match_hash_put(tmap_bwt_match_hash_t *h, uint32_t key, uint8_t c, uint32_t val);

/*!
  @param  h  the hash 
  @param  key  the hash key
  @param  c    the base in integer format
  @param  found  1 if the it was contained in the hash, 0 otherwise
  @return the hash value, or UINT32_MAX if not found
 */
uint32_t
tmap_bwt_match_hash_get(tmap_bwt_match_hash_t *h, uint32_t key, uint8_t c, uint32_t *found);

/*! 
  analagous function to tmap_bwt_occ but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  c        base in two-bit integer format
  @param  next     pointer to the next match structure
  @param  hash     a occurence array hash
  @details         this will not set the upper occurrence of the SA interval
  */
inline void
tmap_bwt_match_hash_occ(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, uint8_t c, tmap_bwt_match_occ_t *next, tmap_bwt_match_hash_t *hash);

/*! 
  analagous function to tmap_bwt_2occ but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  c        base in two-bit integer format
  @param  next     pointer to the next match structure
  @param  hash     a occurence array hash
  @details         this will not set the upper occurrences of the SA interval
  */
inline void
tmap_bwt_match_hash_2occ(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, uint8_t c, tmap_bwt_match_occ_t *next, tmap_bwt_match_hash_t *hash);

/*! 
  analagous function to tmap_bwt_occ4 but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  next     pointer to the next match structure
  @param  hash     a occurence array hash
  */
inline void
tmap_bwt_match_hash_occ4(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, tmap_bwt_match_occ_t next[4], tmap_bwt_match_hash_t *hash);

/*! 
  analagous function to tmap_bwt_2occ4 but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  next     pointer to the next match structure
  @param  hash     a occurence array hash
  */
inline void
tmap_bwt_match_hash_2occ4(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, tmap_bwt_match_occ_t next[4], tmap_bwt_match_hash_t *hash);

/*! 
  calculates a lower bound on the number of mismatches in the string for each interval [i,len-1]
  @param  bwt    the reverse BWT being used to search
  @param  len    the length of the string
  @param  str    the string with bases in integer format
  @param  width  array of widths, one for each interval [i,len-1], for 0 <= i < len
  @param  hash   a occurence array hash
*/
void
tmap_bwt_match_hash_cal_width_forward(const tmap_bwt_t *bwt, int len, const char *str, tmap_bwt_match_width_t *width, tmap_bwt_match_hash_t *hash);

/*! 
  calculates a lower bound on the number of mismatches in the string for each interval [0,i]
  @param  bwt    the reverse BWT being used to search
  @param  len    the length of the string
  @param  str    the string with bases in integer format
  @param  width  array of widths, one for each interval [0,i], for 0 <= i < len
  @param  hash   a occurence array hash
*/
void
tmap_bwt_match_hash_cal_width_reverse(const tmap_bwt_t *bwt, int len, const char *str, tmap_bwt_match_width_t *width, tmap_bwt_match_hash_t *hash);

/*! 
  computes the SA interval for the given sequence (if any), using forward search
  @param  bwt       pointer to the bwt structure 
  @param  len       the length of the sequence
  @param  str       the DNA sequence in 2-bit format
  @param  match_sa  pointer to the match structure to be returned
  @param  hash      a occurence array hash
  @return           the size of the SA interval, 0 if none found
  */
uint32_t
tmap_bwt_match_hash_exact(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa, tmap_bwt_match_hash_t *hash);

/*! 
  computes the SA interval for the given sequence (if any), using reverse search
  @param  bwt       pointer to the bwt structure 
  @param  len       the length of the sequence
  @param  str       the DNA sequence in 2-bit format
  @param  match_sa  pointer to the match structure to be returned
  @param  hash      a occurence array hash
  @return           the size of the SA interval, 0 if none found
  */
uint32_t
tmap_bwt_match_hash_exact_reverse(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa, tmap_bwt_match_hash_t *hash);

/*! 
  computes the SA interval for the given sequence (if any), using forward search
  @param  bwt       pointer to the bwt structure 
  @param  len       the length of the sequence
  @param  str       the DNA sequence in 2-bit format
  @param  match_sa  pointer to the match structure to be returned
  @param  hash      a occurence array hash
  @return           the size of the SA interval, 0 if none found
  @details          the search will be started at SA interval [k0,l0], with the results returned as [k0,l0]
  */
uint32_t
tmap_bwt_match_hash_exact_alt(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa, tmap_bwt_match_hash_t *hash);

/*! 
  computes the SA interval for the given sequence (if any), using reverse search
  @param  bwt       pointer to the bwt structure 
  @param  len       the length of the sequence
  @param  str       the DNA sequence in 2-bit format
  @param  match_sa  pointer to the match structure to be returned
  @param  hash      a occurence array hash
  @return           the size of the SA interval, 0 if none found
  @details          the search will be started at SA interval [k0,l0], with the results returned as [k0,l0]
  */
uint32_t
tmap_bwt_match_hash_exact_alt_reverse(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa, tmap_bwt_match_hash_t *hash);

/*!  
  inverse Psi function
  @param  bwt  pointer to the bwt structure
  @param  sa_intv  the SA interval
  @param  k  the occurrence position 
  @param  s  the loop counter to return
  @param  hash  a occurence array hash
  @return  the suffix array position
  */
uint32_t
tmap_bwt_match_hash_invPsi(const tmap_bwt_t *bwt, uint32_t sa_intv, uint32_t k, uint32_t *s, tmap_bwt_match_hash_t *hash);

#endif
