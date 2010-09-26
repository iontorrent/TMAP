#ifndef FMAP_BWT_MATCH_H_
#define FMAP_BWT_MATCH_H_

/*! 
  API for BWT Index Lookups
  @details  This API facilitates a secondary hash into the BWT
  */

/*! 
  @details hi is set to UINT32_MAX if the offset is greater than the hash width.  If l is set to
  UINT32_MAX then that value is unavailable. 
  */
typedef struct {
    uint32_t offset;  /*!< the number of (read) bases used so far in this search (one-based) */
    uint32_t hi;  /*!< the hash index of the SA interval if the offset is less than or equal to the hash width */
    uint32_t k;  /*!< the lower occurrence of the SA interval */
    uint32_t l;  /*!< the upper occurrence of the SA interval */
} fmap_bwt_match_occ_t;

/*! 
  analagous function to fmap_bwt_occ but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  c        base in two-bit integer format
  @param  next     pointer to the next match structure
  @details         this will not set the upper occurrence of the SA interval
  */
inline void
fmap_bwt_match_occ(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, uint8_t c, fmap_bwt_match_occ_t *next);

/*! 
  analagous function to fmap_bwt_2occ but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  c        base in two-bit integer format
  @param  next     pointer to the next match structure
  @details         this will not set the upper occurrences of the SA interval
  */
inline void
fmap_bwt_match_2occ(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, uint8_t c, fmap_bwt_match_occ_t *next);

/*! 
  analagous function to fmap_bwt_occ4 but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  next     pointer to the next match structure
  */
inline void
fmap_bwt_match_occ4(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, fmap_bwt_match_occ_t next[4]);

/*! 
  analagous function to fmap_bwt_2occ4 but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  next     pointer to the next match structure
  */
inline void
fmap_bwt_match_2occ4(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, fmap_bwt_match_occ_t next[4]);

/*! 
  stores the lower bound of the number of mismatches in the string from [i,len-1].
  */
typedef struct {
    uint32_t w;  /*!< the maximum number of occurrences */
    int32_t bid;  /*!< the minimum number of mismatches */
} fmap_bwt_match_width_t;

/*! 
  calculates a lower bound on the number of mismatches in the string for each interval [i,len-1]
  @param  bwt    the reverse BWT being used to search
  @param  len    the length of the string
  @param  str    the string with bases in integer format
  @param  width  array of widths, one for each interval [i,len-1], for 0 <= i < len
*/
void
fmap_bwt_match_cal_width(const fmap_bwt_t *bwt, int len, const char *str, fmap_bwt_match_width_t *width);
/*! 
  computes the SA interval for the given sequence (if any), using forward search
  @param  bwt       pointer to the bwt structure 
  @param  len       the length of the sequence
  @param  str       the DNA sequence in 2-bit format
  @param  match_sa  pointer to the match structure to be returned
  @return           the size of the SA interval, 0 if none found
  */
int
fmap_bwt_match_exact(const fmap_bwt_t *bwt, int len, const uint8_t *str, fmap_bwt_match_occ_t *match_sa);

/*! 
  computes the SA interval for the given sequence (if any), using forward search
  @param  bwt       pointer to the bwt structure 
  @param  len       the length of the sequence
  @param  str       the DNA sequence in 2-bit format
  @param  match_sa  pointer to the match structure to be returned
  @return           the size of the SA interval, 0 if none found
  @details          the search will be started at SA interval [k0,l0], with the results returned as [k0,l0]
  */
int
fmap_bwt_match_exact_alt(const fmap_bwt_t *bwt, int len, const uint8_t *str, fmap_bwt_match_occ_t *match_sa);

#endif
