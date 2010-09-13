#ifndef FMAP_BWT_MATCH_H_
#define FMAP_BWT_MATCH_H_

/*! @typedef
  @abstract
  @field  offset  the number of (read) bases used so far in this search (one-based)
  @field  hi      the hash index of the SA interval if the offset is less than or equal to the hash width
  @field  k       the lower occurrence of the SA interval
  @field  l       the upper occurrence of the SA interval
  @discussion    hi is set to UINT32_MAX if the offset is greater than the hash width.  If l is set to
                 UINT32_MAX then that value is unavailable. 
  */
typedef struct {
    uint32_t offset;
    uint32_t hi, k, l;
} fmap_bwt_match_occ_t;

/*! @function
  @abstract        analagous function to fmap_bwt_occ but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  c        base in two-bit integer format
  @param  next     pointer to the next match structure
  @discussion      this will not set the upper occurrence of the SA interval
  */
inline void
fmap_bwt_match_occ(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, uint8_t c, fmap_bwt_match_occ_t *next);

/*! @function
  @abstract        analagous function to fmap_bwt_2occ but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  c        base in two-bit integer format
  @param  next     pointer to the next match structure
  @discussion      this will not set the upper occurrences of the SA interval
  */
inline void
fmap_bwt_match_2occ(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, uint8_t c, fmap_bwt_match_occ_t *next);

/*! @function
  @abstract        analagous function to fmap_bwt_occ4 but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  next     pointer to the next match structure
  */
inline void
fmap_bwt_match_occ4(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, fmap_bwt_match_occ_t *next[4]);

/*! @function
  @abstract        analagous function to fmap_bwt_2occ4 but using a hash
  @param  bwt      pointer to the bwt structure 
  @param  prev     pointer to the previous match structure
  @param  next     pointer to the next match structure
  */
inline void
fmap_bwt_match_2occ4(const fmap_bwt_t *bwt, fmap_bwt_match_occ_t *prev, fmap_bwt_match_occ_t *next[4]);

/*! @typedef
  @abstract   stores the lower bound of the number of mismatches in the string from [i,len-1].
  @field  w    the maximum number of occurrences
  @field  bid  the minimum number of mismatches
  */
typedef struct {
    uint32_t w;
    int32_t bid;
} fmap_bwt_match_width_t;

/*! @function
  @abstract  calculates a lower bound on the number of mismatches in the string for each interval [i,len-1]
  @param  bwt    the reverse BWT being used to search
  @param  len    the length of the string
  @param  str    the string with bases in integer format
  @param  width  array of widths, one for each interval [i,len-1], for 0 <= i < len
*/
void
fmap_bwt_match_cal_width(const fmap_bwt_t *bwt, int len, const char *str, fmap_bwt_match_width_t *width);

/*! @function
  @abstract         computes the SA interval for the given sequence (if any), using forward search
  @param  bwt       pointer to the bwt structure 
  @param  len       the length of the sequence
  @param  str       the DNA sequence in 2-bit format
  @param  sa_begin  the beginning of the returned occurence of the SA interval
  @param  sa_end    the end of the returned occurrence of the SA interval
  @return           the size of the SA interval, 0 if none found
  */
int
fmap_bwt_match_exact(const fmap_bwt_t *bwt, int len, const uint8_t *str, uint32_t *sa_begin, uint32_t *sa_end);

/*! @function
  @abstract         computes the SA interval for the given sequence (if any), using forward search
  @param  bwt       pointer to the bwt structure 
  @param  len       the length of the sequence
  @param  str       the DNA sequence in 2-bit format
  @param  k0        the beginning occurence of the SA interval
  @param  l0        the end occurrence of the SA interval
  @return           the size of the SA interval, 0 if none found
  @discussion       the search will be started at SA interval [k0,l0], with the results returned as [k0,l0]
  */
int
fmap_bwt_match_exact_alt(const fmap_bwt_t *bwt, int len, const uint8_t *str, uint32_t *k0, uint32_t *l0);

#endif
