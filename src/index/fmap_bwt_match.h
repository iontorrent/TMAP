#ifndef FMAP_BWT_MATCH_H_
#define FMAP_BWT_MATCH_H_

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
