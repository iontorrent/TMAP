#ifndef FMAP_BWT_MATCH_H_
#define FMAP_BWT_MATCH_H_

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
