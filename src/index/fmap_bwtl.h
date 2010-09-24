#ifndef FMAP_BWT_LITE_H_
#define FMAP_BWT_LITE_H_

#include <stdint.h>

// TODO: document

/*! @typedef
  @abstract  a light-weight BWT structure
  @field  seq_len   sequence length
  @field  bwt_size  size of the bwt in bytes
  @field  n_occ     number of occurences
  @field  primary   S^{-1}(0), or the primary index of BWT
  @field  bwt       burrows-wheeler transform
  @field  occ       the occurence array
  @field  sa        the suffix array
  @field  L2        C(), cumulative count
  @field  cnt_table occurrence array
  @discussion  this structure contains the occurrence array, bwt string, and suffix array
  */
typedef struct {
    uint32_t seq_len, bwt_size, n_occ;
    uint32_t primary;
    uint32_t *bwt, *occ, *sa, L2[5];
    uint32_t cnt_table[256];
} fmap_bwtl_t;

/*! @macro
  @abstract retrieve a character from the $-removed BWT string. 
  @param  b   pointer to the bwt light-weight structure
  @param  k   the zero-based index of the bwt character to retrieve
  @return     the bwt character from the $-removed BWT string.
  @discussion Note that fmap_bwt_t::bwt is not exactly the BWT string 
  and therefore this macro is called fmap_bwtl_B0 instead of fmap_bwtl_B. 
  */
#define fmap_bwtl_B0(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

/*! @function
  @abstract    creates a light-weight bwt from the given sequence
  @param  len  the sequence length
  @param  seq  the sequence
  @return      a pointer to the initalized bwt-light-weight structure
  */
fmap_bwtl_t *
fmap_bwtl_seq2bwtl(int32_t len, const uint8_t *seq);

/*! @function
  @abstract    calculates the next occurrence given the previous occurence and the next base
  @param  bwt  pointer to the bwt structure 
  @param  k    previous occurence
  @param  c    base in two-bit integer format
  @return      the next occurrence given the base
  */
inline uint32_t 
fmap_bwtl_occ(const fmap_bwtl_t *bwtl, uint32_t k, uint8_t c);

/*! @function
  @abstract    calculates the next occurrences given the previous occurence for all four bases
  @param  bwt  pointer to the bwt structure 
  @param  k    previous occurence
  @param  cnt  pointer to the next occurences for all four bases
  */
inline void 
fmap_bwtl_occ4(const fmap_bwtl_t *bwtl, uint32_t k, uint32_t cnt[4]);

/*! @function
  @abstract     calculates the next SA intervals given the previous SA intervals for all four bases
  @param  bwt   pointer to the bwt structure 
  @param  k     previous lower occurence
  @param  l     previous upper occurence
  @param  cntk  next upper occurences
  @param  cntl  next lower occurences
  @discussion   more efficient version of bwt_occ4 but requires that k <= l (not checked)
  */
inline void 
fmap_bwtl_2occ4(const fmap_bwtl_t *bwtl, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4]);

/*! @function
  @abstract
  @param  bwt  pointer to the bwt-light-weight structure
  */
void 
fmap_bwtl_destroy(fmap_bwtl_t *bwtl);

#endif
