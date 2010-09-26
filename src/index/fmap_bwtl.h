#ifndef FMAP_BWT_LITE_H_
#define FMAP_BWT_LITE_H_

#include <stdint.h>

/*! 
  A Lightweight BWT index library; this includes the SA intervals as well.
  */

/*! 
  a light-weight BWT structure
  details  this structure contains the occurrence array, bwt string, and suffix array
  */
typedef struct {
    uint32_t seq_len;  /*!< sequence length */
    uint32_t bwt_size;  /*!< size of the bwt in bytes */
    uint32_t n_occ;  /*!< number of occurences */
    uint32_t primary;  /*!< S^{-1}(0), or the primary index of BWT */
    uint32_t *bwt;  /*!< burrows-wheeler transform */
    uint32_t *occ;  /*!< the occurence array */
    uint32_t *sa;  /*!< the suffix array */
    uint32_t L2[5];  /*!< C(), cumulative count */
    uint32_t cnt_table[256];  /*!< occurrence array */
} fmap_bwtl_t;

/*! 
  @param  b   pointer to the bwt light-weight structure
  @param  k   the zero-based index of the bwt character to retrieve
  @return     the bwt character from the $-removed BWT string.
  details Note that fmap_bwt_t::bwt is not exactly the BWT string 
  and therefore this define is called fmap_bwtl_B0 instead of fmap_bwtl_B. 
  */
#define fmap_bwtl_B0(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

/*! 
  creates a light-weight bwt from the given sequence
  @param  len  the sequence length
  @param  seq  the sequence
  @return      a pointer to the initalized bwt-light-weight structure
  */
fmap_bwtl_t *
fmap_bwtl_seq2bwtl(int32_t len, const uint8_t *seq);

/*! 
  calculates the next occurrence given the previous occurence and the next base
  @param  bwtl  pointer to the bwt structure 
  @param  k     previous occurence
  @param  c     base in two-bit integer format
  @return       the next occurrence given the base
  */
inline uint32_t 
fmap_bwtl_occ(const fmap_bwtl_t *bwtl, uint32_t k, uint8_t c);

/*! 
  calculates the next occurrences given the previous occurence for all four bases
  @param  bwtl  pointer to the bwt structure 
  @param  k     previous occurence
  @param  cnt   pointer to the next occurences for all four bases
  */
inline void 
fmap_bwtl_occ4(const fmap_bwtl_t *bwtl, uint32_t k, uint32_t cnt[4]);

/*! 
  calculates the next SA intervals given the previous SA intervals for all four bases
  @param  bwtl  pointer to the bwt structure 
  @param  k     previous lower occurence
  @param  l     previous upper occurence
  @param  cntk  next upper occurences
  @param  cntl  next lower occurences
  details   more efficient version of bwt_occ4 but requires that k <= l (not checked)
  */
inline void 
fmap_bwtl_2occ4(const fmap_bwtl_t *bwtl, uint32_t k, uint32_t l, uint32_t cntk[4], uint32_t cntl[4]);

/*! 
  @param  bwtl  pointer to the bwt structure 
  */
void 
fmap_bwtl_destroy(fmap_bwtl_t *bwtl);

#endif
