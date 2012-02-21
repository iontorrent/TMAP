/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SA_H
#define TMAP_SA_H

#define TMAP_SA_INTERVAL 32

#include <stdint.h>
#include "tmap_bwt.h"
#include "tmap_bwt_match.h"
#include "tmap_bwt_match_hash.h"

/*! 
  A Suffix Array Library
  */

/*! 
  */
typedef struct {
    tmap_bwt_int_t primary;  /*!< S^{-1}(0), or the primary index of BWT */
    uint32_t sa_intv;  /*!< the suffix array interval (sampled) */
    tmap_bwt_int_t seq_len;  /*!< the length of the reference sequence */
    tmap_bwt_int_t n_sa;  /*!< number of suffix array entries */
    tmap_bwt_int_t *sa;  /*!< pointer to the suffix array entries */
    uint32_t is_shm;  /*!< 1 if loaded from shared memory, 0 otherwise */
    // Not stored in the file
    uint32_t sa_intv_log2;  /*!< the log2 suffix array interval (sampled) */
} tmap_sa_t;

/*! 
  @param  fn_fasta  the FASTA file name
  @return           pointer to the sa structure 
  */
tmap_sa_t *
tmap_sa_read(const char *fn_fasta);

/*! 
  @param  fn_fasta  the FASTA file name
  @param  sa        the sa structure to write
  */
void
tmap_sa_write(const char *fn_fasta, tmap_sa_t *sa);

/*! 
  @param  len      the reference length
  @param  sa_intv  the suffix array interval
  @return          the approximate number of bytes required for this sa in shared memory
  */
size_t
tmap_sa_approx_num_bytes(uint64_t len, uint32_t sa_intv);

/*! 
  @param  sa  the sa structure 
  @return     the number of bytes required for this sa in shared memory
  */
size_t
tmap_sa_shm_num_bytes(tmap_sa_t *sa);

/*! 
  @param  fn_fasta  the FASTA file name
  @return     the number of bytes required for this sa in shared memory
  */
size_t
tmap_sa_shm_read_num_bytes(const char *fn_fasta);

/*! 
  @param  sa  the sa structure to pack 
  @param  buf  the byte array in which to pack the sa data
  @return      a pointer to the next unused byte in memory
  */
uint8_t *
tmap_sa_shm_pack(tmap_sa_t *sa, uint8_t *buf);

/*! 
  @param  buf  the byte array in which to unpack the sa data
  @return      a pointer to the initialized sa structure
  */
tmap_sa_t *
tmap_sa_shm_unpack(uint8_t *buf);

/*! 
  @param  sa  pointer to the suffix array structure
  */
void 
tmap_sa_destroy(tmap_sa_t *sa);

/*! 
  returns the suffix array position given the occurrence position
  @param  sa   the suffix array
  @param  bwt  the bwt structure 
  @param  k    the suffix array position
  @param  hash  the user occurence hash
  @return      the pac position
*/
tmap_bwt_int_t
tmap_sa_pac_pos_hash(const tmap_sa_t *sa, const tmap_bwt_t *bwt, tmap_bwt_int_t k, tmap_bwt_match_hash_t *hash);

/*! 
  returns the suffix array position given the occurrence position
  @param  sa   the suffix array
  @param  bwt  the bwt structure 
  @param  k    the suffix array position
  @return      the pac position
*/
tmap_bwt_int_t
tmap_sa_pac_pos(const tmap_sa_t *sa, const tmap_bwt_t *bwt, tmap_bwt_int_t k);

/*! 
  @param  fn_fasta  the FASTA file name
  @param  intv      the suffix array interval
  */
void
tmap_sa_bwt2sa(const char *fn_fasta, uint32_t intv);

/*! 
  constructs the suffix array of a given string.
  @param T   T[0..n-1] The input string.
  @param SA  SA[0..n] The output array of suffixes.
  @param n   the length of the given string.
  @return    0 if no error occurred
 */
uint32_t 
tmap_sa_gen_short(const uint8_t *T, int32_t *SA, uint32_t n);

/*! 
  main-like function for 'tmap bwt2sa'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
tmap_sa_bwt2sa_main(int argc, char *argv[]);

#define KEY(V, I, p, h)                                 ( V[ I[p] + h ] )
#define INSERT_SORT_NUM_ITEM    16

void QSufSortSuffixSort(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t numChar, const tmap_bwt_sint_t largestInputSymbol,
                                                                        const tmap_bwt_sint_t smallestInputSymbol, const int skipTransform);
void QSufSortGenerateSaFromInverse(const tmap_bwt_sint_t *V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t numChar);

#endif
