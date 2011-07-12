/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/* The MIT License

   Copyright (c) 2003-2006, 2008, by Heng Li <lh3lh3@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING TMAP_SW_FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   */

#ifndef TMAP_SW_H
#define TMAP_SW_H

#include <stdint.h>

#define TMAP_SW_CIGAR_OP(_cigar) (((_cigar) & 0xf))
#define TMAP_SW_CIGAR_LENGTH(_cigar) (((_cigar) >> 4))
#define TMAP_SW_CIGAR_STORE(_cigar, _op, _len) ((_cigar) = ((_len) << 4) | ((_op) & 0xf))
#define TMAP_SW_CIGAR_ADD_LENGTH(_cigar, _add) ((_cigar) += ((_add) << 4))

/*! 
  Functions for Performing Efficient Smith-Waterman
  */

/*!
  From which cell; helps recovert the best scoring path.
  */
enum {
    TMAP_SW_FROM_M = 0, /*!< from a mismatch cell*/
    TMAP_SW_FROM_I = 1, /*!< from an insertion cell */
    TMAP_SW_FROM_D = 2, /*!< from a deletion cell */
    TMAP_SW_FROM_S = 3  /*!< from a start cell */
};

/*! This is the smallest integer for Smith-Waterman alignment. It might be CPU-dependent in very RARE cases. */
#define TMAP_SW_MINOR_INF INT32_MIN/2
//#define TMAP_SW_MINOR_INF -1073741823

#define TMAP_SW_SET_INF(s) (s).match_score = (s).ins_score = (s).del_score = TMAP_SW_MINOR_INF

#define TMAP_SW_SET_FROM(s, from) (s).match_from = (s).ins_from = (s).del_from = from 

/*!
  Stores from which cell the current cell was extended
 */
typedef struct
{
    uint8_t match_from:3; /*!< from cell for match */
    uint8_t ins_from:2; /*!< from cell for insertion */
    uint8_t del_from:2; /*!< from cell for deletion */
} tmap_sw_dpcell_t;

/*!
  Stores the score for the current cell
  */
typedef struct
{
    int32_t match_score; /*!< match score */
    int32_t ins_score; /*!< insertion score */
    int32_t del_score; /*! <deletion score */
} tmap_sw_dpscore_t;

/*!
  Parameters for the Smith-Waterman alignment.
  @details  The gap penalties should be positive (they will be subtracted)
  while the substitution score should be positive for matches, and negative for
  mismatches (they will be added).
  */
typedef struct
{
  int32_t gap_open; /*!< gap open penalty (positive) */
  int32_t gap_ext; /*!< gap extension penalty (positive) */
  int32_t gap_end; /*!< gap end penalty (positive */
  int32_t *matrix; /*!< substitution matrix (see details)*/
  int32_t row; /*!< the alphabet size */  
  int32_t band_width; /*!< for Smith-Waterman banding */
} tmap_sw_param_t;

/*!
  The best-scoring alignment path.
  */
typedef struct
{
  int32_t i; /*!< the seq1 index (1-based) */
  int32_t j; /*!< the seq2 index (1-based) */
  uint8_t ctype; /*!< the edit operator applied */
} tmap_sw_path_t;

/*!
  Stores a Smith Waterman alignment.
  */
typedef struct
{
  tmap_sw_path_t *path; /*!< Smith-Waterman alignment path */
  int32_t path_len; /*!< The path length */
  int32_t start1; /*!< the start of the first sequence (1-based) */
  int32_t end1; /*!< the end of the first sequence (1-based) */
  int32_t start2; /*!< the start of the second sequence (1-based) */
  int32_t end2; /*!< the end of the second sequence (1-based) */
  int32_t score; /*!< the alignment score */
  int32_t subo; /*!< the sub-optimal alignment score (next best) */
  char *out1; /*!< sequence 1 alignment string */ 
  char *out2; /*!< sequence 2 alignment string */ 
  char *outm; /*!< match/indel alignment string */
  int32_t n_cigar; /*!< the number of cigar operators */
  uint32_t *cigar32; /*!< the cigar operators */
} tmap_sw_aln_t;

/*!
  Initialize an alignment
  @return  a pointer to the initialized memory
  */
tmap_sw_aln_t *
tmap_sw_aln_init();

/*!
  Destroy an alignment
  @param  aa  pointer to the alignment
  */
void
tmap_sw_aln_destroy(tmap_sw_aln_t *aa);

/*!
  Performs the global Smith-Waterman alignment.
  @details          actually, it performs it with banding.
  @param  seq1      the first DNA sequence (in 2-bit format)
  @param  len1      the length of the first sequence
  @param  seq2      the second DNA sequence (in 2-bit format)
  @param  len2      the length of the second sequence
  @param  ap        the alignment parameters
  @param  path      the Smith-Waterman alignment path
  @param  path_len  the Smith-Waterman alignment path length
  @param  right_j   0 if we are to left-justify indels, 1 otherwise
  @return           the alignment score, 0 if none was found
  */
int32_t 
tmap_sw_global_core(uint8_t *seq1, int32_t len1, 
                    uint8_t *seq2, int32_t len2, 
                    const tmap_sw_param_t *ap,
                    tmap_sw_path_t *path, int32_t *path_len,
                    int32_t right_j);

/*!
  Performs the local Smith-Waterman alignment.
  @details          actually, it performs it with banding.
  @param  seq1      the first DNA sequence (in 2-bit format)
  @param  len1      the length of the first sequence
  @param  seq2      the second DNA sequence (in 2-bit format)
  @param  len2      the length of the second sequence
  @param  ap        the alignment parameters
  @param  path      the Smith-Waterman alignment path
  @param  path_len  the Smith-Waterman alignment path length
  @param  _thres    the scoring threshold for local alignment only (the absolute value will be taken); a value zero or negative value will cause no path to be filled
  @param  _subo     the sub-optimal alignment score (next best) 
  @return           the alignment score, 0 if none was found
  */
int32_t 
tmap_sw_local_core(uint8_t *seq1, int32_t len1, 
                   uint8_t *seq2, int32_t len2, 
                   const tmap_sw_param_t *ap,
                   tmap_sw_path_t *path, int32_t *path_len, 
                   int32_t _thres, int32_t *_subo);

/*!
  Performs banded global Smith-Waterman alignment.
  @details          the band width is determined by the expected alignment score 
  @param  seq1      the first DNA sequence (in 2-bit format)
  @param  len1      the length of the first sequence
  @param  seq2      the second DNA sequence (in 2-bit format)
  @param  len2      the length of the second sequence
  @param  ap        the alignment parameters
  @param  path      the Smith-Waterman alignment path
  @param  path_len  the Smith-Waterman alignment path length
  @param  right_j   0 if we are to left-justify indels, 1 otherwise
  @return           the alignment score, 0 if none was found
  */
int32_t
tmap_sw_global_banded_core(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
                           int32_t score, tmap_sw_path_t *path, int32_t *path_len, int32_t right_j);

/*!
  Extens an alignment with the local Smith-Waterman.
  @details            actually, it performs it with banding.
  @param  seq1        the first DNA sequence (in 2-bit format)
  @param  len1        the length of the first sequence
  @param  seq2        the second DNA sequence (in 2-bit format)
  @param  len2        the length of the second sequence
  @param  ap          the alignment parameters
  @param  path        the Smith-Waterman alignment path
  @param  path_len    the Smith-Waterman alignment path length
  @param  prev_score  the initial alignment score
  @param  _mem        allocated memory with size of (len1+2)*(ap->row+1)*4
  @param  right_j     0 if we are to left-justify indels, 1 otherwise
  @return             the alignment score, 0 if none was found
  */
int32_t 
tmap_sw_extend_core(uint8_t *seq1, int32_t len1, 
                    uint8_t *seq2, int32_t len2, 
                    const tmap_sw_param_t *ap,
                    tmap_sw_path_t *path, int32_t *path_len, 
                    int32_t prev_score, uint8_t *_mem, 
                    int32_t right_j);

/*!
  Extens an alignment with the local Smith-Waterman and aligns the entire seq2.
  @param  seq1        the first DNA sequence (in 2-bit format)
  @param  len1        the length of the first sequence
  @param  seq2        the second DNA sequence (in 2-bit format)
  @param  len2        the length of the second sequence
  @param  ap          the alignment parameters
  @param  path        the Smith-Waterman alignment path
  @param  path_len    the Smith-Waterman alignment path length
  @param  prev_score  the initial alignment score
  @param  _mem        allocated memory with size of (len1+2)*(ap->row+1)*4
  @param  right_j     0 if we are to left-justify indels, 1 otherwise
  @return             the alignment score, 0 if none was found
  @details            actually, it performs it with banding.
  */
int32_t 
tmap_sw_extend_fitting_core(uint8_t *seq1, int32_t len1, 
                    uint8_t *seq2, int32_t len2, 
                    const tmap_sw_param_t *ap,
                    tmap_sw_path_t *path, int32_t *path_len, 
                    int32_t prev_score, uint8_t *_mem,
                    int32_t right_j);

/*!
  Performs a fitting aligment, whereby seq2 is fit into seq1
  @param  seq1             the first DNA sequence (in 2-bit format)
  @param  len1             the length of the first sequence
  @param  seq2             the second DNA sequence (in 2-bit format)
  @param  len2             the length of the second sequence
  @param  ap               the alignment parameters
  @param  seq2_start_clip  1 to allow clipping at the start of seq2, 0 otherwise
  @param  seq2_end_clip    1 to allow clipping at the end of seq2, 0 otherwise
  @param  path             the Smith-Waterman alignment path
  @param  path_len         the Smith-Waterman alignment path length
  @param  right_j          0 if we are to left-justify indels, 1 otherwise
  @return                  the alignment score, 0 if none was found
  */
int32_t
tmap_sw_clipping_core(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
                      int32_t seq2_start_clip, int32_t seq2_end_clip,
                      tmap_sw_path_t *path, int32_t *path_len,
                      int32_t right_j);

/*!
  Performs a fitting aligment, whereby seq2 is fit into seq1
  @param  seq1      the first DNA sequence (in 2-bit format)
  @param  len1      the length of the first sequence
  @param  seq2      the second DNA sequence (in 2-bit format)
  @param  len2      the length of the second sequence
  @param  ap        the alignment parameters
  @param  path      the Smith-Waterman alignment path
  @param  path_len  the Smith-Waterman alignment path length
  @param  right_j   0 if we are to left-justify indels, 1 otherwise
  @return           the alignment score, 0 if none was found
  */
int32_t
tmap_sw_fitting_core(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
                                tmap_sw_path_t *path, int32_t *path_len, int32_t right_j);

/*!
  Creates a cigar array from an alignment path
  @param  path      the Smith-Waterman alignment path
  @param  path_len  the Smith-Waterman alignment path length
  @param  n_cigar   pointer to the returned number of cigar operations
  @return           the cigar array, NULL if the path is NULL or the path length is zero 
  */
uint32_t *
tmap_sw_path2cigar(const tmap_sw_path_t *path, int32_t path_len, int32_t *n_cigar);


/********************
 * global variables *
 ********************/

/*!
  alignment parameters for short read alignment [ACGTN]
  */
extern tmap_sw_param_t tmap_sw_param_short; /* = { 13,  2,  2, tmap_sw_sm_short, 5, 50 }; */
/*!
  alignment parameters for blast read alignment [ACGTN]
  */
extern tmap_sw_param_t tmap_sw_param_blast; /* = {  5,  2,  2, tmap_sw_sm_blast, 5, 50 }; */
/*!
  alignment parameters for ... 
  */
extern tmap_sw_param_t tmap_sw_param_nt2nt; /* = {  8,  2,  2, tmap_sw_sm_nt, 16, 75 }; */
/*!
  alignment parameters for ...
  */
extern tmap_sw_param_t tmap_sw_param_rd2rd; /* = {  1, 19, 19, tmap_sw_sm_read, 16, 75 }; */
/*!
  alignment parameters for ...
  */
extern tmap_sw_param_t tmap_sw_param_aa2aa; /* = { 10,  2,  2, tmap_sw_sm_blosum62, 22, 50 }; */

/*!
  substitution matrix for short read alignment 
  */
extern int32_t tmap_sw_sm_short[];

/*!
  substitution matrix for blast
  */
extern int32_t tmap_sw_sm_blast[];

/*!
  substitution matrix for ...
  */
extern int32_t tmap_sw_sm_nt[];

/*!
  substitution matrix for ...
  */
extern int32_t tmap_sw_sm_read[];

/*!
  substitution matrix for BLOSUM62
  */
extern int32_t tmap_sw_sm_blosum62[];

/*!
  substitution matrix for BLOSUM45
  */
extern int32_t tmap_sw_sm_blosum45[];

/*!
  substitution matrix for human to mouse
  */
extern int32_t           tmap_sw_sm_hs[];

#endif
