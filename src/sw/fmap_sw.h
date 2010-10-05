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
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FMAP_SW_FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   */

#ifndef FMAP_SW_H_
#define FMAP_SW_H_

#include <stdint.h>

/*! 
  Functions for Performing Efficient Smith-Waterman
  */

/*!
  From which cell; helps recovert the best scoring path.
  */
enum {
    FMAP_SW_FROM_M = 0, /*!< from a mismatch cell*/
    FMAP_SW_FROM_I = 1, /*!< from an insertion cell */
    FMAP_SW_FROM_D = 2, /*!< from a deletion cell */
    FMAP_SW_FROM_S = 3  /*!< from a start cell */
};

/*!
  The type of Smith-Waterman to perform.
  */
enum {
    FMAP_SW_TYPE_LOCAL  = 0, /*!< local alignment */
    FMAP_SW_TYPE_GLOBAL = 1, /*!< global alignment */
    FMAP_SW_TYPE_EXTEND = 2  /*!< extend an alignment */
};

/*! This is the smallest integer for Smith-Waterman alignment. It might be CPU-dependent in very RARE cases. */
#define FMAP_SW_MINOR_INF INT32_MIN/2
//#define FMAP_SW_MINOR_INF -1073741823

#define FMAP_SW_SET_INF(s) (s).match_score = (s).ins_score = (s).del_score = FMAP_SW_MINOR_INF


/*!
  Stores from which cell the current cell was extended
 */
typedef struct
{
    uint8_t match_from:3; /*!< from cell for match */
    uint8_t ins_from:2; /*!< from cell for insertion */
    uint8_t del_from:2; /*!< from cell for deletion */
} fmap_sw_dpcell_t;

/*!
  Stores the score for the current cell
  */
typedef struct
{
    int32_t match_score; /*!< match score */
    int32_t ins_score; /*!< insertion score */
    int32_t del_score; /*! <deletion score */
} fmap_sw_dpscore_t;

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
} fmap_sw_param_t;

/*!
  The best-scoring alignment path.
  */
typedef struct
{
  int32_t i; /*!< the seq1 index (1-based) */
  int32_t j; /*!< the seq2 index (1-based) */
  uint8_t ctype; /*!< the edit operator applied */
} fmap_sw_path_t;

/*!
  Stores a Smith Waterman alignment.
  */
typedef struct
{
  fmap_sw_path_t *path; /*!< Smith-Waterman alignment path */
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
} fmap_sw_aln_t;

/*!
  Initialize an alignment
  @return  a pointer to the initialized memory
  */
fmap_sw_aln_t *
fmap_sw_aln_init();

/*!
  Destroy an alignment
  @param  aa  pointer to the alignment
  */
void
fmap_sw_aln_destroy(fmap_sw_aln_t *aa);

/*!
  Performs Smith-Waterman alignment.
  @details  actually, it performs banded Smith-Waterman.
  @param  seq1   the first DNA sequence (character array)
  @param  seq2   the second DNA sequence (character array)
  @param  ap     the alignment parameters
  @param  type   the Smith-Waterman type (global, local, extend)
  @param  thres  the scoring threshold for local alignment only (the absolute value will be taken); a value zero or negative value will cause no path to be filled
  @param  len1   the length of the first sequence to align, or -1 to align the full sequence
  @param  len2   the length of the second sequence to align, or -1 to align the full sequence
  @return        the Smith-Waterman alignment
  */
fmap_sw_aln_t *
fmap_sw_stdaln_aux(const char *seq1, const char *seq2, const fmap_sw_param_t *ap,
                   int32_t type, int32_t thres , int32_t len1, int32_t len2);

/*!
  Performs Smith-Waterman alignment.
  @details  actually, it performs banded Smith-Waterman.
  @param  seq1   the first DNA sequence (in 2-bit format)
  @param  seq2   the second DNA sequence (in 2-bit format)
  @param  ap     the alignment parameters
  @param  type   the Smith-Waterman type (global, local, extend)
  @param  thres  the scoring threshold for local alignment only (the absolute value will be taken); a value zero or negative value will cause no path to be filled
  @return        the Smith-Waterman alignment
  */
fmap_sw_aln_t *
fmap_sw_stdaln(const char *seq1, const char *seq2, const fmap_sw_param_t *ap, int32_t type, int32_t thres);

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
  @return           the alignment score, 0 if none was found
  */
int32_t 
fmap_sw_global_core(uint8_t *seq1, int32_t len1, 
                    uint8_t *seq2, int32_t len2, 
                    const fmap_sw_param_t *ap,
                    fmap_sw_path_t *path, int32_t *path_len);

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
  @param  seq1_fit   1 if the whole of seq1 should be aligned, 0 otherwise
  @param  seq2_fit   1 if the whole of seq2 should be aligned, 0 otherwise
  @param  _thres    the scoring threshold for local alignment only (the absolute value will be taken); a value zero or negative value will cause no path to be filled
  @param  _subo     the sub-optimal alignment score (next best) 
  @return           the alignment score, 0 if none was found
  */
int32_t 
fmap_sw_local_core(uint8_t *seq1, int32_t len1, 
                   uint8_t *seq2, int32_t len2, 
                   const fmap_sw_param_t *ap,
                   fmap_sw_path_t *path, int32_t *path_len, 
                   int32_t seq1_fit, int32_t seq2_fit,
                   int32_t _thres, int32_t *_subo);

/*!
  Extens an alignment with the local Smith-Waterman.
  @details          actually, it performs it with banding.
  @param  seq1      the first DNA sequence (in 2-bit format)
  @param  len1      the length of the first sequence
  @param  seq2      the second DNA sequence (in 2-bit format)
  @param  len2      the length of the second sequence
  @param  ap        the alignment parameters
  @param  path      the Smith-Waterman alignment path
  @param  path_len  the Smith-Waterman alignment path length
  @param  seq1_fit   1 if the whole of seq1 should be aligned, 0 otherwise
  @param  seq2_fit   1 if the whole of seq2 should be aligned, 0 otherwise
  @param  G0        the initial alignment score
  @param  _mem      allocated memory with size of (len1+2)*(ap->row+1)*4
  @return           the alignment score, 0 if none was found
  */
int32_t 
fmap_sw_extend_core(uint8_t *seq1, int32_t len1, 
                    uint8_t *seq2, int32_t len2, 
                    const fmap_sw_param_t *ap,
                    fmap_sw_path_t *path, int32_t *path_len, 
                    int32_t seq1_fit, int32_t seq2_fit,
                    int32_t G0, uint8_t *_mem);

/*!
  Creates a cigar array from an alignment path
  @param  path      the Smith-Waterman alignment path
  @param  path_len  the Smith-Waterman alignment path length
  @param  n_cigar   pointer to the returned number of cigar operations
  @return           the cigar array, NULL if the path is NULL or the path length is zero 
  */
uint32_t *
fmap_sw_path2cigar(const fmap_sw_path_t *path, int32_t path_len, int32_t *n_cigar);


/********************
 * global variables *
 ********************/

/*!
  alignment parameters for short read alignment [ACGTN]
  */
extern fmap_sw_param_t fmap_sw_param_short; /* = { 13,  2,  2, fmap_sw_sm_short, 5, 50 }; */
/*!
  alignment parameters for blast read alignment [ACGTN]
  */
extern fmap_sw_param_t fmap_sw_param_blast; /* = {  5,  2,  2, fmap_sw_sm_blast, 5, 50 }; */
/*!
  alignment parameters for ... 
  */
extern fmap_sw_param_t fmap_sw_param_nt2nt; /* = {  8,  2,  2, fmap_sw_sm_nt, 16, 75 }; */
/*!
  alignment parameters for ...
  */
extern fmap_sw_param_t fmap_sw_param_rd2rd; /* = {  1, 19, 19, fmap_sw_sm_read, 16, 75 }; */
/*!
  alignment parameters for ...
  */
extern fmap_sw_param_t fmap_sw_param_aa2aa; /* = { 10,  2,  2, fmap_sw_sm_blosum62, 22, 50 }; */

/*!
  substitution matrix for short read alignment 
  */
extern int32_t fmap_sw_sm_short[];

/*!
  substitution matrix for blast
  */
extern int32_t fmap_sw_sm_blast[];

/*!
  substitution matrix for ...
  */
extern int32_t fmap_sw_sm_nt[];

/*!
  substitution matrix for ...
  */
extern int32_t fmap_sw_sm_read[];

/*!
  substitution matrix for BLOSUM62
  */
extern int32_t fmap_sw_sm_blosum62[];

/*!
  substitution matrix for BLOSUM45
  */
extern int32_t fmap_sw_sm_blosum45[];

/*!
  substitution matrix for human to mouse
  */
extern int32_t           fmap_sw_sm_hs[];

#endif
