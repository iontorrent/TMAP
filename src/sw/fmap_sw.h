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

// TODO: document

#define FMAP_SW_FROM_M 0
#define FMAP_SW_FROM_I 1
#define FMAP_SW_FROM_D 2
#define FMAP_SW_FROM_S 3

#define FMAP_SW_TYPE_LOCAL  0
#define FMAP_SW_TYPE_GLOBAL 1
#define FMAP_SW_TYPE_EXTEND 2

/* This is the smallest integer. It might be CPU-dependent in very RARE cases. */
//#define FMAP_SW_MINOR_INF -1073741823
#define FMAP_SW_MINOR_INF INT32_MIN/2

typedef struct
{
  int32_t gap_open;
  int32_t gap_ext;
  int32_t gap_end;

  int32_t *matrix;
  int32_t row;
  int32_t band_width;
} fmap_sw_param_t;

typedef struct
{
  int32_t i, j;
  unsigned char ctype;
} fmap_sw_path_t;

typedef struct
{
  fmap_sw_path_t *path; /* for advanced users... :-) */
  int32_t path_len; /* for advanced users... :-) */
  int32_t start1, end1; /* start and end of the first sequence, coordinations are 1-based */
  int32_t start2, end2; /* start and end of the second sequence, coordinations are 1-based */
  int32_t score, subo; /* score */

  char *out1, *out2; /* print32_t them, and then you will know */
  char *outm;

  int32_t n_cigar;
  uint32_t *cigar32;
} fmap_sw_aln_t;

fmap_sw_aln_t *
fmap_sw_stdaln_aux(const char *seq1, const char *seq2, const fmap_sw_param_t *ap,
                       int32_t type, int32_t do_align, int32_t len1, int32_t len2);

fmap_sw_aln_t *
fmap_sw_stdaln(const char *seq1, const char *seq2, const fmap_sw_param_t *ap, int32_t type, int32_t do_align);

void fmap_sw_free_fmap_sw_aln_t(fmap_sw_aln_t *aa);

int32_t fmap_sw_global_core(unsigned char *seq1, int32_t len1, unsigned char *seq2, int32_t len2, const fmap_sw_param_t *ap,
                    fmap_sw_path_t *path, int32_t *path_len);
int32_t fmap_sw_local_core(unsigned char *seq1, int32_t len1, unsigned char *seq2, int32_t len2, const fmap_sw_param_t *ap,
                   fmap_sw_path_t *path, int32_t *path_len, int32_t _thres, int32_t *_subo);
int32_t fmap_sw_extend_core(unsigned char *seq1, int32_t len1, unsigned char *seq2, int32_t len2, const fmap_sw_param_t *ap,
                    fmap_sw_path_t *path, int32_t *path_len, int32_t G0, uint8_t *_mem);
uint16_t *fmap_sw_path2cigar(const fmap_sw_path_t *path, int32_t path_len, int32_t *n_cigar);
uint32_t *fmap_sw_path2cigar32(const fmap_sw_path_t *path, int32_t path_len, int32_t *n_cigar);


/********************
 * global variables *
 ********************/

extern fmap_sw_param_t fmap_sw_param_bwa;   /* = { 37,  9,  0, fmap_sw_sm_maq, 5, 50 }; */
extern fmap_sw_param_t fmap_sw_param_blast; /* = {  5,  2,  2, fmap_sw_sm_blast, 5, 50 }; */
extern fmap_sw_param_t fmap_sw_param_nt2nt; /* = { 10,  2,  2, fmap_sw_sm_nt, 16, 75 }; */
extern fmap_sw_param_t fmap_sw_param_aa2aa; /* = { 20, 19, 19, fmap_sw_sm_read, 16, 75 }; */
extern fmap_sw_param_t fmap_sw_param_rd2rd; /* = { 12,  2,  2, fmap_sw_sm_blosum62, 22, 50 }; */

/* common nucleotide score matrix for 16 bases */
extern int32_t           fmap_sw_sm_nt[], fmap_sw_sm_bwa[];

/* BLOSUM62 and BLOSUM45 */
extern int32_t           fmap_sw_sm_blosum62[], fmap_sw_sm_blosum45[];

/* common read for 16 bases. note that read alignment is quite different from common nucleotide alignment */
extern int32_t           fmap_sw_sm_read[];

/* human-mouse score matrix for 4 bases */
extern int32_t           fmap_sw_sm_hs[];

#endif
