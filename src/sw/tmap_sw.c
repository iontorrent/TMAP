/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/* The MIT License

   Copyright (c) 2003-2006, 2008, 2009, by Heng Li <lh3lh3@gmail.com>

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "tmap_sw.h"

/* char -> 17 (=16+1) nucleotides */
uint8_t tmap_sw_nt16_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,16 /*'-'*/,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15, 1,14, 4, 11,15,15, 2, 13,15,15,10, 15, 5,15,15,
    15,15, 3, 6,  8,15, 7, 9,  0,12,15,15, 15,15,15,15,
    15, 1,14, 4, 11,15,15, 2, 13,15,15,10, 15, 5,15,15,
    15,15, 3, 6,  8,15, 7, 9,  0,12,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};
char *tmap_sw_nt16_rev_table = "XAGRCMSVTWKDYHBN-";

/* char -> 5 (=4+1) nucleotides */
uint8_t tmap_sw_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
char *tmap_sw_nt4_rev_table = "AGCTN-";

/* char -> 22 (=20+1+1) amino acids */
uint8_t tmap_sw_aa_table[256] = {
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
    21,21,21,21, 21,21,21,21, 21,21,20,21, 21,22 /*'-'*/,21,21,
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
    21, 0,21, 4,  3, 6,13, 7,  8, 9,21,11, 10,12, 2,21,
    14, 5, 1,15, 16,21,19,17, 21,18,21,21, 21,21,21,21,
    21, 0,21, 4,  3, 6,13, 7,  8, 9,21,11, 10,12, 2,21,
    14, 5, 1,15, 16,21,19,17, 21,18,21,21, 21,21,21,21,
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
    21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21
};
char *tmap_sw_aa_rev_table = "ARNDCQEGHILKMFPSTWYV*X-";
/* 01234567890123456789012 */

/* translation table. They are useless in stdaln.c, but when you realize you need it, you need not write the table again. */
uint8_t tmap_sw_trans_table_eu[66] = {
    11,11, 2, 2,  1, 1,15,15, 16,16,16,16,  9,12, 9, 9,
    6, 6, 3, 3,  7, 7, 7, 7,  0, 0, 0, 0, 19,19,19,19,
    5, 5, 8, 8,  1, 1, 1, 1, 14,14,14,14, 10,10,10,10,
    20,20,18,18, 20,17, 4, 4, 15,15,15,15, 10,10,13,13, 21, 22
};
char *tmap_sw_trans_table_eu_char = "KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFFX";
/* 01234567890123456789012345678901234567890123456789012345678901234 */

int32_t tmap_sw_sm_short[] = {
    11, -19, -19, -19, -13,
    -19, 11, -19, -19, -13,
    -19, -19, 11, -19, -13,
    -19, -19, -19, 11, -13,
    -13, -13, -13, -13, -13
};

int32_t tmap_sw_sm_blast[] = {
    1, -3, -3, -3, -2,
    -3, 1, -3, -3, -2,
    -3, -3, 1, -3, -2,
    -3, -3, -3, 1, -2,
    -2, -2, -2, -2, -2
};

int32_t tmap_sw_sm_nt[] = {
    /*	 X  A  G  R  C  M  S  V  T  W  K  D  Y  H  B  N */
    -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
    -2, 2,-1, 1,-2, 1,-2, 0,-2, 1,-2, 0,-2, 0,-2, 0,
    -2,-1, 2, 1,-2,-2, 1, 0,-2,-2, 1, 0,-2,-2, 0, 0,
    -2, 1, 1, 1,-2,-1,-1, 0,-2,-1,-1, 0,-2, 0, 0, 0,
    -2,-2,-2,-2, 2, 1, 1, 0,-1,-2,-2,-2, 1, 0, 0, 0,
    -2, 1,-2,-1, 1, 1,-1, 0,-2,-1,-2, 0,-1, 0, 0, 0,
    -2,-2, 1,-1, 1,-1, 1, 0,-2,-2,-1, 0,-1, 0, 0, 0,
    -2, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0,
    -2,-2,-2,-2,-1,-2,-2,-2, 2, 1, 1, 0, 1, 0, 0, 0,
    -2, 1,-2,-1,-2,-1,-2, 0, 1, 1,-1, 0,-1, 0, 0, 0,
    -2,-2, 1,-1,-2,-2,-1, 0, 1,-1, 1, 0,-1, 0, 0, 0,
    -2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -2,-2,-2,-2, 1,-1,-1, 0, 1,-1,-1, 0, 1, 0, 0, 0,
    -2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

int32_t tmap_sw_sm_read[] = {
    /*	  X   A   G   R   C   M   S   V   T   W   K   D   Y   H   B   N  */
    -17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,
    -17,  2,-17,  1,-17,  1,-17,  0,-17,  1,-17,  0,-17,  0,-17,  0,
    -17,-17,  2,  1,-17,-17,  1,  0,-17,-17,  1,  0,-17,-17,  0,  0,
    -17,  1,  1,  1,-17,-17,-17,  0,-17,-17,-17,  0,-17,  0,  0,  0,
    -17,-17,-17,-17,  2,  1,  1,  0,-17,-17,-17,-17,  1,  0,  0,  0,
    -17,  1,-17,-17,  1,  1,-17,  0,-17,-17,-17,  0,-17,  0,  0,  0,
    -17,-17,  1,-17,  1,-17,  1,  0,-17,-17,-17,  0,-17,  0,  0,  0,
    -17,  0,  0,  0,  0,  0,  0,  0,-17,  0,  0,  0,  0,  0,  0,  0,
    -17,-17,-17,-17,-17,-17,-17,-17,  2,  1,  1,  0,  1,  0,  0,  0,
    -17,  1,-17,-17,-17,-17,-17,  0,  1,  1,-17,  0,-17,  0,  0,  0,
    -17,-17,  1,-17,-17,-17,-17,  0,  1,-17,  1,  0,-17,  0,  0,  0,
    -17,  0,  0,  0,-17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    -17,-17,-17,-17,  1,-17,-17,  0,  1,-17,-17,  0,  1,  0,  0,  0,
    -17,  0,-17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    -17,-17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    -17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

int32_t tmap_sw_sm_blosum62[] = {
    /*	 A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  *  X */
    4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-4, 0,
    -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-4,-1,
    -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,-4,-1,
    -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,-4,-1,
    0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-4,-2,
    -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,-4,-1,
    -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,-4,-1,
    0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-4,-1,
    -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,-4,-1,
    -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-4,-1,
    -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-1,
    -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,-4,-1,
    -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-4,-1,
    -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-4,-1,
    -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-4,-2,
    1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,-4, 0,
    0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-4, 0,
    -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,
    -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-4,-1,
    0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-4,-1,
    -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,-4,
    0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-4,-1
};

int32_t tmap_sw_sm_blosum45[] = {
    /*	 A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  *  X */
    5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-2,-2, 0,-5, 0,
    -2, 7, 0,-1,-3, 1, 0,-2, 0,-3,-2, 3,-1,-2,-2,-1,-1,-2,-1,-2,-5,-1,
    -1, 0, 6, 2,-2, 0, 0, 0, 1,-2,-3, 0,-2,-2,-2, 1, 0,-4,-2,-3,-5,-1,
    -2,-1, 2, 7,-3, 0, 2,-1, 0,-4,-3, 0,-3,-4,-1, 0,-1,-4,-2,-3,-5,-1,
    -1,-3,-2,-3,12,-3,-3,-3,-3,-3,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1,-5,-2,
    -1, 1, 0, 0,-3, 6, 2,-2, 1,-2,-2, 1, 0,-4,-1, 0,-1,-2,-1,-3,-5,-1,
    -1, 0, 0, 2,-3, 2, 6,-2, 0,-3,-2, 1,-2,-3, 0, 0,-1,-3,-2,-3,-5,-1,
    0,-2, 0,-1,-3,-2,-2, 7,-2,-4,-3,-2,-2,-3,-2, 0,-2,-2,-3,-3,-5,-1,
    -2, 0, 1, 0,-3, 1, 0,-2,10,-3,-2,-1, 0,-2,-2,-1,-2,-3, 2,-3,-5,-1,
    -1,-3,-2,-4,-3,-2,-3,-4,-3, 5, 2,-3, 2, 0,-2,-2,-1,-2, 0, 3,-5,-1,
    -1,-2,-3,-3,-2,-2,-2,-3,-2, 2, 5,-3, 2, 1,-3,-3,-1,-2, 0, 1,-5,-1,
    -1, 3, 0, 0,-3, 1, 1,-2,-1,-3,-3, 5,-1,-3,-1,-1,-1,-2,-1,-2,-5,-1,
    -1,-1,-2,-3,-2, 0,-2,-2, 0, 2, 2,-1, 6, 0,-2,-2,-1,-2, 0, 1,-5,-1,
    -2,-2,-2,-4,-2,-4,-3,-3,-2, 0, 1,-3, 0, 8,-3,-2,-1, 1, 3, 0,-5,-1,
    -1,-2,-2,-1,-4,-1, 0,-2,-2,-2,-3,-1,-2,-3, 9,-1,-1,-3,-3,-3,-5,-1,
    1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-3,-1,-2,-2,-1, 4, 2,-4,-2,-1,-5, 0,
    0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1, 2, 5,-3,-1, 0,-5, 0,
    -2,-2,-4,-4,-5,-2,-3,-2,-3,-2,-2,-2,-2, 1,-3,-4,-3,15, 3,-3,-5,-2,
    -2,-1,-2,-2,-3,-1,-2,-3, 2, 0, 0,-1, 0, 3,-3,-2,-1, 3, 8,-1,-5,-1,
    0,-2,-3,-3,-1,-3,-3,-3,-3, 3, 1,-2, 1, 0,-3,-1, 0,-3,-1, 5,-5,-1,
    -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5, 1,-5,
    0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-2,-1,-1,-5,-1
};

int32_t tmap_sw_sm_hs[] = {
    /*     A    G    C    T    N */
    91, -31,-114,-123, -44,
    -31, 100,-125,-114, -42,
    -123,-125, 100, -31, -42,
    -114,-114, -31,  91, -42,
    -44, -42, -42, -42, -43
};

/********************/
/* START OF align.c */
/********************/

tmap_sw_param_t tmap_sw_param_short   = { 13,  2,  2, tmap_sw_sm_short, 5, 50 }; /* e=1.21%; T=2.18 */
tmap_sw_param_t tmap_sw_param_blast   = {  5,  2,  2, tmap_sw_sm_blast, 5, 50 };
tmap_sw_param_t tmap_sw_param_nt2nt   = {  8,  2,  2, tmap_sw_sm_nt, 16, 75 };
tmap_sw_param_t tmap_sw_param_rd2rd   = {  1, 19, 19, tmap_sw_sm_read, 16, 75 };
tmap_sw_param_t tmap_sw_param_aa2aa   = { 10,  2,  2, tmap_sw_sm_blosum62, 22, 50 };

tmap_sw_aln_t *
tmap_sw_aln_init()
{
  tmap_sw_aln_t *aa;
  aa = tmap_malloc(sizeof(tmap_sw_aln_t), "aa");
  aa->path = 0;
  aa->out1 = aa->out2 = aa->outm = 0;
  aa->path_len = 0;
  return aa;
}

void 
tmap_sw_aln_destroy(tmap_sw_aln_t *aa)
{
  free(aa->path); free(aa->cigar32);
  free(aa->out1); free(aa->out2); free(aa->outm);
  free(aa);
}

/***************************/
/* START OF common_align.c */
/***************************/

#define TMAP_SW_LOCAL_OVERFLOW_THRESHOLD 32000
#define TMAP_SW_LOCAL_OVERFLOW_REDUCE 16000
#define TMAP_SW_NT_LOCAL_SCORE int32_t
#define TMAP_SW_NT_LOCAL_SHIFT 16
#define TMAP_SW_NT_LOCAL_MASK 0xffff

#define tmap_sw_set_match(MM, cur, p, sc, right_j) \
{ \
  if((p)->match_score > (p)->ins_score \
     || (0 == right_j && (p)->match_score == (p)->ins_score)) { \
      if((p)->match_score > (p)->del_score \
          || (0 == right_j && (p)->match_score == (p)->del_score)) { \
          (MM) = (p)->match_score + (sc); (cur)->match_from = TMAP_SW_FROM_M; \
      } else { \
          (MM) = (p)->del_score + (sc); (cur)->match_from = TMAP_SW_FROM_D; \
      } \
  } else { \
      if((p)->ins_score > (p)->del_score) { \
          (MM) = (p)->ins_score + (sc); (cur)->match_from = TMAP_SW_FROM_I; \
      } else { \
          (MM) = (p)->del_score + (sc); (cur)->match_from = TMAP_SW_FROM_D; \
      } \
  } \
}

#define tmap_sw_set_ins(II, cur, p, right_j) \
{ \
  if((p)->match_score - gap_open > (p)->ins_score \
     || (1 == right_j && (p)->match_score - gap_open == (p)->ins_score)) { \
      (cur)->ins_from = TMAP_SW_FROM_M; \
      (II) = (p)->match_score - gap_open - gap_ext; \
  } else { \
      (cur)->ins_from = TMAP_SW_FROM_I; \
      (II) = (p)->ins_score - gap_ext; \
  } \
}

#define tmap_sw_set_end_ins(II, cur, p, right_j) \
{ \
  if(gap_end >= 0) { \
      if((p)->match_score - gap_open > (p)->ins_score \
         || (1 == right_j && (p)->match_score - gap_open == (p)->ins_score)) { \
          (cur)->ins_from = TMAP_SW_FROM_M; \
          (II) = (p)->match_score - gap_open - gap_end; \
      } else { \
          (cur)->ins_from = TMAP_SW_FROM_I; \
          (II) = (p)->ins_score - gap_end; \
      } \
  } else tmap_sw_set_ins(II, cur, p, right_j); \
}

#define tmap_sw_set_del(DD, cur, p, right_j) \
{ \
  if((p)->match_score - gap_open > (p)->del_score \
     || (1 == right_j && (p)->match_score - gap_open == (p)->del_score)) { \
      (cur)->del_from = TMAP_SW_FROM_M; \
      (DD) = (p)->match_score - gap_open - gap_ext; \
  } else { \
      (cur)->del_from = TMAP_SW_FROM_D; \
      (DD) = (p)->del_score - gap_ext; \
  } \
}

#define tmap_sw_set_end_del(DD, cur, p, right_j) \
{ \
  if(gap_end >= 0) { \
      if((p)->match_score - gap_open > (p)->del_score \
         || (1 == right_j && (p)->match_score - gap_open == (p)->del_score)) { \
          (cur)->del_from = TMAP_SW_FROM_M; \
          (DD) = (p)->match_score - gap_open - gap_end; \
      } else { \
          (cur)->del_from = TMAP_SW_FROM_D; \
          (DD) = (p)->del_score - gap_end; \
      } \
  } else tmap_sw_set_del(DD, cur, p, right_j); \
}

/* build score profile for accelerating alignment, in theory */
void 
tmap_sw_score_array_init(uint8_t *seq, int32_t len, int32_t row, int32_t *score_matrix, int32_t **s_array)
{
  int32_t *tmp, *tmp2, i, k;
  for(i = 0; i != row; ++i) {
      tmp = score_matrix + i * row;
      tmp2 = s_array[i];
      for(k = 0; k != len; ++k)
        tmp2[k] = tmp[seq[k]];
  }
}

/***************************
 * banded global alignment *
 ***************************/
int32_t 
tmap_sw_global_core(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
                    tmap_sw_path_t *path, int32_t *path_len, int32_t right_j)
{
  register int32_t i, j;
  tmap_sw_dpcell_t **dpcell, *q;
  tmap_sw_dpscore_t *curr, *last, *s;
  tmap_sw_path_t *p;
  int32_t b1, b2, tmp_end;
  int32_t *mat, end, max;
  uint8_t type, ctype;

  int32_t gap_open, gap_ext, gap_end, b;
  int32_t *score_matrix, N_MATRIX_ROW;

  /* initialize some align-related parameters. just for compatibility */
  gap_open = ap->gap_open;
  gap_ext = ap->gap_ext;
  gap_end = ap->gap_end;
  b = ap->band_width;
  score_matrix = ap->matrix;
  N_MATRIX_ROW = ap->row;

  if(len1 == 0 || len2 == 0) {
      *path_len = 0;
      return TMAP_SW_MINOR_INF;
  }
  /* calculate b1 and b2 */
  if(len1 > len2) {
      b1 = len1 - len2 + b;
      b2 = b;
  } else {
      b1 = b;
      b2 = len2 - len1 + b;
  }
  if(b1 > len1) b1 = len1;
  if(b2 > len2) b2 = len2;
  --seq1; --seq2;

  /* allocate memory */
  end = (b1 + b2 <= len1)? (b1 + b2 + 1) : (len1 + 1);
  dpcell = tmap_malloc(sizeof(tmap_sw_dpcell_t*) * (len2 + 1), "dpcell");
  for(j = 0; j <= len2; ++j)
    dpcell[j] = tmap_malloc(sizeof(tmap_sw_dpcell_t) * end, "dpcell[j]");
  for(j = b2 + 1; j <= len2; ++j)
    dpcell[j] -= j - b2;
  curr = tmap_malloc(sizeof(tmap_sw_dpscore_t) * (len1 + 1), "curr");
  last = tmap_malloc(sizeof(tmap_sw_dpscore_t) * (len1 + 1), "last");

  /* set first row */
  TMAP_SW_SET_INF(*curr); curr->match_score = 0;
  for(i = 1, s = curr + 1; i < b1; ++i, ++s) {
      TMAP_SW_SET_INF(*s);
      tmap_sw_set_end_del(s->del_score, dpcell[0] + i, s - 1, right_j);
  }
  s = curr; curr = last; last = s;

  /* core dynamic programming, part 1 */
  tmp_end = (b2 < len2)? b2 : len2 - 1;
  for(j = 1; j <= tmp_end; ++j) {
      q = dpcell[j]; s = curr; TMAP_SW_SET_INF(*s);
      tmap_sw_set_end_ins(s->ins_score, q, last, right_j);
      end = (j + b1 <= len1 + 1)? (j + b1 - 1) : len1;
      mat = score_matrix + seq2[j] * N_MATRIX_ROW;
      ++s; ++q;
      for(i = 1; i != end; ++i, ++s, ++q) {
          tmap_sw_set_match(s->match_score, q, last + i - 1, mat[seq1[i]], right_j); /* this will change s->match_score ! */
          tmap_sw_set_ins(s->ins_score, q, last + i, right_j);
          tmap_sw_set_del(s->del_score, q, s - 1, right_j);
      }
      tmap_sw_set_match(s->match_score, q, last + i - 1, mat[seq1[i]], right_j);
      tmap_sw_set_del(s->del_score, q, s - 1, right_j);
      if(j + b1 - 1 > len1) { /* bug fixed, 040227 */
          tmap_sw_set_end_ins(s->ins_score, q, last + i, right_j);
      } else s->ins_score = TMAP_SW_MINOR_INF;
      s = curr; curr = last; last = s;
  }
  /* last row for part 1, use tmap_sw_set_end_del() instead of tmap_sw_set_del() */
  if(j == len2 && b2 != len2 - 1) {
      q = dpcell[j]; s = curr; TMAP_SW_SET_INF(*s);
      tmap_sw_set_end_ins(s->ins_score, q, last, right_j);
      end = (j + b1 <= len1 + 1)? (j + b1 - 1) : len1;
      mat = score_matrix + seq2[j] * N_MATRIX_ROW;
      ++s; ++q;
      for(i = 1; i != end; ++i, ++s, ++q) {
          tmap_sw_set_match(s->match_score, q, last + i - 1, mat[seq1[i]], right_j); /* this will change s->match_score ! */
          tmap_sw_set_ins(s->ins_score, q, last + i, right_j);
          tmap_sw_set_end_del(s->del_score, q, s - 1, right_j);
      }
      tmap_sw_set_match(s->match_score, q, last + i - 1, mat[seq1[i]], right_j);
      tmap_sw_set_end_del(s->del_score, q, s - 1, right_j);
      if(j + b1 - 1 > len1) { /* bug fixed, 040227 */
          tmap_sw_set_end_ins(s->ins_score, q, last + i, right_j);
      } else s->ins_score = TMAP_SW_MINOR_INF;
      s = curr; curr = last; last = s;
      ++j;
  }

  /* core dynamic programming, part 2 */
  for(; j <= len2 - b2 + 1; ++j) {
      TMAP_SW_SET_INF(curr[j - b2]);
      mat = score_matrix + seq2[j] * N_MATRIX_ROW;
      end = j + b1 - 1;
      for(i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i != end; ++i, ++s, ++q) {
          tmap_sw_set_match(s->match_score, q, last + i - 1, mat[seq1[i]], right_j);
          tmap_sw_set_ins(s->ins_score, q, last + i, right_j);
          tmap_sw_set_del(s->del_score, q, s - 1, right_j);
      }
      tmap_sw_set_match(s->match_score, q, last + i - 1, mat[seq1[i]], right_j);
      tmap_sw_set_del(s->del_score, q, s - 1, right_j);
      s->ins_score = TMAP_SW_MINOR_INF;
      s = curr; curr = last; last = s;
  }

  /* core dynamic programming, part 3 */
  for(; j < len2; ++j) {
      TMAP_SW_SET_INF(curr[j - b2]);
      mat = score_matrix + seq2[j] * N_MATRIX_ROW;
      for(i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i < len1; ++i, ++s, ++q) {
          tmap_sw_set_match(s->match_score, q, last + i - 1, mat[seq1[i]], right_j);
          tmap_sw_set_ins(s->ins_score, q, last + i, right_j);
          tmap_sw_set_del(s->del_score, q, s - 1, 0 );
      }
      tmap_sw_set_match(s->match_score, q, last + len1 - 1, mat[seq1[i]], right_j);
      tmap_sw_set_end_ins(s->ins_score, q, last + i, right_j);
      tmap_sw_set_del(s->del_score, q, s - 1, right_j);
      s = curr; curr = last; last = s;
  }
  /* last row */
  if(j == len2) {
      TMAP_SW_SET_INF(curr[j - b2]);
      mat = score_matrix + seq2[j] * N_MATRIX_ROW;
      for(i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i < len1; ++i, ++s, ++q) {
          tmap_sw_set_match(s->match_score, q, last + i - 1, mat[seq1[i]], right_j);
          tmap_sw_set_ins(s->ins_score, q, last + i, right_j);
          tmap_sw_set_end_del(s->del_score, q, s - 1, right_j);
      }
      tmap_sw_set_match(s->match_score, q, last + len1 - 1, mat[seq1[i]], right_j);
      tmap_sw_set_end_ins(s->ins_score, q, last + i, right_j);
      tmap_sw_set_end_del(s->del_score, q, s - 1, right_j);
      s = curr; curr = last; last = s;
  }

  /* backtrace */
  i = len1; j = len2;
  q = dpcell[j] + i;
  s = last + len1;
  max = s->match_score; type = q->match_from; ctype = TMAP_SW_FROM_M;
  if(s->ins_score > max) { max = s->ins_score; type = q->ins_from; ctype = TMAP_SW_FROM_I; }
  if(s->del_score > max) { max = s->del_score; type = q->del_from; ctype = TMAP_SW_FROM_D; }

  p = path;
  p->ctype = ctype; p->i = i; p->j = j; 
  ++p;
  do {
      switch(ctype) {
        case TMAP_SW_FROM_M: --i; --j; break;
        case TMAP_SW_FROM_I: --j; break;
        case TMAP_SW_FROM_D: --i; break;
      }
      q = dpcell[j] + i;
      ctype = type;
      switch(type) {
        case TMAP_SW_FROM_M: type = q->match_from; break;
        case TMAP_SW_FROM_I: type = q->ins_from; break;
        case TMAP_SW_FROM_D: type = q->del_from; break;
      }
      p->ctype = ctype; p->i = i; p->j = j;
      ++p;
  } while (i || j);
  (*path_len) = p - path - 1;

  /* free memory */
  for(j = b2 + 1; j <= len2; ++j)
    dpcell[j] += j - b2;
  for(j = 0; j <= len2; ++j)
    free(dpcell[j]);
  free(dpcell);
  free(curr); free(last);

  return max;
}

/*************************************************
 * local alignment combined with banded strategy *
 *************************************************/
int32_t 
tmap_sw_local_core(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
                   tmap_sw_path_t *path, int32_t *path_len, int32_t _thres, int32_t *_subo)
{
  register TMAP_SW_NT_LOCAL_SCORE *s;
  register int32_t i;
  int32_t q, r, qr, qr_shift;
  int32_t **s_array, *score_array;
  int32_t e, f;
  int32_t is_overflow, of_base;
  TMAP_SW_NT_LOCAL_SCORE *eh, curr_h, last_h, curr_last_h;
  int32_t j, start_i, start_j, end_i, end_j;
  tmap_sw_path_t *p;
  int32_t score_f, score_r, score_g;
  int32_t start, end, max_score;
  int32_t thres, *suba, *ss;

  int32_t gap_open, gap_ext, b;
  int32_t *score_matrix, N_MATRIX_ROW;

  /* initialize some align-related parameters. just for compatibility */
  gap_open = ap->gap_open;
  gap_ext = ap->gap_ext;
  b = ap->band_width;
  score_matrix = ap->matrix;
  N_MATRIX_ROW = ap->row;
  thres = _thres > 0? _thres : -_thres;

  if(len1 == 0 || len2 == 0) return TMAP_SW_MINOR_INF;

  /* allocate memory */
  suba = tmap_malloc(sizeof(int32_t) * (len2 + 1), "suba");
  eh = tmap_malloc(sizeof(TMAP_SW_NT_LOCAL_SCORE) * (len1 + 1), "eh");
  s_array = tmap_malloc(sizeof(int32_t*) * N_MATRIX_ROW, "s_array");
  for(i = 0; i != N_MATRIX_ROW; ++i)
    s_array[i] = tmap_malloc(sizeof(int32_t) * len1, "s_array[i]");
  /* initialization */
  tmap_sw_score_array_init(seq1, len1, N_MATRIX_ROW, score_matrix, s_array);
  q = gap_open;
  r = gap_ext;
  qr = q + r;
  qr_shift = (qr+1) << TMAP_SW_NT_LOCAL_SHIFT;
  start_i = start_j = end_i = end_j = 0;
  for(i = 0, max_score = 0; i != N_MATRIX_ROW * N_MATRIX_ROW; ++i)
    if(max_score < score_matrix[i]) max_score = score_matrix[i];
  /* convert the coordinate */
  --seq1; --seq2;
  for(i = 0; i != N_MATRIX_ROW; ++i) --s_array[i];

  /* forward dynamic programming */
  for(i = 0, s = eh; i != len1 + 1; ++i, ++s) *s = 0;
  score_f = 0;
  is_overflow = of_base = 0;
  suba[0] = 0;
  for(j = 1, ss = suba + 1; j <= len2; ++j, ++ss) {
      int32_t subo = 0;
      last_h = f = 0;
      score_array = s_array[seq2[j]];
      if(is_overflow) { /* adjust eh[] array if overflow occurs. */
          /* If TMAP_SW_LOCAL_OVERFLOW_REDUCE is too small, optimal alignment might be missed.
           * If it is too large, this block will be excuted frequently and therefore
           * slow down the whole program.
           * Acually, smaller TMAP_SW_LOCAL_OVERFLOW_REDUCE might also help to reduce the
           * number of assignments because it sets some cells to zero when overflow
           * happens. */
          int32_t tmp, tmp2;
          score_f -= TMAP_SW_LOCAL_OVERFLOW_REDUCE;
          of_base += TMAP_SW_LOCAL_OVERFLOW_REDUCE;
          is_overflow = 0;
          for(i = 1, s = eh; i <= len1 + 1; ++i, ++s) {
              tmp = *s >> TMAP_SW_NT_LOCAL_SHIFT; tmp2 = *s & TMAP_SW_NT_LOCAL_MASK;
              if(tmp2 < TMAP_SW_LOCAL_OVERFLOW_REDUCE) tmp2 = 0;
              else tmp2 -= TMAP_SW_LOCAL_OVERFLOW_REDUCE;
              if(tmp < TMAP_SW_LOCAL_OVERFLOW_REDUCE) tmp = 0;
              else tmp -= TMAP_SW_LOCAL_OVERFLOW_REDUCE;
              *s = (tmp << TMAP_SW_NT_LOCAL_SHIFT) | tmp2;
          }
      }
      for(i = 1, s = eh; i != len1 + 1; ++i, ++s) {
          /* prepare for calculate current h */
          curr_h = (*s >> TMAP_SW_NT_LOCAL_SHIFT) + score_array[i];
          if(curr_h < 0) curr_h = 0;
          if(last_h > 0) { /* initialize f */
              f = (f > last_h - q)? f - r : last_h - qr;
              if(curr_h < f) curr_h = f;
          }
          if(TMAP_SW_UNLIKELY(*(s+1) >= qr_shift)) { /* initialize e */
              curr_last_h = *(s+1) >> TMAP_SW_NT_LOCAL_SHIFT;
              e = ((*s & TMAP_SW_NT_LOCAL_MASK) > curr_last_h - q)? (*s & TMAP_SW_NT_LOCAL_MASK) - r : curr_last_h - qr;
              if(curr_h < e) curr_h = e;
              *s = (last_h << TMAP_SW_NT_LOCAL_SHIFT) | e;
          } else *s = last_h << TMAP_SW_NT_LOCAL_SHIFT; /* e = 0 */
          last_h = curr_h;
          if(TMAP_SW_UNLIKELY(subo < curr_h)) subo = curr_h;
          if(TMAP_SW_UNLIKELY(score_f < curr_h)) {
              score_f = curr_h; end_i = i; end_j = j;
              if(score_f > TMAP_SW_LOCAL_OVERFLOW_THRESHOLD) is_overflow = 1;
          }
      }
      *s = last_h << TMAP_SW_NT_LOCAL_SHIFT;
      *ss = subo + of_base;
  }
  score_f += of_base;

  if(score_f < thres) { /* no matching residue at all, 090218 */
      *path_len = 0;
      goto end_func;
  }
  if(path == 0) goto end_func; /* skip path-filling */

  /* reverse dynamic programming */
  for(i = end_i, s = eh + end_i; i >= 0; --i, --s) *s = 0;
  if(end_i == 0 || end_j == 0) goto end_func; /* no local match */
  score_r = score_matrix[seq1[end_i] * N_MATRIX_ROW + seq2[end_j]];
  is_overflow = of_base = 0;
  start_i = end_i; start_j = end_j;
  eh[end_i] = ((TMAP_SW_NT_LOCAL_SCORE)(qr + score_r)) << TMAP_SW_NT_LOCAL_SHIFT; /* in order to initialize f and e, 040408 */
  start = end_i - 1;
  end = end_i - 3;
  if(end <= 0) end = 0;

  /* second pass DP can be done in a band, speed will thus be enhanced */
  for(j = end_j - 1; j != 0; --j) {
      last_h = f = 0;
      score_array = s_array[seq2[j]];
      if(is_overflow) { /* adjust eh[] array if overflow occurs. */
          int32_t tmp, tmp2;
          score_r -= TMAP_SW_LOCAL_OVERFLOW_REDUCE;
          of_base += TMAP_SW_LOCAL_OVERFLOW_REDUCE;
          is_overflow = 0;
          for(i = start, s = eh + start + 1; i >= end; --i, --s) {
              tmp = *s >> TMAP_SW_NT_LOCAL_SHIFT; tmp2 = *s & TMAP_SW_NT_LOCAL_MASK;
              if(tmp2 < TMAP_SW_LOCAL_OVERFLOW_REDUCE) tmp2 = 0;
              else tmp2 -= TMAP_SW_LOCAL_OVERFLOW_REDUCE;
              if(tmp < TMAP_SW_LOCAL_OVERFLOW_REDUCE) tmp = 0;
              else tmp -= TMAP_SW_LOCAL_OVERFLOW_REDUCE;
              *s = (tmp << TMAP_SW_NT_LOCAL_SHIFT) | tmp2;
          }
      }
      for(i = start, s = eh + start + 1; i != end; --i, --s) {
          /* prepare for calculate current h */
          curr_h = (*s >> TMAP_SW_NT_LOCAL_SHIFT) + score_array[i];
          if(curr_h < 0) curr_h = 0;
          if(last_h > 0) { /* initialize f */
              f = (f > last_h - q)? f - r : last_h - qr;
              if(curr_h < f) curr_h = f;
          }
          curr_last_h = *(s-1) >> TMAP_SW_NT_LOCAL_SHIFT;
          e = ((*s & TMAP_SW_NT_LOCAL_MASK) > curr_last_h - q)? (*s & TMAP_SW_NT_LOCAL_MASK) - r : curr_last_h - qr;
          if(e < 0) e = 0;
          if(curr_h < e) curr_h = e;
          *s = (last_h << TMAP_SW_NT_LOCAL_SHIFT) | e;
          last_h = curr_h;
          if(score_r < curr_h) {
              score_r = curr_h; start_i = i; start_j = j;
              if(score_r + of_base - qr == score_f) {
                  j = 1; break;
              }
              if(score_r > TMAP_SW_LOCAL_OVERFLOW_THRESHOLD) is_overflow = 1;
          }
      }
      *s = last_h << TMAP_SW_NT_LOCAL_SHIFT;
      /* recalculate start and end, the boundaries of the band */
      if((eh[start] >> TMAP_SW_NT_LOCAL_SHIFT) <= qr) --start;
      if(start <= 0) start = 0;
      end = start_i - (start_j - j) - (score_r + of_base + (start_j - j) * max_score) / r - 1;
      if(end <= 0) end = 0;
  }

  if(_subo) {
      int32_t tmp2 = 0, tmp = (int32_t)(start_j - .33 * (end_j - start_j) + .499);
      for(j = 1; j <= tmp; ++j)
        if(tmp2 < suba[j]) tmp2 = suba[j];
      tmp = (int32_t)(end_j + .33 * (end_j - start_j) + .499);
      for(j = tmp; j <= len2; ++j)
        if(tmp2 < suba[j]) tmp2 = suba[j];
      *_subo = tmp2;
  }

  if(path_len == 0) {
      path[0].i = start_i; path[0].j = start_j;
      path[1].i = end_i; path[1].j = end_j;
      goto end_func;
  }

  score_r += of_base;
  score_r -= qr;

#ifdef DEBUG
  /* this seems not to be a bug */
  if(score_f != score_r)
    fprintf(stderr, "[tmap_sw_local_core] unknown flaw occurs: score_f(%d) != score_r(%d)\n", score_f, score_r);
#endif

  if(_thres > 0) { /* call global alignment to fill the path */
      score_g = 0;
      j = (end_i - start_i > end_j - start_j)? end_i - start_i : end_j - start_j;
      ++j; /* j is the maximum band_width */
      for(i = ap->band_width;; i <<= 1) {
          tmap_sw_param_t ap_real = *ap;
          ap_real.gap_end = -1;
          ap_real.band_width = i;
          score_g = tmap_sw_global_core(seq1 + start_i, end_i - start_i + 1, seq2 + start_j,
                                        end_j - start_j + 1, &ap_real, path, path_len, 0);
          if(score_g == score_r || score_f == score_g) break;
          if(i > j) break;
      }
#ifdef DEBUG
      // this seems not to be a bug
      if(score_r > score_g && score_f > score_g) {
          fprintf(stderr, "[tmap_sw_local_core] Potential bug: (%d,%d) > %d\n", score_f, score_r, score_g);
          score_f = score_r = -1;
      } else score_f = score_g;
#endif
      score_f = score_g;

      /* convert coordinate */
      for(p = path + *path_len - 1; p >= path; --p) {
          p->i += start_i - 1;
          p->j += start_j - 1;
      }
  } else { /* just store the start and end */
      *path_len = 2;
      path[1].i = start_i; path[1].j = start_j;
      path->i = end_i; path->j = end_j;
  }

end_func:
  /* free */
  free(eh); free(suba);
  for(i = 0; i != N_MATRIX_ROW; ++i) {
      ++s_array[i];
      free(s_array[i]);
  }
  free(s_array);
  return score_f;
}

/*
static int32_t 
tmap_sw_extend_aux(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
                    tmap_sw_path_t *path, int32_t *path_len, int32_t prev_score, uint8_t *_mem)
{
  int32_t q, r, qr;
  int32_t **s_array, *score_array;
  int32_t is_overflow, of_base;
  uint32_t *eh;
  int32_t i, j, end_i, end_j;
  int32_t score, start, end;
  int32_t *score_matrix, N_MATRIX_ROW;
  uint8_t *mem, *_p;

  // initialize some align-related parameters. just for compatibility 
  q = ap->gap_open;
  r = ap->gap_ext;
  qr = q + r;
  score_matrix = ap->matrix;
  N_MATRIX_ROW = ap->row;

  if(len1 == 0 || len2 == 0) return TMAP_SW_MINOR_INF;

  // allocate memory 
  mem = _mem ? _mem : tmap_calloc((len1 + 2) * (N_MATRIX_ROW + 1), sizeof(uint32_t), "mem");
  _p = mem;
  eh = (uint32_t*)_p, _p += 4 * (len1 + 2);
  s_array = tmap_calloc(N_MATRIX_ROW, sizeof(void*), "s_array");
  for(i = 0; i != N_MATRIX_ROW; ++i)
    s_array[i] = (int32_t*)_p, _p += 4 * len1;
  // initialization 
  tmap_sw_score_array_init(seq1, len1, N_MATRIX_ROW, score_matrix, s_array);
  start = 1; end = 2;
  end_i = end_j = 0;
  score = 0;
  is_overflow = of_base = 0;
  // convert the coordinate 
  --seq1; --seq2;
  for(i = 0; i != N_MATRIX_ROW; ++i) --s_array[i];
  // dynamic programming 
  memset(eh, 0, 4 * (len1 + 2));
  eh[1] = (uint32_t)prev_score<<16;
  for(j = 1; j <= len2; ++j) {
      int32_t _start, _end;
      int32_t h1 = 0, f = 0;
      score_array = s_array[seq2[j]];
      // set start and end 
      _start = j - ap->band_width;
      if(_start < 1) _start = 1;
      if(_start > start) start = _start;
      _end = j + ap->band_width;
      if(_end > len1 + 1) _end = len1 + 1;
      if(_end < end) end = _end;
      if(start == end) break;
      // adjust eh[] array if overflow occurs. 
      if(is_overflow) {
          int32_t tmp, tmp2;
          score -= TMAP_SW_LOCAL_OVERFLOW_REDUCE;
          of_base += TMAP_SW_LOCAL_OVERFLOW_REDUCE;
          is_overflow = 0;
          for(i = start; i <= end; ++i) {
              uint32_t *s = &eh[i];
              tmp = *s >> 16; tmp2 = *s & 0xffff;
              if(tmp2 < TMAP_SW_LOCAL_OVERFLOW_REDUCE) tmp2 = 0;
              else tmp2 -= TMAP_SW_LOCAL_OVERFLOW_REDUCE;
              if(tmp < TMAP_SW_LOCAL_OVERFLOW_REDUCE) tmp = 0;
              else tmp -= TMAP_SW_LOCAL_OVERFLOW_REDUCE;
              *s = (tmp << 16) | tmp2;
          }
      }
      _start = _end = 0;
      // the inner loop 
      for(i = start; i < end; ++i) {
          // At the beginning of each cycle:
          // eh[i] -> h[j-1,i-1]<<16 | e[j,i]
          // f     -> f[j,i]
          // h1    -> h[j,i-1]
             
          uint32_t *s = &eh[i];
          int32_t h = (int32_t)(*s >> 16);
          int32_t e = *s & 0xffff; // this is e[j,i] 
          *s = (uint32_t)h1 << 16; // eh[i] now stores h[j,i-1]<<16 
          h += h? score_array[i] : 0; // this is left_core() specific 
          // calculate h[j,i]; don't need to test 0, as {e,f}>=0 
          h = h > e? h : e;
          h = h > f? h : f; // h now is h[j,i] 
          h1 = h;
          if(h > 0) {
              if(_start == 0) _start = i;
              _end = i;
              if(score < h) {
                  score = h; end_i = i; end_j = j;
                  if(score > TMAP_SW_LOCAL_OVERFLOW_THRESHOLD) is_overflow = 1;
              }
          }
          // calculate e[j+1,i] and f[j,i+1] 
          h -= qr;
          h = h > 0? h : 0;
          e -= r;
          e = e > h? e : h;
          f -= r;
          f = f > h? f : h;
          *s |= e;
      }			
      eh[end] = h1 << 16;
      // recalculate start and end, the boundaries of the band 
      if(_end <= 0) break; // no cell in this row has a positive score 
      start = _start;
      end = _end + 3;
  }

  score += of_base - 1;

  if(NULL == path) {
      // do nothing
  }
  else if(NULL == path_len) {
      path[0].i = end_i; path[0].j = end_j;
  }
  else if(score <= 0) {
      (*path_len) = 0;
  }
  else {
      // call global alignment to fill the path 
      int32_t score_g = 0;
      j = (end_i - 1 > end_j - 1)? end_i - 1 : end_j - 1;
      ++j; // j is the maximum band_width 
      for(i = ap->band_width;; i <<= 1) {
          tmap_sw_param_t ap_real = *ap;
          ap_real.gap_end = -1;
          ap_real.band_width = i;
          score_g = tmap_sw_global_core(seq1 + 1, end_i, seq2 + 1, end_j, &ap_real, path, path_len, 0);
          if(score - prev_score == score_g) break; // TODO: is this correct?
          //if(score == score_g) break;
          if(i > j) break;
      }
      if(score > score_g) {
          tmap_error("no suitable band width", Exit, OutOfRange);
      }
      score = score_g + prev_score;
  }

  free(s_array);
  if(!_mem) free(mem);
  return score;
}
*/

static void
tmap_sw_kmp_create_table(uint8_t *seq, int32_t len, int32_t *table)
{
  int32_t cur, next;

  table[0] = -1;
  table[1] = 0;

  cur = 2;
  next = 0;
  while(cur < len) {
      if(seq[cur-1] == seq[next]) {
          table[cur] = next+1;
          cur++; next++;
      }
      else if(0 < next) {
          next = table[next];
      }
      else {
          table[cur] = 0;
          cur++;
      }
  }
}

static int32_t
tmap_sw_kmp_search(uint8_t *seq1, int32_t len1, uint8_t *seq2, uint32_t len2) 
{
  int32_t *table, i, m, r;

  table = tmap_malloc(sizeof(int32_t) * len2, "table");

  // create the table
  tmap_sw_kmp_create_table(seq2, len2, table);

  // substring match
  i = m = 0;
  r = -1;
  while(i + m < len1) {
      if(seq1[m+i] == seq2[i]) {
          i++;
          if(i == (int32_t)len2) {
              r = m;
              break;
          }
      }
      else {
          m += i - table[i];
          if(0 < i) {
              i = table[i];
          }
      }
  }

  // free
  free(table);

  return r;
}

int32_t 
tmap_sw_global_banded_core(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
                           int32_t score, tmap_sw_path_t *path, int32_t *path_len, int32_t right_j)
{
  int32_t i, j, max_bw, score_max, score_gb, len;
  int32_t *mat=NULL, *score_matrix=NULL, N_MATRIX_ROW;
  tmap_sw_param_t ap_real;
  tmap_sw_path_t *p=NULL;

  score_matrix = ap->matrix;
  N_MATRIX_ROW = ap->row;

  // Step 1: calculate the band width
  // get the maximum score
  score_max = 0; 
  for(i=0;i<N_MATRIX_ROW;i++) {
      mat = score_matrix + 0 * N_MATRIX_ROW;
      if(score_max < mat[i]) {
          score_max = mat[i];
      }
  }
  // the maximum alignment score
  score_max = (len1 < len2) ? (len1 * score_max) : (len2 * score_max);
  // the difference between the max and observed
  score_max -= score; 
  max_bw = 0;
  if(ap->gap_open + ap->gap_ext <= score_max) { // gap open allowed
      max_bw++;
      score_max -= ap->gap_open + ap->gap_ext;
      if(ap->gap_ext <= score_max) { // gap extension allowed
          max_bw += (score_max + ap->gap_ext - 1) / ap->gap_ext;
      }
  }
  
  // Step 2: if we have a perfect match, run the perfect match algorithm
  if(0 == max_bw && len1 == len2) { // implies tlen == qlen, not checked 
      int32_t best_i, best_j;
      best_i = len2;
      best_j = len1;
      if(0 < score_max) { // NB: potential for IUPAC bases
          score = 0;
          for(i=j=0;i<len1 && j<len2;i++,j++) {
              mat = score_matrix + seq2[j] * N_MATRIX_ROW;
              score += mat[seq1[i]];
          }
      }
      // retrieve the path if necessary
      if(NULL == path) {
          // do nothing
      }
      else if(NULL == path_len) {
          path[0].i = best_i; path[0].j = best_j;
      }
      else {
          i = best_i; j = best_j; p = path;
          while(0 < i && 0 < j) {
              // add to the path
              p->ctype = TMAP_SW_FROM_M;
              p->i = i; p->j = j;
              ++p;
              i--; j--;
          }
          (*path_len) = p - path;
      }
      return score;
  }

  // Step 3: run a banded alignment
  ap_real = *ap;
  ap_real.gap_end = -1;
  ap_real.band_width = ap->band_width;
  len = (len1 < len2) ? len1 : len2;
  do {
      score_gb = tmap_sw_global_core(seq1, len1, seq2, len2, &ap_real, path, path_len, right_j);
      ap_real.band_width <<= 1; // double it
  } while(score != score_gb && ap_real.band_width <= max_bw && ap_real.band_width <= len);
  // check if we need to run it at hte maximum band width
  if(score != score_gb 
     && max_bw < ap_real.band_width
     && ap_real.band_width <= (max_bw << 1) 
     && ap_real.band_width <= len) {
      ap_real.band_width = max_bw;
      score_gb = tmap_sw_global_core(seq1, len1, seq2, len2, &ap_real, path, path_len, right_j);
  }
  if(score != score_gb) {
      // NB: the vectorized smith waterman sometimes considers deletions then
      // insertions, so sometimes it can produce a better alignment score than
      // the global banded smith waterman.  So unfortunately we cannot ensure
      // that the two scores will match.  Nonetheless, this should happen
      // rarely, and even more rarely with soft clipping turned on.
      /*
      fprintf(stderr, "score=%d score_gb=%d\n", score, score_gb); 
      fprintf(stderr, "ap_real.band_width=%d\n", ap_real.band_width);
      for(i=0;i<len2;i++) {
          fputc("ACGTN"[seq2[i]], stderr);
      }
      fputc('\n', stderr);
      for(i=0;i<len1;i++) {
          fputc("ACGTN"[seq1[i]], stderr);
      }
      fputc('\n', stderr);
      tmap_bug();
      */
      return score_gb;
  }

  return score;
}

//#define TMAP_SW_CLIPPING_CORE_DEBUG 1
// TODO: optimize similar to tmap_sw_local
// - local align within a band
int32_t 
tmap_sw_clipping_core2(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
             int32_t seq1_start_skip, int32_t seq1_end_skip, int32_t seq2_start_clip, int32_t seq2_end_clip,
             tmap_sw_path_t *path, int32_t *path_len, int32_t right_justify)
{
  register int32_t i, j;

  tmap_sw_path_t *p=NULL;
  tmap_sw_dpcell_t **dpcell=NULL;
  tmap_sw_dpscore_t *curr=NULL, *last=NULL, *s=NULL;

  uint8_t ctype, ctype_next = 0;

  int32_t gap_open, gap_ext, gap_end;
  int32_t *mat=NULL, *score_matrix=NULL, N_MATRIX_ROW;

  int32_t best_i=-1, best_j=-1;
  uint8_t best_ctype=0;
  int32_t best_score = TMAP_SW_MINOR_INF;
  
#ifdef TMAP_SW_CLIPPING_CORE_DEBUG 
  char aln1[1024], aln2[1024];
  int32_t k = 0;
#endif

  gap_open = ap->gap_open;
  gap_ext = ap->gap_ext;
  gap_end = ap->gap_end;
  score_matrix = ap->matrix;
  N_MATRIX_ROW = ap->row;

#ifdef TMAP_SW_CLIPPING_CORE_DEBUG 
  fprintf(stdout, "\nHERE %s [%d,%d]\n", __func__, seq2_start_clip, seq2_end_clip);
  for(i=0;i<len2;i++) {
      fputc("ACGTN"[seq2[i]], stdout);
  }
  fputc('\n', stdout);
  for(i=0;i<len1;i++) {
      fputc("ACGTN"[seq1[i]], stdout);
  }
  fputc('\n', stdout);
#endif

  // check if seq2 is a sub-sequence of seq1
  if(len2 <= len1) {
      if(0 == right_justify) {
          for(i=j=0;j<len2;j++) { // exact from the start
              if(seq1[j] != seq2[j]) {
                  i = -1;
                  break;
              }
          }
      }
      else {
          for(i=len1-len2,j=0;j<len2;j++) { // exact from the end
              if(seq1[j+i] != seq2[j]) {
                  i = -1;
                  break;
              }
          }
      }
      if(i < 0 && 1 < len1 && 1 < len2) { // try the Knuth Morris Pratt algorithm
          if(1 == seq1_start_skip && 1 == seq1_end_skip) {
              // Here we use the Knuth-Morris-Pratt algorithms, but we could use
              // others, such as Boyer-Moore (small alphabet 2007).  KMP is easy
              // to implement
              if(0 == right_justify) {
                  i = tmap_sw_kmp_search(seq1, len1, seq2, len2);
              }
              else {
                  tmap_reverse_int(seq1, len1);
                  tmap_reverse_int(seq2, len2);
                  i = tmap_sw_kmp_search(seq1, len1, seq2, len2);
                  if(0 <= i) {
                      i = len1 - len2 - i; // zero-based
                  }
                  tmap_reverse_int(seq1, len1);
                  tmap_reverse_int(seq2, len2);
              }
          }
          else {
              i = -1;
          }
      }
      if(0 <= i) {
          best_i = i + len2; best_j = len2;
          // get the best score
          best_score = 0; 
          for(j=0;i<best_i && j<best_j;i++,j++) {
              // Note: i and j are zero-based here
              mat = score_matrix + seq1[i] * N_MATRIX_ROW;
              best_score += mat[seq2[j]];
          }
          //fprintf(stderr, "1 best_i=%d best_j=%d best_score=%d\n", best_i, best_j, best_score);
          // retrieve the path if necessary
          if(NULL == path) {
              // do nothing
          }
          else if(NULL == path_len) {
              path[0].i = best_i; path[0].j = best_j;
          }
          else {
              i = best_i; j = best_j; p = path;
              while(0 < i && 0 < j) {
                  // add to the path
                  p->ctype = TMAP_SW_FROM_M;
                  p->i = i; p->j = j;
                  ++p;
                  i--; j--;
              }
              (*path_len) = p - path;
          }
          return best_score;
      }
  }

  // allocate memory for the main cells
  dpcell = tmap_malloc(sizeof(tmap_sw_dpcell_t*) * (len1 + 1), "dpcell");
  for(i=0;i<=len1;i++) {
      dpcell[i] = tmap_malloc(sizeof(tmap_sw_dpcell_t) * (len2 + 1), "dpcell[i]");
  }
  curr = tmap_malloc(sizeof(tmap_sw_dpscore_t) * (len2 + 1), "curr");
  last = tmap_malloc(sizeof(tmap_sw_dpscore_t) * (len2 + 1), "last");

  best_score = TMAP_SW_MINOR_INF;

  // set first row
  TMAP_SW_SET_INF(curr[0]); 
  dpcell[0][0].match_from = TMAP_SW_FROM_S;
  curr[0].match_score = 0;
  if(1 == seq2_start_clip) { // start anywhere in seq2
      for(j=1;j<=len2;j++) { // for each col
          TMAP_SW_SET_INF(curr[j]);
          TMAP_SW_SET_FROM(dpcell[0][j], TMAP_SW_FROM_S);
          //dpcell[0][j].match_from = TMAP_SW_FROM_S;
          curr[j].match_score = 0;
      }
  }
  else { // start at the first base in seq2
      for(j=1;j<=len2;j++) { // for each col
          TMAP_SW_SET_INF(curr[j]);
          tmap_sw_set_end_ins(curr[j].ins_score, dpcell[0]+j, curr+j-1, right_justify);
      }
  }
  // swap curr and last
  s = curr; curr = last; last = s;

  for(i=1;i<=len1;i++) { // for each row (seq1)
      // set first column
      TMAP_SW_SET_INF(curr[0]); 
      if(1 == seq1_start_skip) { // start anywhere in seq 1
          TMAP_SW_SET_FROM(dpcell[i][0], TMAP_SW_FROM_S);
          curr[0].match_score = 0; 
      }

      mat = score_matrix + seq1[i-1] * N_MATRIX_ROW;

      for(j=1;j<=len2;j++) { // for each col (seq2)
          tmap_sw_set_match(curr[j].match_score, dpcell[i] + j,
                            last + j - 1, mat[seq2[j-1]], right_justify);
          tmap_sw_set_del(curr[j].del_score, dpcell[i] + j, last + j, right_justify);
          tmap_sw_set_ins(curr[j].ins_score, dpcell[i] + j, curr + j - 1, right_justify);
          // deal with starting anywhere in seq2
          if(1 == seq2_start_clip && curr[j].match_score < 0) {
              curr[j].match_score = 0;
              dpcell[i][j].match_from = TMAP_SW_FROM_S;
          }

          if((1 == seq1_end_skip || i == len1) && (1 == seq2_end_clip || j == len2)) { 
              if(best_score < curr[j].del_score
                 || (1 == right_justify && best_score == curr[j].del_score)) {
                  best_i = i;
                  best_j = j; 
                  best_score = curr[j].del_score;
                  best_ctype = TMAP_SW_FROM_D;
              }
              if(best_score < curr[j].ins_score
                 || (1 == right_justify && best_score == curr[j].ins_score)) {
                  best_i = i;
                  best_j = j; 
                  best_score = curr[j].ins_score;
                  best_ctype = TMAP_SW_FROM_I;
              }
              if(best_score < curr[j].match_score
                 || (1 == right_justify && best_score == curr[j].match_score)) {
                  best_i = i;
                  best_j = j; 
                  best_score = curr[j].match_score;
                  best_ctype = TMAP_SW_FROM_M;
              }
          }
      }
      // swap curr and last
      s = curr; curr = last; last = s;
  }

  // get best scoring end cell
  i = best_i; j = best_j; p = path;
  //fprintf(stderr, "2 best_i=%d best_j=%d best_score=%d\n", best_i, best_j, best_score);
  ctype = best_ctype;

  if(NULL == path) {
      // do nothing
  }
  else if(NULL == path_len) {
      path[0].i = best_i; path[0].j = best_j;
  }
  else if(best_i < 0 || best_j < 0) { // was not updated
      if(NULL != path_len) {
          (*path_len) = 0;
      }
      // do nothing
  }
  else {
      while(TMAP_SW_FROM_S != ctype
            && ((0 < i && 0 < j) || TMAP_SW_FROM_I == ctype)) {

          // get:
          // - # of read bases called from the flow
          // - the next cell type 
          // - the current offset
          switch(ctype) {
            case TMAP_SW_FROM_M:
              ctype_next = dpcell[i][j].match_from;
              break;
            case TMAP_SW_FROM_I:
              ctype_next = dpcell[i][j].ins_from;
              break;
            case TMAP_SW_FROM_D:
              ctype_next = dpcell[i][j].del_from;
              break;
            default:
              tmap_error(NULL, Exit, OutOfRange);
          }
          
          if(0 < i && 0 < j && ctype_next == TMAP_SW_FROM_S) {
              break;
          }

          // add to the path
          p->ctype = ctype;
          p->i = i; p->j = j;
          ++p;
#ifdef TMAP_SW_CLIPPING_CORE_DEBUG 
          fprintf(stdout, "[%d,%d,%d,%d]\n", i, j, ctype, ctype_next);
#endif

          // move the row and column (as necessary)
          switch(ctype) {
            case TMAP_SW_FROM_M:
#ifdef TMAP_SW_CLIPPING_CORE_DEBUG 
              aln1[k]=seq1[i-1]; aln2[k]=seq2[j-1]; k++;
#endif
              --i; --j; break;
            case TMAP_SW_FROM_I:
#ifdef TMAP_SW_CLIPPING_CORE_DEBUG 
              aln1[k]=5; aln2[k]=seq2[j-1]; k++;
#endif
              --j; break;
            case TMAP_SW_FROM_D:
#ifdef TMAP_SW_CLIPPING_CORE_DEBUG 
              aln1[k]=seq1[i-1]; aln2[k]=5; k++;
#endif
              --i; break;
            default:
              tmap_error(NULL, Exit, OutOfRange);
          }

          // move to the next cell type
          ctype = ctype_next;
      }
      (*path_len) = p - path;
#ifdef TMAP_SW_CLIPPING_CORE_DEBUG 
      fprintf(stdout, "i=%d j=%d\n", i, j);
      fprintf(stdout, "best_i=%d best_j=%d best_score=%d len1=%d len2=%d\n", best_i, best_j, best_score, len1, len2);
      fprintf(stdout, "path_len=%d best_score=%d\n", (*path_len), best_score);
#endif
  }

  // free memory for the main cells
  for(i=0;i<=len1;i++) {
      free(dpcell[i]);
  }
  free(dpcell);
  free(curr);
  free(last);
      
#ifdef TMAP_SW_CLIPPING_CORE_DEBUG 
  aln1[k]=aln2[k]='\0';
  for(i=0;i<k>>1;i++) {
      char c;
      c = aln1[i];
      aln1[i] = aln1[k-i-1];
      aln1[k-i-1] = c;
      c = aln2[i];
      aln2[i] = aln2[k-i-1];
      aln2[k-i-1] = c;
  }
  for(i=0;i<k;i++) {
      aln1[i] = "ACGTN-"[(int)aln1[i]];
      aln2[i] = "ACGTN-"[(int)aln2[i]];
  }
  fprintf(stdout, "%s\n%s\n", aln2, aln1);
  fflush(stdout);
#endif

  return best_score;
}

int32_t 
tmap_sw_extend_core(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
                    tmap_sw_path_t *path, int32_t *path_len, int32_t right_justify)
{
  return tmap_sw_clipping_core2(seq1, len1, seq2, len2, ap, 0, 1, 0, 1, path, path_len, right_justify);
}

int32_t 
tmap_sw_extend_fitting_core(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
                            tmap_sw_path_t *path, int32_t *path_len, int32_t right_justify)
{
  return tmap_sw_clipping_core2(seq1, len1, seq2, len2, ap, 0, 1, 0, 0, path, path_len, right_justify);
}

int32_t 
tmap_sw_fitting_core(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
                     tmap_sw_path_t *path, int32_t *path_len, int32_t right_justify)
{
  return tmap_sw_clipping_core2(seq1, len1, seq2, len2, ap, 1, 1, 0, 0, path, path_len, right_justify);
}

int32_t 
tmap_sw_clipping_core(uint8_t *seq1, int32_t len1, uint8_t *seq2, int32_t len2, const tmap_sw_param_t *ap,
             int32_t seq2_start_clip, int32_t seq2_end_clip,
             tmap_sw_path_t *path, int32_t *path_len, int32_t right_justify)
{
  return tmap_sw_clipping_core2(seq1, len1, seq2, len2, ap, 1, 1, seq2_start_clip, seq2_end_clip, path, path_len, right_justify);
}

uint32_t *
tmap_sw_path2cigar(const tmap_sw_path_t *path, int32_t path_len, int32_t *n_cigar)
{
  int32_t i, n;
  uint32_t *cigar;
  uint8_t last_type;

  if(path_len == 0 || path == 0) {
      *n_cigar = 0;
      return 0;
  }

  last_type = path->ctype;
  for(i = n = 1; i < path_len; ++i) {
      if(last_type != path[i].ctype) ++n;
      last_type = path[i].ctype;
  }
  *n_cigar = n;
  cigar = tmap_malloc(*n_cigar * 4, "cigar");

  TMAP_SW_CIGAR_STORE(cigar[0], path[path_len-1].ctype, 1u);
  last_type = path[path_len-1].ctype;
  for(i = path_len - 2, n = 0; i >= 0; --i) {
      if(path[i].ctype == last_type) {
          TMAP_SW_CIGAR_ADD_LENGTH(cigar[n], 1u);
      }
      else {
          TMAP_SW_CIGAR_STORE(cigar[++n], path[i].ctype, 1u);
          last_type = path[i].ctype;
      }
  }

  return cigar;
}
