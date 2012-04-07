// Coder: ngthuydiem
// Submission: 23 
// URL: http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=23054943&subnum=23

/*
     Copyright 2006, by Michael Farrar.  All rights reserved. The SWSSE2
     program and documentation may not be sold or incorporated into a
     commercial product, in whole or in part, without written consent of
     Michael Farrar.

     For further information regarding permission for use or reproduction,
     please contact Michael Farrar at:

         farrar.michael@gmail.com


     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
     CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
     TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
/* Permission was received from Michael Farrar to use the SWSSE2 code as long as 
it was released under GPL v2. He passed suddenly in late 2010.  RIP. */ 

#include <stdlib.h>
#include <stdio.h>
#include <emmintrin.h> // use SSE2
#include <cstring>
#include <sstream>
#include "Solution10.h"
using namespace std;

#define max(a, b) ((a)>(b)?a:b)
#define check_bit(var,pos) ((var) & (1<<(pos)))
#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

const int DNA_VALUE[20] = {
    0, -1,  1, -1, -1, 
    -1, 2, -1, -1, -1, 
    -1, -1, -1, -1, -1,
    -1, -1, -1, -1, 3
};

int flag = 1;
unsigned short gapOpen=0, gapExtend=0;
signed char * matrix;        
unsigned char querySeq[512], targetSeq[1024];

void buildMatrix (int match, int mismatch)
{
  matrix = (signed char *) malloc (16);
  int offset;

  for (int i = 0; i < 4; ++i) {
      offset = (i << 2);
      for (int j = 0; j < 4; ++j) 
        matrix[offset+j] = (j==i) ? match : mismatch;
  }                
}

typedef struct {
    __m128i        *pvbQueryProf;
    __m128i        *pvsQueryProf;
    __m128i        *pvH1;
    __m128i        *pvH2;
    __m128i        *pvE;
    unsigned char  *pData;
    unsigned short  bias;
} SwStripedData;

SwStripedData * swStripedInitByte(unsigned char   *querySeq,
                                  int              queryLength,
                                  signed char *    matrix)
{
  int i, j, k;

  int segSize;
  int nCount;

  int bias;

  int lenQryByte;
  int lenQryShort;

  int weight;

  //short *ps;
  char *pc;

  signed char *matrixRow;

  size_t aligned;

  SwStripedData *pSwData;

  lenQryByte = (queryLength + 15) >> 4;
  lenQryShort = (queryLength + 7) >> 3;

  pSwData = (SwStripedData *) malloc (sizeof (SwStripedData));
  if (!pSwData) {
      fprintf (stderr, "Unable to allocate memory for SW data\n");
      exit (-1);
  }

  nCount = 64 +                             /* slack bytes */
    (lenQryByte << 2) +        /* query profile byte */
    (lenQryShort << 2) +       /* query profile short */
    (lenQryShort * 3);               /* vH1, vH2 and vE */

  pSwData->pData = (unsigned char *) calloc (nCount, sizeof (__m128i));
  if (!pSwData->pData) {
      fprintf (stderr, "Unable to allocate memory for SW data buffers\n");
      exit (-1);
  }

  /* since we might port this to another platform, lets align the data */
  /* to 16 byte boundries ourselves */
  aligned = ((size_t) pSwData->pData + 15) & ~(0x0f);

  pSwData->pvbQueryProf = (__m128i *) aligned;
  pSwData->pvsQueryProf = pSwData->pvbQueryProf + (lenQryByte << 2);

  pSwData->pvH1 = pSwData->pvsQueryProf + (lenQryShort << 2);
  pSwData->pvH2 = pSwData->pvH1 + lenQryShort;
  pSwData->pvE  = pSwData->pvH2 + lenQryShort;

  /* Find the bias to use in the substitution matrix */
  bias = matrix[1]; // bias is the mismatch score

  /* Fill in the byte query profile */
  pc = (char *) pSwData->pvbQueryProf;
  segSize = (queryLength + 15) >> 4;
  nCount = (segSize << 4);
  for (i = 0; i < 4; ++i) {
      matrixRow = matrix + (i << 2);
      for (j = 0; j < segSize; ++j) {
          for (k = j; k < nCount; k += segSize) {
              if (k >= queryLength) 
                weight = 0;
              else 
                weight = matrixRow[*(querySeq + k)];

              *pc++ = (char) (weight - bias);
          }
      }
  }

  pSwData->bias = (unsigned short) -bias;

  return pSwData;
}

SwStripedData * swStripedInitWord(unsigned char   *querySeq,
                                  int              queryLength,
                                  signed char *    matrix)
{
  int i, j, k;

  int segSize;
  int nCount;
  int lenQryShort;
  int weight;

  short *ps;    
  signed char *matrixRow;

  size_t aligned;

  SwStripedData *pSwData;

  lenQryShort = (queryLength + 7) >> 3;

  pSwData = (SwStripedData *) malloc (sizeof (SwStripedData));
  if (!pSwData) {
      fprintf (stderr, "Unable to allocate memory for SW data\n");
      exit (-1);
  }

  nCount = 64 +                             /* slack bytes */
    (lenQryShort << 2) +       /* query profile short */
    (lenQryShort * 3);               /* vH1, vH2 and vE */

  pSwData->pData = (unsigned char *) calloc (nCount, sizeof (__m128i));
  if (!pSwData->pData) {
      fprintf (stderr, "Unable to allocate memory for SW data buffers\n");
      exit (-1);
  }

  /* since we might port this to another platform, lets align the data */
  /* to 16 byte boundries ourselves */
  aligned = ((size_t) pSwData->pData + 15) & ~(0x0f);

  pSwData->pvsQueryProf = (__m128i *) aligned;
  pSwData->pvH1 = pSwData->pvsQueryProf + (lenQryShort << 2);
  pSwData->pvH2 = pSwData->pvH1 + lenQryShort;
  pSwData->pvE  = pSwData->pvH2 + lenQryShort;

  /* Fill in the short query profile */
  ps = (short *) pSwData->pvsQueryProf;
  segSize = ((queryLength + 7) >> 3);
  nCount = (segSize << 3);
  for (i = 0; i < 4; ++i) {
      matrixRow = (matrix + (i << 2));
      for (j = 0; j < segSize; ++j) {
          for (k = j; k < nCount; k += segSize) {
              if (k >= queryLength) 
                weight = 0;
              else 
                weight = matrixRow[*(querySeq + k)];

              *ps++ = (unsigned short) weight;
          }
      }
  }

  return pSwData;
}

void swStripedInitWordAfterByte( SwStripedData *pSwData,
                                unsigned char   *querySeq,
                                int              queryLength,
                                signed char *    matrix)
{
  int i, j, k;

  int segSize, nCount, weight;

  short *ps;

  signed char *matrixRow;

  /* Fill in the short query profile */
  ps = (short *) pSwData->pvsQueryProf;
  segSize = ((queryLength + 7) >> 3);
  nCount = (segSize << 3);
  for (i = 0; i < 4; ++i) {
      matrixRow = (matrix + (i << 2));
      for (j = 0; j < segSize; ++j) {
          for (k = j; k < nCount; k += segSize) {
              if (k >= queryLength) 
                weight = 0;
              else 
                weight = matrixRow[*(querySeq + k)];

              *ps++ = (unsigned short) weight;
          }
      }
  }
}

int swStripedByte3(int              queryLength,
                   unsigned char   *dbSeq,
                   int              dbLength,
                   unsigned short   gapOpen,
                   unsigned short   gapExtend,
                   __m128i         *pvQueryProf,
                   __m128i         *pvHLoad,
                   __m128i         *pvHStore,
                   __m128i         *pvE,
                   unsigned short   bias,
                   int                dir,
                   int *_opt, int *_te, int *_qe, int *_n_best)
{
  int     i, j, temp;    
  int     score = 0;

  int     n_best = 0, query_end = -1, target_end = -1;
  int     greaterMask, equalMask, cmp;
  int     iter = (queryLength + 15) >> 4;    

  __m128i *pv;

  __m128i vE, vF, vH;

  __m128i vMaxScore;
  __m128i vBias;
  __m128i vGapOpen;
  __m128i vGapExtend;

  __m128i vTemp;
  __m128i vZero = _mm_setzero_si128();

  __m128i *pvScore;

  vBias = _mm_set1_epi8(bias);
  vGapOpen = _mm_set1_epi8(gapOpen);
  vGapExtend = _mm_set1_epi8(gapExtend);
  vMaxScore = vZero;

  /* Zero out the storage vector */
  for (i = 0; i < iter; ++i)
    {
      _mm_store_si128 (pvE + i, vZero);
      _mm_store_si128 (pvHStore + i, vZero);
    }

  for (i = 0; i < dbLength; ++i)
    {
      /* fetch first data asap. */
      pvScore = pvQueryProf + dbSeq[i] * iter;

      /* zero out F. */
      vF = vZero;

      /* load the next h value */
      vH = _mm_load_si128 (pvHStore + iter - 1);
      vH = _mm_slli_si128 (vH, 1);

      pv = pvHLoad;
      pvHLoad = pvHStore;
      pvHStore = pv;

      for (j = 0; j < iter; ++j)
        {
          /* load values of vF and vH from previous row (one unit up) */
          vE = _mm_load_si128 (pvE + j);

          /* add score to vH */
          vH = _mm_adds_epu8 (vH, *(pvScore + j));
          vH = _mm_subs_epu8 (vH, vBias);

          /* Update highest score encountered this far */               
          vTemp = _mm_cmpeq_epi8 (vH, vMaxScore);
          equalMask = _mm_movemask_epi8 (vTemp);

          vTemp = _mm_subs_epu8 (vH, vMaxScore);
          vTemp = _mm_cmpgt_epi8 (vTemp, vZero);            
          greaterMask = _mm_movemask_epi8(vTemp);

          if (greaterMask != 0x0000)
            {                
              vMaxScore = vH;
              /* find largest score in the vMaxScore vector */
              vTemp = _mm_srli_si128 (vMaxScore, 8);
              vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
              vTemp = _mm_srli_si128 (vMaxScore, 4);
              vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
              vTemp = _mm_srli_si128 (vMaxScore, 2);
              vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);
              vTemp = _mm_srli_si128 (vMaxScore, 1);
              vMaxScore = _mm_max_epu8 (vMaxScore, vTemp);

              /* store in temporary variable */
              score = _mm_extract_epi16 (vMaxScore, 0);
              score = score & 0x00ff;                        
              if (score + bias >= 255)
                return 255;

              vMaxScore = _mm_set1_epi8(score);

              vTemp = _mm_cmpeq_epi8 (vH, vMaxScore);
              equalMask = _mm_movemask_epi8 (vTemp);

              target_end = i;
              n_best = 0;
              query_end = -1;
              for (int k = 0; k <= 15; ++k)
                {                    
                  if (check_bit(equalMask, k)) {
                      temp = k * iter + j;
                      if (temp < queryLength) {
                          ++n_best;
                          if ( dir == 0 && query_end == -1) 
                            query_end = temp;            
                          if (dir == 1)
                            query_end = temp;            
                      }
                  }
                }                        
            } 
          else if (equalMask != 0x0000)
            {                
              for (int k = 0; k <= 15; ++k)
                {
                  if (check_bit(equalMask, k)) {
                      temp = k * iter + j;
                      if (temp < queryLength) {
                          ++n_best;                                        
                          if (dir == 1 && temp >= query_end) {
                              query_end = temp;
                              target_end = i;
                          }
                      }
                  }
                }
            }                                  
          /* get max from vH, vE and vF */
          vH = _mm_max_epu8 (vH, vE);
          vH = _mm_max_epu8 (vH, vF);

          /* save vH values */
          _mm_store_si128 (pvHStore + j, vH);

          /* update vE value */
          vH = _mm_subs_epu8 (vH, vGapOpen);
          vE = _mm_subs_epu8 (vE, vGapExtend);
          vE = _mm_max_epu8 (vE, vH);

          /* update vF value */
          vF = _mm_subs_epu8 (vF, vGapExtend);
          vF = _mm_max_epu8 (vF, vH);

          /* save vE values */
          _mm_store_si128 (pvE + j, vE);

          /* load the next h value */
          vH = _mm_load_si128 (pvHLoad + j);
        }

      /* reset pointers to the start of the saved data */
      j = 0;
      vH = _mm_load_si128 (pvHStore + j);

      /*  the computed vF value is for the given column.  since */
      /*  we are at the end, we need to shift the vF value over */
      /*  to the next column. */
      vF = _mm_slli_si128 (vF, 1);
      vTemp = _mm_subs_epu8 (vH, vGapOpen);
      vTemp = _mm_subs_epu8 (vF, vTemp);
      vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
      cmp  = _mm_movemask_epi8 (vTemp);

      while (cmp != 0xffff) 
        {
          vE = _mm_load_si128 (pvE + j);

          vH = _mm_max_epu8 (vH, vF);

          /* save vH values */
          _mm_store_si128 (pvHStore + j, vH);

          /*  update vE incase the new vH value would change it */
          vH = _mm_subs_epu8 (vH, vGapOpen);
          vE = _mm_max_epu8 (vE, vH);
          _mm_store_si128 (pvE + j, vE);

          /* update vF value */
          vF = _mm_subs_epu8 (vF, vGapExtend);

          ++j;
          if (j >= iter)
            {
              j = 0;
              vF = _mm_slli_si128 (vF, 1);
            }

          vH = _mm_load_si128 (pvHStore + j);

          vTemp = _mm_subs_epu8 (vH, vGapOpen);
          vTemp = _mm_subs_epu8 (vF, vTemp);
          vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
          cmp  = _mm_movemask_epi8 (vTemp);
        }
    }        

  (*_opt) = score;
  (*_te) = target_end;
  (*_qe) = query_end;
  (*_n_best) = n_best;

  /* return largest score */
  return score;
}

void swStripedWord3(int              queryLength,
                    unsigned char   *dbSeq,
                    int              dbLength,
                    unsigned short   gapOpen,
                    unsigned short   gapExtend,
                    __m128i         *pvQueryProf,
                    __m128i         *pvHLoad,
                    __m128i         *pvHStore,
                    __m128i         *pvE,
                    int               dir,
                    int *_opt, int *_te, int *_qe, int *_n_best)
{
  int     i, j, temp;
  int     score = 0;

  int     cmp;
  int     iter = (queryLength + 7) >> 3;

  __m128i *pv;

  __m128i vE, vF, vH;

  __m128i vMaxScore, vGapOpen, vGapExtend;

  __m128i vMin;
  __m128i vTemp;
  __m128i *pvScore;
  __m128i vZero = _mm_setzero_si128();
  __m128i vInf = _mm_set1_epi16(0x8000);

  int equalMask, greaterMask;
  int n_best = 0, query_end = -1, target_end = -1;

  /* Load gap opening penalty to all elements of a constant */
  vGapOpen = _mm_set1_epi16 (gapOpen);

  /* Load gap extension penalty to all elements of a constant */
  vGapExtend = _mm_set1_epi16 (gapExtend);

  /*  load vMaxScore with the vZero.  since we are using signed */
  /*  math, we will bias the maxscore to -32768 so we have the */
  /*  full range of the short. */

  vMin = _mm_insert_epi16 (vZero, 0x8000, 0);
  vMaxScore = vInf;

  /* Zero out the storage vector */
  for (i = 0; i < iter; ++i)
    {
      _mm_store_si128 (pvE + i, vMaxScore);
      _mm_store_si128 (pvHStore + i, vMaxScore);
    }

  for (i = 0; i < dbLength; ++i)
    {
      /* fetch first data asap. */
      pvScore = pvQueryProf + dbSeq[i] * iter;

      /* zero out F. */
      vF = vInf;

      /* load the next h value */
      vH = _mm_load_si128 (pvHStore + iter - 1);
      vH = _mm_slli_si128 (vH, 2);
      vH = _mm_or_si128 (vH, vMin);

      pv = pvHLoad;
      pvHLoad = pvHStore;
      pvHStore = pv;

      for (j = 0; j < iter; ++j)
        {
          /* load values of vF and vH from previous row (one unit up) */
          vE = _mm_load_si128 (pvE + j);

          /* add score to vH */
          vH = _mm_adds_epi16 (vH, *pvScore++);

          /* Update highest score encountered this far */           
          vTemp = _mm_cmpeq_epi16(vH, vMaxScore);
          equalMask = _mm_movemask_epi8(vTemp);
          vTemp = _mm_cmpgt_epi16(vH, vMaxScore);
          greaterMask = _mm_movemask_epi8(vTemp);

          if (greaterMask != 0x0000)
            {
              vMaxScore = vH;
              // find largest score in the vMaxScore vector 
              vTemp = _mm_srli_si128 (vMaxScore, 8);
              vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
              vTemp = _mm_srli_si128 (vMaxScore, 4);
              vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
              vTemp = _mm_srli_si128 (vMaxScore, 2);
              vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);

              // store in temporary variable 
              score = (short)_mm_extract_epi16 (vMaxScore, 0);
              vMaxScore = _mm_set1_epi16(score);

              vTemp = _mm_cmpeq_epi16(vH, vMaxScore);
              equalMask = _mm_movemask_epi8(vTemp);

              target_end = i;
              n_best = 0;
              query_end = -1;
              for (int k = 0; k <= 7; ++k)
                {                    
                  if (check_bit(equalMask, k+k)) {
                      temp = k * iter + j;
                      if (temp < queryLength) {
                          ++n_best;
                          if ( dir == 0 && query_end == -1) 
                            query_end = temp;            
                          if (dir == 1)
                            query_end = temp;            
                      }
                  }
                }            
            } 
          else if (equalMask != 0x0000)
            {                
              for (int k = 0; k <= 7; ++k)
                {
                  if (check_bit(equalMask, k+k)) {
                      temp = k * iter + j;
                      if (temp < queryLength) {
                          ++n_best;                                        
                          if (dir == 1 && temp >= query_end) {
                              query_end = temp;
                              target_end = i;
                          }
                      }
                  }
                }
            }                                                                                

          /* get max from vH, vE and vF */ 
          vH = _mm_max_epi16 (vH, vE);
          vH = _mm_max_epi16 (vH, vF);                    

          /* save vH values */
          _mm_store_si128 (pvHStore + j, vH);        

          /* update vE value */
          vH = _mm_subs_epi16 (vH, vGapOpen);
          vE = _mm_subs_epi16 (vE, vGapExtend);
          vE = _mm_max_epi16 (vE, vH);

          /* update vF value */
          vF = _mm_subs_epi16 (vF, vGapExtend);
          vF = _mm_max_epi16 (vF, vH);

          /* save vE values */
          _mm_store_si128 (pvE + j, vE);

          /* load the next h value */
          vH = _mm_load_si128 (pvHLoad + j);
        }

      /* reset pointers to the start of the saved data */
      j = 0;
      vH = _mm_load_si128 (pvHStore + j);

      /*  the computed vF value is for the given column.  since */
      /*  we are at the end, we need to shift the vF value over */
      /*  to the next column. */
      vF = _mm_slli_si128 (vF, 2);
      vF = _mm_or_si128 (vF, vMin);
      vTemp = _mm_subs_epi16 (vH, vGapOpen);
      vTemp = _mm_cmpgt_epi16 (vF, vTemp);
      cmp  = _mm_movemask_epi8 (vTemp);
      while (cmp != 0x0000) 
        {                               
          vH = _mm_max_epi16 (vH, vF);

          /* save vH values */
          _mm_store_si128 (pvHStore + j, vH);

          /*  update vE incase the new vH value would change it */
          vH = _mm_subs_epi16 (vH, vGapOpen);

          vE = _mm_load_si128 (pvE + j);     
          vE = _mm_max_epi16 (vE, vH);
          _mm_store_si128 (pvE + j, vE);

          /* update vF value */
          vF = _mm_subs_epi16 (vF, vGapExtend);

          ++j;
          if (j >= iter)
            {
              j = 0;
              vF = _mm_slli_si128 (vF, 2);
              vF = _mm_or_si128 (vF, vMin);
            }

          vH = _mm_load_si128 (pvHStore + j);

          vTemp = _mm_subs_epi16 (vH, vGapOpen);
          vTemp = _mm_cmpgt_epi16 (vF, vTemp);
          cmp  = _mm_movemask_epi8 (vTemp);
        }
    }   
  (*_opt) = score+0x8000;
  (*_te) = target_end;
  (*_qe) = query_end;
  (*_n_best) = n_best;
}

SwStripedData *        swData3;

void swStripedScan3 (unsigned char* querySeq,
                     int              queryLength,
                     unsigned char   *targetSeq,
                     int                targetLength,    
                     unsigned short    gapOpen,
                     unsigned short    gapExtend,        
                     int                dir,
                     int *opt, int *te, int *qe, int *n_best,
                     int            buildProfile)
{
  int score;              

  if (buildProfile == 1)    {
      if (swData3)
        {
          free(swData3->pData);
          free(swData3);
        }
      swData3 = swStripedInitByte(querySeq, queryLength, matrix);
  }
  score = swStripedByte3 (queryLength, 
                          targetSeq, targetLength,
                          gapOpen, gapExtend,
                          swData3->pvbQueryProf,
                          swData3->pvH1,
                          swData3->pvH2,
                          swData3->pvE,
                          swData3->bias,
                          dir, opt, te, qe, n_best);

  /* check if needs a run with higher precision */
  if (score >= 255) 
    {        
      swStripedInitWordAfterByte(swData3, querySeq, queryLength, matrix);
      swStripedWord3 (queryLength, 
                      targetSeq, targetLength,
                      gapOpen, gapExtend,
                      swData3->pvsQueryProf,
                      swData3->pvH1,
                      swData3->pvH2,
                      swData3->pvE,
                      dir, opt, te, qe, n_best);
    }           
}

void swStripedWord2(int              queryLength,
                    unsigned char   *dbSeq,
                    int              dbLength,
                    unsigned short   gapOpen,
                    unsigned short   gapExtend,
                    __m128i         *pvQueryProf,
                    __m128i         *pvHLoad,
                    __m128i         *pvHStore,
                    __m128i         *pvE,
                    int               dir,
                    int *_opt, int *_te, int *_qe, int *_n_best)
{
  int     i, j;
  int     score = 0;

  int     cmp;
  int     iter = (queryLength + 7) >> 3;

  __m128i *pv;

  __m128i vE, vF, vH;

  __m128i vMaxScore, vGapOpen, vGapExtend;

  __m128i vMin, vTemp;
  __m128i *pvScore;
  __m128i vZero = _mm_setzero_si128();
  __m128i vInf = _mm_set1_epi16(0x8000);

  int n_best = 0, query_end = -1;

  /* Load gap opening penalty to all elements of a constant */
  vGapOpen = _mm_set1_epi16 (gapOpen);

  /* Load gap extension penalty to all elements of a constant */
  vGapExtend = _mm_set1_epi16 (gapExtend);

  /*  load vMaxScore with the vZero.  since we are using signed */
  /*  math, we will bias the maxscore to -32768 so we have the */
  /*  full range of the short. */

  vMin = _mm_insert_epi16 (vZero, 0x8000, 0);
  vMaxScore = vInf;    

  /* Zero out the storage vector */
  for (i = 0; i < iter; ++i)
    {
      _mm_store_si128 (pvE + i, vInf);
      _mm_store_si128 (pvHStore + i, vInf);
    }

  for (i = 0; i < dbLength; ++i)
    {
      /* fetch first data asap. */
      pvScore = pvQueryProf + dbSeq[i] * iter;

      /* zero out F. */
      vF = vInf;

      /* load the next h value */
      vH = _mm_load_si128 (pvHStore + iter - 1);
      vH = _mm_slli_si128 (vH, 2);
      vH = _mm_or_si128 (vH, vMin);

      pv = pvHLoad;
      pvHLoad = pvHStore;
      pvHStore = pv;

      for (j = 0; j < iter; ++j)
        {
          /* load values of vF and vH from previous row (one unit up) */
          vE = _mm_load_si128 (pvE + j);

          /* add score to vH */
          vH = _mm_adds_epi16 (vH, *pvScore++);           

          /* get max from vH, vE and vF */ 
          vH = _mm_max_epi16 (vH, vE);
          vH = _mm_max_epi16 (vH, vF);                 

          /* save vH values */
          _mm_store_si128 (pvHStore + j, vH);        

          /* update vE value */
          vH = _mm_subs_epi16 (vH, vGapOpen);
          vE = _mm_subs_epi16 (vE, vGapExtend);
          vE = _mm_max_epi16 (vE, vH);

          /* update vF value */
          vF = _mm_subs_epi16 (vF, vGapExtend);
          vF = _mm_max_epi16 (vF, vH);

          /* save vE values */
          _mm_store_si128 (pvE + j, vE);

          /* load the next h value */
          vH = _mm_load_si128 (pvHLoad + j);
        }

      /* reset pointers to the start of the saved data */
      j = 0;
      vH = _mm_load_si128 (pvHStore + j);

      /*  the computed vF value is for the given column.  since */
      /*  we are at the end, we need to shift the vF value over */
      /*  to the next column. */
      vF = _mm_slli_si128 (vF, 2);
      vF = _mm_or_si128 (vF, vMin);
      vTemp = _mm_subs_epi16 (vH, vGapOpen);
      vTemp = _mm_cmpgt_epi16 (vF, vTemp);
      cmp  = _mm_movemask_epi8 (vTemp);
      while (cmp != 0x0000) 
        {                               
          vH = _mm_max_epi16 (vH, vF);

          /* save vH values */
          _mm_store_si128 (pvHStore + j, vH);

          /*  update vE incase the new vH value would change it */
          vH = _mm_subs_epi16 (vH, vGapOpen);

          vE = _mm_load_si128 (pvE + j);     
          vE = _mm_max_epi16 (vE, vH);
          _mm_store_si128 (pvE + j, vE);

          /* update vF value */
          vF = _mm_subs_epi16 (vF, vGapExtend);

          ++j;
          if (j >= iter)
            {
              j = 0;
              vF = _mm_slli_si128 (vF, 2);
              vF = _mm_or_si128 (vF, vMin);
            }

          vH = _mm_load_si128 (pvHStore + j);

          vTemp = _mm_subs_epi16 (vH, vGapOpen);
          vTemp = _mm_cmpgt_epi16 (vF, vTemp);
          cmp  = _mm_movemask_epi8 (vTemp);
        }
    }   

  short * pH = (short*) pvHStore;

  score = -0x7fff;        
  int offset = 0, counter = 0;
  short temp;
  for (i = 0; i < iter; ++i)     
    {    
      counter = 0;
      for (j = i; j < queryLength; j += iter)     
        {                
          temp = pH[offset + counter];
          if (score == temp) 
            {
              ++n_best;
              if (dir == 1 && j >= query_end) 
                query_end = j;
            }
          else if (score < temp) 
            {                        
              score = temp;
              n_best = 1;
              query_end = j;
            }     

          ++counter;
        }
      offset += 8;
    }

  /* return largest score */     
  (*_opt) = score+0x8000;
  (*_te) = dbLength-1;
  (*_qe) = query_end;
  (*_n_best) = n_best;
}


void swStripedWord1(int              queryLength,
                    unsigned char   *dbSeq,
                    int              dbLength,
                    unsigned short   gapOpen,
                    unsigned short   gapExtend,
                    __m128i         *pvQueryProf,
                    __m128i         *pvHLoad,
                    __m128i         *pvHStore,
                    __m128i         *pvE,
                    int               dir,
                    int *_opt, int *_te, int *_qe, int *_n_best)
{
  int     i, j, temp;
  int     score = 0;

  int     cmp;
  int     iter = (queryLength + 7) >> 3;

  __m128i *pv;

  __m128i vE, vF, vH;

  __m128i vMaxScore, vGapOpen, vGapExtend;

  __m128i vMin;
  __m128i vTemp;
  __m128i *pvScore;
  __m128i vZero = _mm_setzero_si128();
  __m128i vInf = _mm_set1_epi16(0x8000);

  int equalMask, greaterMask;
  int n_best = 0, query_end = -1, target_end = -1;

  /* Load gap opening penalty to all elements of a constant */
  vGapOpen = _mm_set1_epi16 (gapOpen);

  /* Load gap extension penalty to all elements of a constant */
  vGapExtend = _mm_set1_epi16 (gapExtend);

  /*  load vMaxScore with the vZero.  since we are using signed */
  /*  math, we will bias the maxscore to -32768 so we have the */
  /*  full range of the short. */

  vMin = _mm_insert_epi16 (vZero, 0x8000, 0);
  vMaxScore = vInf;    

  /* Zero out the storage vector */
  for (i = 0; i < iter; ++i) {        
      _mm_store_si128 (pvHStore + i, vZero);
      _mm_store_si128 (pvE + i, vInf);
  }    

  for (i = 0; i < dbLength; ++i)
    {
      /* fetch first data asap. */
      pvScore = pvQueryProf + dbSeq[i] * iter;

      /* zero out F. */
      vF = vInf;            

      /* load the next h value */
      vH = _mm_load_si128 (pvHStore + iter - 1);
      vH = _mm_slli_si128 (vH, 2);    
      if (i > 0)
        vH = _mm_insert_epi16 (vH, -gapOpen - (i-1) * gapExtend, 0);

      pv = pvHLoad;
      pvHLoad = pvHStore;
      pvHStore = pv;

      for (j = 0; j < iter; ++j)
        {
          /* load values of vF and vH from previous row (one unit up) */
          vE = _mm_load_si128 (pvE + j);

          if (i == 0)                 
            vE = _mm_subs_epi16(vZero, vGapOpen);
          else
            vE = _mm_load_si128 (pvE + j);

          /* add score to vH */
          vH = _mm_adds_epi16 (vH, *pvScore++);  

          /* get max from vH, vE and vF */ 
          vH = _mm_max_epi16 (vH, vE);
          vH = _mm_max_epi16 (vH, vF);                    

          /* save vH values */
          _mm_store_si128 (pvHStore + j, vH);        

          /* Update  highest score encountered this far */           
          vTemp = _mm_cmpeq_epi16(vH, vMaxScore);
          equalMask = _mm_movemask_epi8(vTemp);
          vTemp = _mm_cmpgt_epi16(vH, vMaxScore);
          greaterMask = _mm_movemask_epi8(vTemp);

          if (greaterMask != 0x0000)
            {
              vMaxScore = vH;
              // find largest score in the vMaxScore vector 
              vTemp = _mm_srli_si128 (vMaxScore, 8);
              vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
              vTemp = _mm_srli_si128 (vMaxScore, 4);
              vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);
              vTemp = _mm_srli_si128 (vMaxScore, 2);
              vMaxScore = _mm_max_epi16 (vMaxScore, vTemp);

              // store in temporary variable 
              score = (short)_mm_extract_epi16 (vMaxScore, 0);
              vMaxScore = _mm_set1_epi16(score);

              vTemp = _mm_cmpeq_epi16(vH, vMaxScore);
              equalMask = _mm_movemask_epi8(vTemp);

              target_end = i;
              n_best = 0;
              query_end = -1;
              for (int k = 0; k <= 7; ++k)
                {                    
                  if (check_bit(equalMask, k+k)) {
                      temp = k * iter + j;
                      if (temp < queryLength) {
                          ++n_best;
                          if ( dir == 0 && query_end == -1) 
                            query_end = temp;            
                          if (dir == 1)
                            query_end = temp;            
                      }
                  }
                }                
            } 
          else if (equalMask != 0x0000)
            {                
              for (int k = 0; k <= 7; ++k)
                {
                  if (check_bit(equalMask, k+k)) {
                      temp = k * iter + j;
                      if (temp < queryLength) {
                          ++n_best;                                        
                          if (dir == 1 && temp >= query_end) {
                              query_end = temp;
                              target_end = i;
                          }
                      }
                  }
                }
            }                                                                                            

          /* update vE value */
          vH = _mm_subs_epi16 (vH, vGapOpen);
          vE = _mm_subs_epi16 (vE, vGapExtend);
          vE = _mm_max_epi16 (vE, vH);

          /* update vF value */
          vF = _mm_subs_epi16 (vF, vGapExtend);
          vF = _mm_max_epi16 (vF, vH);

          /* save vE values */
          _mm_store_si128 (pvE + j, vE);

          /* load the next h value */
          vH = _mm_load_si128 (pvHLoad + j);
        }

      /* reset pointers to the start of the saved data */
      j = 0;
      vH = _mm_load_si128 (pvHStore + j);

      /*  the computed vF value is for the given column.  since */
      /*  we are at the end, we need to shift the vF value over */
      /*  to the next column. */        

      vF = _mm_slli_si128 (vF, 2);
      vF = _mm_or_si128 (vF, vMin);

      vTemp = _mm_subs_epi16 (vH, vGapOpen);
      vTemp = _mm_cmpgt_epi16 (vF, vTemp);
      cmp  = _mm_movemask_epi8 (vTemp);
      while (cmp != 0x0000) 
        {                                                           
          vH = _mm_max_epi16 (vH, vF);

          /* save vH values */
          _mm_store_si128 (pvHStore + j, vH);

          /*  update vE incase the new vH value would change it */
          vH = _mm_subs_epi16 (vH, vGapOpen);

          vE = _mm_load_si128 (pvE + j);     
          vE = _mm_max_epi16 (vE, vH);
          _mm_store_si128 (pvE + j, vE);

          /* update vF value */
          vF = _mm_subs_epi16 (vF, vGapExtend);

          ++j;
          if (j >= iter)
            {
              j = 0;            
              vF = _mm_slli_si128 (vF, 2);
              vF = _mm_or_si128 (vF, vMin);
            }            

          vH = _mm_load_si128 (pvHStore + j);

          vTemp = _mm_subs_epi16 (vH, vGapOpen);
          vTemp = _mm_cmpgt_epi16 (vF, vTemp);
          cmp  = _mm_movemask_epi8 (vTemp);
        }
    }   

  (*_opt) = score;
  (*_te) = target_end;
  (*_qe) = query_end;
  (*_n_best) = n_best;
}

void swStripedWord0(int              queryLength,
                    unsigned char   *dbSeq,
                    int              dbLength,
                    unsigned short   gapOpen,
                    unsigned short   gapExtend,
                    __m128i         *pvQueryProf,
                    __m128i         *pvHLoad,
                    __m128i         *pvHStore,
                    __m128i         *pvE,
                    int               dir,
                    int *_opt, int *_te, int *_qe, int *_n_best)
{
  int     i, j;
  int     score = 0;

  int     cmp;
  int     iter = (queryLength + 7) >> 3;

  __m128i *pv;

  __m128i vE, vF, vH;

  __m128i vMaxScore, vGapOpen, vGapExtend;

  __m128i vMin, vTemp;
  __m128i *pvScore;
  __m128i vZero = _mm_setzero_si128();
  __m128i vInf = _mm_set1_epi16(0x8000);

  int n_best = 0, query_end = -1;

  /* Load gap opening penalty to all elements of a constant */
  vGapOpen = _mm_set1_epi16 (gapOpen);

  /* Load gap extension penalty to all elements of a constant */
  vGapExtend = _mm_set1_epi16 (gapExtend);

  /*  load vMaxScore with the vZero.  since we are using signed */
  /*  math, we will bias the maxscore to -32768 so we have the */
  /*  full range of the short. */

  vMin = _mm_insert_epi16 (vZero, 0x8000, 0);
  vMaxScore = vInf;    

  /* Zero out the storage vector */        
  for (i = 0; i < iter; ++i) {                
      _mm_store_si128 (pvHStore + i, vZero);
      _mm_store_si128 (pvE + i, vInf);
  }    

  for (i = 0; i < dbLength; ++i)
    {
      /* fetch first data asap. */
      pvScore = pvQueryProf + dbSeq[i] * iter;

      /* zero out F. */
      vF = vInf;    

      /* load the next h value */
      vH = _mm_load_si128 (pvHStore + iter - 1);
      vH = _mm_slli_si128 (vH, 2);    
      if (i > 0)
        vH = _mm_insert_epi16 (vH, -gapOpen - (i-1) * gapExtend, 0);

      pv = pvHLoad;
      pvHLoad = pvHStore;
      pvHStore = pv;

      for (j = 0; j < iter; ++j)
        {
          /* load values of vF and vH from previous row (one unit up) */                        
          if (i == 0)                 
            vE = _mm_subs_epi16(vZero, vGapOpen);
          else
            vE = _mm_load_si128 (pvE + j);

          /* add score to vH */
          vH = _mm_adds_epi16 (vH, *pvScore++);           

          /* get max from vH, vE and vF */ 
          vH = _mm_max_epi16 (vH, vE);
          vH = _mm_max_epi16 (vH, vF);                 

          /* save vH values */
          _mm_store_si128 (pvHStore + j, vH);        

          /* update vE value */
          vH = _mm_subs_epi16 (vH, vGapOpen);
          vE = _mm_subs_epi16 (vE, vGapExtend);
          vE = _mm_max_epi16 (vE, vH);

          /* update vF value */
          vF = _mm_subs_epi16 (vF, vGapExtend);
          vF = _mm_max_epi16 (vF, vH);

          /* save vE values */
          _mm_store_si128 (pvE + j, vE);

          /* load the next h value */
          vH = _mm_load_si128 (pvHLoad + j);
        }

      // reset pointers to the start of the saved data 
      j = 0;
      vH = _mm_load_si128 (pvHStore + j);

      //  the computed vF value is for the given column.  since 
      //  we are at the end, we need to shift the vF value over 
      //  to the next column. 
      vF = _mm_slli_si128 (vF, 2);
      vF = _mm_or_si128 (vF, vMin);
      vTemp = _mm_subs_epi16 (vH, vGapOpen);
      vTemp = _mm_cmpgt_epi16 (vF, vTemp);
      cmp  = _mm_movemask_epi8 (vTemp);
      while (cmp != 0x0000) 
        {                               
          vH = _mm_max_epi16 (vH, vF);

          // save vH values 
          _mm_store_si128 (pvHStore + j, vH);

          //  update vE incase the new vH value would change it 
          vH = _mm_subs_epi16 (vH, vGapOpen);

          vE = _mm_load_si128 (pvE + j);     
          vE = _mm_max_epi16 (vE, vH);
          _mm_store_si128 (pvE + j, vE);

          // update vF value 
          vF = _mm_subs_epi16 (vF, vGapExtend);

          ++j;
          if (j >= iter)
            {
              j = 0;
              vF = _mm_slli_si128 (vF, 2);
              vF = _mm_or_si128 (vF, vMin);            
            }

          vH = _mm_load_si128 (pvHStore + j);

          vTemp = _mm_subs_epi16 (vH, vGapOpen);
          vTemp = _mm_cmpgt_epi16 (vF, vTemp);
          cmp  = _mm_movemask_epi8 (vTemp);
        }        
    }   

  short * pH = (short*) pvHStore;

  score = -0x7fff;       
  int offset = 0, counter = 0;
  short temp;
  for (i = 0; i < iter; ++i)     
    {    
      counter = 0;
      for (j = i; j < queryLength; j += iter)     
        {                
          temp = pH[offset + counter];
          if (score == temp) 
            {
              ++n_best;
              if (dir == 1 && j >= query_end) 
                query_end = j;
            }
          else if (score < temp) 
            {                        
              score = temp;
              n_best = 1;
              query_end = j;
            }                 
          counter++;
        }
      offset += 8;
    }

  (*_opt) = score;
  (*_te) = dbLength-1;
  (*_qe) = query_end;
  (*_n_best) = n_best;
}

string currentQuery="";

Solution10::Solution10() {
}
  
Solution10::~Solution10() {
}

// b: target (1-1024), a: query (1-512)
// qsc: query start clip (bool), qec: query end clip (bool)
// mm: match score (1,10), mi: mismatch score (-10 -1)
// o: gap openning (-10 -2), e: gap extension o+1 to -1
// dir: direction (bool)
int Solution10::process(string& b, string& a, int qsc, int qec,
                           int mm, int mi, int o, int e, int dir,
                           int *opt, int *te, int *qe, int *n_best) {
    int n = b.size(), m = a.size();

    int id = ((qsc << 1) | qec);         

    int i;

    SwStripedData *swData;    

    if (unlikely(flag)) {                
        gapOpen = -(o+e);
        gapExtend = -e;
        buildMatrix(mm,mi);
        flag=0;
    }

    for (i = 0; i < m; ++i)
      querySeq[i] = DNA_VALUE[a[i]-65];
    for (i = 0; i < n; ++i)
      targetSeq[i] = DNA_VALUE[b[i]-65];            
    if (id == 3) {    
        int buildProfile = 0;
        if (a!=currentQuery) {
            buildProfile = 1;
            currentQuery = a;
        }
        else
          buildProfile = 0;
        swStripedScan3 (querySeq, m, targetSeq, n, 
                        gapOpen, gapExtend,
                        dir, opt, te, qe, n_best, buildProfile);                  
    } else if (id == 2) {
        swData = swStripedInitWord(targetSeq, n, matrix);        
        swStripedWord2(n, querySeq, m ,
                       gapOpen, gapExtend,
                       swData->pvsQueryProf,
                       swData->pvH1,
                       swData->pvH2,
                       swData->pvE,
                       dir, opt, te, qe, n_best);                    
    } else if (id == 1) {
        swData = swStripedInitWord(targetSeq, n, matrix);    
        swStripedWord1 (n, querySeq, m,
                        gapOpen, gapExtend,
                        swData->pvsQueryProf,
                        swData->pvH1,
                        swData->pvH2,
                        swData->pvE,
                        dir, opt, te, qe, n_best);                    
    } else if (id == 0) {
        swData = swStripedInitWord(targetSeq, n, matrix);    
        swStripedWord0 (n, querySeq, m ,
                        gapOpen, gapExtend,
                        swData->pvsQueryProf,
                        swData->pvH1,
                        swData->pvH2,
                        swData->pvE,
                        dir, opt, te, qe, n_best);                
    }

    if (id != 3) {
        free (swData->pData);
        free (swData);
    }
    return (*opt);
}
