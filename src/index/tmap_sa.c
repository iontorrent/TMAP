/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_definitions.h"
#include "../io/tmap_file.h"
#include "tmap_bwt.h"
#include "tmap_bwt_gen.h"
#include "tmap_sa.h"

tmap_sa_t *
tmap_sa_read(const char *fn_fasta, uint32_t is_rev)
{
  char *fn_sa = NULL;
  tmap_file_t *fp_sa = NULL;
  tmap_sa_t *sa = NULL;

  fn_sa = tmap_get_file_name(fn_fasta, (0 == is_rev) ? TMAP_SA_FILE : TMAP_REV_SA_FILE);
  fp_sa = tmap_file_fopen(fn_sa, "rb", (0 == is_rev) ? TMAP_SA_COMPRESSION : TMAP_REV_SA_COMPRESSION);

  sa = tmap_calloc(1, sizeof(tmap_sa_t), "sa");

  if(1 != tmap_file_fread(&sa->primary, sizeof(tmap_bwt_int_t), 1, fp_sa)
     || 1 != tmap_file_fread(&sa->sa_intv, sizeof(tmap_bwt_int_t), 1, fp_sa)
     || 1 != tmap_file_fread(&sa->seq_len, sizeof(tmap_bwt_int_t), 1, fp_sa)
     || 1 != tmap_file_fread(&sa->is_rev, sizeof(uint32_t), 1, fp_sa)) {
      tmap_error(NULL, Exit, ReadFileError);
  }

  if(is_rev != sa->is_rev) {
      tmap_error("is_rev != sa->is_rev", Exit, OutOfRange);
  }

  sa->n_sa = (sa->seq_len + sa->sa_intv) / sa->sa_intv;
  sa->sa = tmap_calloc(sa->n_sa, sizeof(tmap_bwt_int_t), "sa->sa");
  sa->sa[0] = -1;

  if(sa->n_sa-1 != tmap_file_fread(sa->sa + 1, sizeof(tmap_bwt_int_t), sa->n_sa - 1, fp_sa)) {
      tmap_error(NULL, Exit, ReadFileError);
  }

  tmap_file_fclose(fp_sa);
  free(fn_sa);

  sa->is_shm = 0;

  return sa;
}

void 
tmap_sa_write(const char *fn_fasta, tmap_sa_t *sa, uint32_t is_rev)
{
  char *fn_sa = NULL;
  tmap_file_t *fp_sa = NULL;
  
  if(is_rev != sa->is_rev) {
      tmap_error("is_rev != sa->is_rev", Exit, OutOfRange);
  }

  fn_sa = tmap_get_file_name(fn_fasta, (0 == is_rev) ? TMAP_SA_FILE : TMAP_REV_SA_FILE);
  fp_sa = tmap_file_fopen(fn_sa, "wb", (0 == is_rev) ? TMAP_SA_COMPRESSION : TMAP_REV_SA_COMPRESSION);

  if(1 != tmap_file_fwrite(&sa->primary, sizeof(tmap_bwt_int_t), 1, fp_sa)
     || 1 != tmap_file_fwrite(&sa->sa_intv, sizeof(tmap_bwt_int_t), 1, fp_sa) 
     || 1 != tmap_file_fwrite(&sa->seq_len, sizeof(tmap_bwt_int_t), 1, fp_sa)
     || 1 != tmap_file_fwrite(&sa->is_rev, sizeof(uint32_t), 1, fp_sa)
     || sa->n_sa-1 != tmap_file_fwrite(sa->sa+1, sizeof(tmap_bwt_int_t), sa->n_sa-1, fp_sa)) {
      tmap_error(NULL, Exit, WriteFileError);
  }

  tmap_file_fclose(fp_sa);
  free(fn_sa);
}

size_t
tmap_sa_shm_num_bytes(tmap_sa_t *sa)
{
  // returns the number of bytes to allocate for shared memory
  size_t n = 0;
  
  n += sizeof(tmap_bwt_int_t); // primary
  n += sizeof(tmap_bwt_int_t); // sa_intv
  n += sizeof(tmap_bwt_int_t); // seq_len
  n += sizeof(uint32_t); // is_rev
  n += sizeof(tmap_bwt_int_t); // n_sa
  n += sizeof(tmap_bwt_int_t)*sa->n_sa; // sa

  return n;
}

size_t
tmap_sa_shm_read_num_bytes(const char *fn_fasta, uint32_t is_rev)
{
  size_t n = 0;
  char *fn_sa = NULL;
  tmap_file_t *fp_sa = NULL;
  tmap_sa_t *sa = NULL;

  fn_sa = tmap_get_file_name(fn_fasta, (0 == is_rev) ? TMAP_SA_FILE : TMAP_REV_SA_FILE);
  fp_sa = tmap_file_fopen(fn_sa, "rb", (0 == is_rev) ? TMAP_SA_COMPRESSION : TMAP_REV_SA_COMPRESSION);

  sa = tmap_calloc(1, sizeof(tmap_sa_t), "sa");

  if(1 != tmap_file_fread(&sa->primary, sizeof(tmap_bwt_int_t), 1, fp_sa)
     || 1 != tmap_file_fread(&sa->sa_intv, sizeof(tmap_bwt_int_t), 1, fp_sa)
     || 1 != tmap_file_fread(&sa->seq_len, sizeof(tmap_bwt_int_t), 1, fp_sa)
     || 1 != tmap_file_fread(&sa->is_rev, sizeof(uint32_t), 1, fp_sa)) {
      tmap_error(NULL, Exit, ReadFileError);
  }

  if(is_rev != sa->is_rev) {
      tmap_error("is_rev != sa->is_rev", Exit, OutOfRange);
  }

  sa->n_sa = (sa->seq_len + sa->sa_intv) / sa->sa_intv;

  // No need to read in sa->sa
  sa->sa = NULL;

  tmap_file_fclose(fp_sa);
  free(fn_sa);

  sa->is_shm = 0;

  // get the number of bytes
  n = tmap_sa_shm_num_bytes(sa);

  tmap_sa_destroy(sa);

  return n;
}

uint8_t *
tmap_sa_shm_pack(tmap_sa_t *sa, uint8_t *buf)
{
  // fixed length data
  memcpy(buf, &sa->primary, sizeof(tmap_bwt_int_t)); buf += sizeof(tmap_bwt_int_t);
  memcpy(buf, &sa->sa_intv, sizeof(tmap_bwt_int_t)); buf += sizeof(tmap_bwt_int_t);
  memcpy(buf, &sa->seq_len, sizeof(tmap_bwt_int_t)); buf += sizeof(tmap_bwt_int_t);
  memcpy(buf, &sa->is_rev, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(buf, &sa->n_sa, sizeof(tmap_bwt_int_t)); buf += sizeof(tmap_bwt_int_t);
  // variable length data
  memcpy(buf, sa->sa, sa->n_sa*sizeof(tmap_bwt_int_t)); buf += sa->n_sa*sizeof(tmap_bwt_int_t);

  return buf;
}

tmap_sa_t *
tmap_sa_shm_unpack(uint8_t *buf)
{
  tmap_sa_t *sa = NULL;

  if(NULL == buf) return NULL;

  sa = tmap_calloc(1, sizeof(tmap_sa_t), "sa");

  // fixed length data
  memcpy(&sa->primary, buf, sizeof(tmap_bwt_int_t)); buf += sizeof(tmap_bwt_int_t);
  memcpy(&sa->sa_intv, buf, sizeof(tmap_bwt_int_t)); buf += sizeof(tmap_bwt_int_t);
  memcpy(&sa->seq_len, buf, sizeof(tmap_bwt_int_t)); buf += sizeof(tmap_bwt_int_t);
  memcpy(&sa->is_rev, buf, sizeof(uint32_t)); buf += sizeof(uint32_t);
  memcpy(&sa->n_sa, buf, sizeof(tmap_bwt_int_t)); buf += sizeof(tmap_bwt_int_t);
  // variable length data
  sa->sa = (tmap_bwt_int_t*)buf;
  buf += sa->n_sa*sizeof(tmap_bwt_int_t);
  
  sa->is_shm = 1;

  return sa;
}

void 
tmap_sa_destroy(tmap_sa_t *sa)
{
  if(1 == sa->is_shm) {
  free(sa);
  }
  else {
  free(sa->sa);
  free(sa);
  }
}

tmap_bwt_int_t 
tmap_sa_pac_pos(const tmap_sa_t *sa, const tmap_bwt_t *bwt, tmap_bwt_int_t k)
{
  tmap_bwt_int_t s = 0;

  while (k % sa->sa_intv != 0) {
      ++s;
      k = tmap_bwt_invPsi(bwt, k);
  }
  /* without setting bwt->sa[0] = -1, the following line should be
     changed to (s + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1)
     */
  return s + sa->sa[k/sa->sa_intv];
}

extern int32_t debug_on;

void
tmap_sa_bwt2sa(const char *fn_fasta, uint32_t intv)
{
  int64_t isa, s, i; // S(isa) = sa
  tmap_bwt_t *bwt = NULL;
  tmap_sa_t *sa = NULL;
  uint32_t is_rev;

  for(is_rev=0;is_rev<=1;is_rev++) {
      if(is_rev == 0) {
          tmap_progress_print("constructing the SA from the BWT string");
      }
      else {
          tmap_progress_print("constructing the reverse SA from the reverse BWT string");
      }

      bwt = tmap_bwt_read(fn_fasta, is_rev);

      sa = tmap_calloc(1, sizeof(tmap_sa_t), "sa");

      sa->primary = bwt->primary;
      sa->sa_intv = intv;
      sa->seq_len = bwt->seq_len;
      sa->is_rev = is_rev;
      sa->n_sa = (bwt->seq_len + intv) / intv;
      
      // calculate SA value
      sa->sa = tmap_calloc(sa->n_sa, sizeof(tmap_bwt_int_t), "sa->sa");
      isa = 0; s = bwt->seq_len;
      for(i = 0; i < bwt->seq_len; ++i) {
          if(isa % intv == 0) sa->sa[isa/intv] = s;
          --s;
          isa = tmap_bwt_invPsi(bwt, isa);
      }
      if(isa % intv == 0) sa->sa[isa/intv] = s;
      sa->sa[0] = (tmap_bwt_int_t)-1; // before this line, bwt->sa[0] = bwt->seq_len

      tmap_sa_write(fn_fasta, sa, is_rev);

      tmap_bwt_destroy(bwt);
      tmap_sa_destroy(sa);
      sa=NULL;
      bwt=NULL;

      if(is_rev == 0) {
          tmap_progress_print2("constructed the SA from the BWT string");
      }
      else {
          tmap_progress_print2("constructed the reverse SA from the reverse BWT string");
      }
  }
}

/* 
   Original source from qsufsort.c

   Copyright 1999, N. Jesper Larsson, all rights reserved.

   This file contains an implementation of the algorithm presented in "Faster
   Suffix Sorting" by N. Jesper Larsson (jesper@cs.lth.se) and Kunihiko
   Sadakane (sada@is.s.u-tokyo.ac.jp).

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

   Modified by Wong Chi-Kwong, 2004

   Changes summary:	- Used long variable and function names
   - Removed global variables
   - Replace pointer references with array references
   - Used insertion sort in place of selection sort and increased insertion sort threshold
   - Reconstructing suffix array from inverse becomes an option
   - Add handling where end-of-text symbol is not necessary < all characters
   - Removed codes for supporting alphabet size > number of characters

   No warrenty is given regarding the quality of the modifications.

*/

#define BITS_IN_WORD 32

#define min(value1, value2)						( ((value1) < (value2)) ? (value1) : (value2) )
#define max(value1, value2)						( ((value1) > (value2)) ? (value1) : (value2) )
#define med3(a, b, c)							( a<b ? (b<c ? b : a<c ? c : a) : (b>c ? b : a>c ? c : a))
#define swap(a, b, t);							t = a; a = b; b = t;

// Static functions
static void QSufSortSortSplit(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t lowestPos, 
                              const tmap_bwt_sint_t highestPos, const tmap_bwt_sint_t numSortedChar);
static tmap_bwt_sint_t QSufSortChoosePivot(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t lowestPos, 
                                   const tmap_bwt_sint_t highestPos, const tmap_bwt_sint_t numSortedChar);
static void QSufSortInsertSortSplit(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t lowestPos, 
                                    const tmap_bwt_sint_t highestPos, const tmap_bwt_sint_t numSortedChar);
static void QSufSortBucketSort(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t numChar, const tmap_bwt_sint_t alphabetSize);
static tmap_bwt_sint_t QSufSortTransform(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t numChar, const tmap_bwt_sint_t largestInputSymbol, 
                                 const tmap_bwt_sint_t smallestInputSymbol, const tmap_bwt_sint_t maxNewAlphabetSize, tmap_bwt_sint_t *numSymbolAggregated);

// from MiscUtilities.c
static tmap_bwt_sint_t leadingZero(const tmap_bwt_sint_t input) {

    tmap_bwt_sint_t l;
    const static tmap_bwt_sint_t leadingZero8bit[256] = {8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
        2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    if (input & 0xFFFF0000) {
        if (input & 0xFF000000) {
            l = leadingZero8bit[input >> 24];
        } else {
            l = 8 + leadingZero8bit[input >> 16];
        }
    } else {
        if (input & 0x0000FF00) {
            l = 16 + leadingZero8bit[input >> 8];
        } else {
            l = 24 + leadingZero8bit[input];
        }
    }
    return l;

}

/* Makes suffix array p of x. x becomes inverse of p. p and x are both of size
   n+1. Contents of x[0...n-1] are integers in the range l...k-1. Original
   contents of x[n] is disregarded, the n-th symbol being regarded as
   end-of-string smaller than all other symbols.*/
void QSufSortSuffixSort(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t numChar, const tmap_bwt_sint_t largestInputSymbol, 
                        const tmap_bwt_sint_t smallestInputSymbol, const tmap_bwt_sint_t skipTransform) {

    tmap_bwt_sint_t i, j;
    tmap_bwt_sint_t s, negatedSortedGroupLength;
    tmap_bwt_sint_t numSymbolAggregated;
    tmap_bwt_sint_t maxNumInputSymbol;
    tmap_bwt_sint_t numSortedPos = 1;
    tmap_bwt_sint_t newAlphabetSize;

    maxNumInputSymbol = largestInputSymbol - smallestInputSymbol + 1;

    if (!skipTransform) {
        /* bucketing possible*/
        newAlphabetSize = QSufSortTransform(V, I, numChar, largestInputSymbol, smallestInputSymbol, 
                                            numChar, &numSymbolAggregated);
        QSufSortBucketSort(V, I, numChar, newAlphabetSize);
        I[0] = -1;
        V[numChar] = 0;
        numSortedPos = numSymbolAggregated;
    }

    while ((tmap_bwt_sint_t)(I[0]) >= -(tmap_bwt_sint_t)numChar) {
        i = 0;
        negatedSortedGroupLength = 0;
        do {
            s = I[i];
            if (s < 0) {
                i -= s;						/* skip over sorted group.*/
                negatedSortedGroupLength += s;
            } else {
                if (negatedSortedGroupLength) {
                    I[i+negatedSortedGroupLength] = negatedSortedGroupLength;	/* combine preceding sorted groups */
                    negatedSortedGroupLength = 0;
                }
                j = V[s] + 1;
                QSufSortSortSplit(V, I, i, j - 1, numSortedPos);
                i = j;
            }
        } while (i <= numChar);
        if (negatedSortedGroupLength) {
            /* array ends with a sorted group.*/
            I[i+negatedSortedGroupLength] = negatedSortedGroupLength;	/* combine sorted groups at end of I.*/
        }
        numSortedPos *= 2;	/* double sorted-depth.*/
    }

}
void QSufSortGenerateSaFromInverse(const tmap_bwt_sint_t* V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t numChar) {

    tmap_bwt_sint_t i;
    for (i=0; i<=numChar; i++) {
        I[V[i]] = i + 1;
    }

}

/* Sorting routine called for each unsorted group. Sorts the array of integers
   (suffix numbers) of length n starting at p. The algorithm is a ternary-split
   quicksort taken from Bentley & McIlroy, "Engineering a Sort Function",
   Software -- Practice and Experience 23(11), 1249-1265 (November 1993). This
   function is based on Program 7.*/
static void QSufSortSortSplit(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t lowestPos, 
                              const tmap_bwt_sint_t highestPos, const tmap_bwt_sint_t numSortedChar) {

    tmap_bwt_sint_t a, b, c, d;
    tmap_bwt_sint_t l, m;
    tmap_bwt_sint_t f, v, s, t;
    tmap_bwt_sint_t tmp;
    tmap_bwt_sint_t numItem;

#ifdef DEBUG
    if (lowestPos > highestPos) {
        tmap_error("QSufSortSortSplit(): lowestPos > highestPos!\n", Exit, OutOfRange)
    }
#endif

    numItem = highestPos - lowestPos + 1;

    if (numItem <= INSERT_SORT_NUM_ITEM) {
        QSufSortInsertSortSplit(V, I, lowestPos, highestPos, numSortedChar);
        return;
    }

    v = QSufSortChoosePivot(V, I, lowestPos, highestPos, numSortedChar);

    a = b = lowestPos;
    c = d = highestPos;

    while (1) {
        while (c >= b && (f = KEY(V, I, b, numSortedChar)) <= v) {
            if (f == v) {
                swap(I[a], I[b], tmp);
                a++;
            }
            b++;
        }
        while (c >= b && (f = KEY(V, I, c, numSortedChar)) >= v) {
            if (f == v) {
                swap(I[c], I[d], tmp);
                d--;
            }
            c--;
        }
        if (b > c) {
            break;
        }
        swap(I[b], I[c], tmp);
        b++;
        c--;
    }

    s = a - lowestPos;
    t = b - a;
    s = min(s, t);
    for (l = lowestPos, m = b - s; m < b; l++, m++) {
        swap(I[l], I[m], tmp);
    }

    s = d - c;
    t = highestPos - d;
    s = min(s, t);
    for (l = b, m = highestPos - s + 1; m <= highestPos; l++, m++) {
        swap(I[l], I[m], tmp);
    }

    s = b - a;
    t = d - c;
    if (s > 0) {
        QSufSortSortSplit(V, I, lowestPos, lowestPos + s - 1, numSortedChar);
    }

    // Update group number for equal portion
    a = lowestPos + s;
    b = highestPos - t;
    if (a == b) {
        // Sorted group
        V[I[a]] = a;
        I[a] = -1;
    } else {
        // Unsorted group
        for (c=a; c<=b; c++) {
            V[I[c]] = b;
        }
    }

    if (t > 0) {
        QSufSortSortSplit(V, I, highestPos - t + 1, highestPos, numSortedChar);
    }

}

/* Algorithm by Bentley & McIlroy.*/
static tmap_bwt_sint_t QSufSortChoosePivot(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t lowestPos, 
                                   const tmap_bwt_sint_t highestPos, const tmap_bwt_sint_t numSortedChar) {

    tmap_bwt_sint_t m;
    tmap_bwt_sint_t keyl, keym, keyn;
    tmap_bwt_sint_t key1, key2, key3;
    tmap_bwt_sint_t s;
    tmap_bwt_sint_t numItem;

#ifdef DEBUG
    if (lowestPos > highestPos) {
        tmap_error("QSufSortChoosePivot(): lowestPos > highestPos!\n", Exit, OutOfRange)
    }
#endif

    numItem = highestPos - lowestPos + 1;

#ifdef DEBUG
    if (numItem <= INSERT_SORT_NUM_ITEM) {
        tmap_error("QSufSortChoosePivot(): number of items <= INSERT_SORT_NUM_ITEM!\n", Exit, OutOfRange)
    }
#endif

    m = lowestPos + numItem / 2;

    s = numItem / 8;
    key1 = KEY(V, I, lowestPos, numSortedChar);
    key2 = KEY(V, I, lowestPos+s, numSortedChar);
    key3 = KEY(V, I, lowestPos+2*s, numSortedChar);
    keyl = med3(key1, key2, key3);
    key1 = KEY(V, I, m-s, numSortedChar);
    key2 = KEY(V, I, m, numSortedChar);
    key3 = KEY(V, I, m+s, numSortedChar);
    keym = med3(key1, key2, key3);
    key1 = KEY(V, I, highestPos-2*s, numSortedChar);
    key2 = KEY(V, I, highestPos-s, numSortedChar);
    key3 = KEY(V, I, highestPos, numSortedChar);
    keyn = med3(key1, key2, key3);

    return med3(keyl, keym, keyn);


}

/* Quadratic sorting method to use for small subarrays. */
static void QSufSortInsertSortSplit(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t lowestPos, 
                                    const tmap_bwt_sint_t highestPos, const tmap_bwt_sint_t numSortedChar) {

    tmap_bwt_sint_t i, j;
    tmap_bwt_sint_t tmpKey, tmpPos;
    tmap_bwt_sint_t numItem;
    tmap_bwt_sint_t key[INSERT_SORT_NUM_ITEM], pos[INSERT_SORT_NUM_ITEM];
    tmap_bwt_sint_t negativeSortedLength;
    tmap_bwt_sint_t groupNum;

#ifdef DEBUG
    if (lowestPos > highestPos) {
        tmap_error("QSufSortInsertSortSplit(): lowestPos > highestPos!\n", Exit, OutOfRange)
    }
#endif

    numItem = highestPos - lowestPos + 1;

#ifdef DEBUG
    if (numItem > INSERT_SORT_NUM_ITEM) {
        tmap_error("QSufSortInsertSortSplit(): number of items > INSERT_SORT_NUM_ITEM!\n", Exit, OutOfRange)
    }
#endif

    for (i=0; i<numItem; i++) {
#ifdef DEBUG
        if (I[lowestPos + i] < 0) {
            tmap_error("QSufSortInsertSortSplit(): I < 0 in unsorted region!\n", Exit, OutOfRange)
        }
#endif
        pos[i] = I[lowestPos + i];
        key[i] = V[pos[i] + numSortedChar];
    }

    for (i=1; i<numItem; i++) {
        tmpKey = key[i];
        tmpPos = pos[i];
        for (j=i; j>0 && key[j-1] > tmpKey; j--) {
            key[j] = key[j-1];
            pos[j] = pos[j-1];
        }
        key[j] = tmpKey;
        pos[j] = tmpPos;
    }

    negativeSortedLength = -1;

    i = numItem - 1;
    groupNum = highestPos;
    while (i > 0) {
        I[i+lowestPos] = pos[i];
        V[I[i+lowestPos]] = groupNum;
        if (key[i-1] == key[i]) {
            negativeSortedLength = 0;
        } else {
            if (negativeSortedLength < 0) {
                I[i+lowestPos] = negativeSortedLength;
            }
            groupNum = i + lowestPos - 1;
            negativeSortedLength--;
        }
        i--;
    }

    I[lowestPos] = pos[0];
    V[I[lowestPos]] = groupNum;
    if (negativeSortedLength < 0) {
        I[lowestPos] = negativeSortedLength;
    }	

}

/* Bucketsort for first iteration.

Input: x[0...n-1] holds integers in the range 1...k-1, all of which appear
at least once. x[n] is 0. (This is the corresponding output of transform.) k
must be at most n+1. p is array of size n+1 whose contents are disregarded.

Output: x is V and p is I after the initial sorting stage of the refined
suffix sorting algorithm.*/

static void QSufSortBucketSort(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t numChar, const tmap_bwt_sint_t alphabetSize) {

    tmap_bwt_sint_t i, c;
    tmap_bwt_sint_t d;
    tmap_bwt_sint_t groupNum;
    tmap_bwt_sint_t currentIndex;

    // mark linked list empty
    for (i=0; i<alphabetSize; i++) {
        I[i] = -1;
    }

    // insert to linked list
    for (i=0; i<=numChar; i++) {
        c = V[i];
        V[i] = (tmap_bwt_sint_t)(I[c]);
        I[c] = i;
    }

    currentIndex = numChar;
    for (i=alphabetSize; i>0; i--) {
        c = I[i-1];
        d = (tmap_bwt_sint_t)(V[c]);
        groupNum = currentIndex;
        V[c] = groupNum;
        if (d >= 0) {
            I[currentIndex] = c;
            while (d >= 0) {
                c = d;
                d = V[c];
                V[c] = groupNum;
                currentIndex--;
                I[currentIndex] = c;
            }
        } else {
            // sorted group
            I[currentIndex] = -1;
        }
        currentIndex--;
    }

}

/* Transforms the alphabet of x by attempting to aggregate several symbols into
   one, while preserving the suffix order of x. The alphabet may also be
   compacted, so that x on output comprises all integers of the new alphabet
   with no skipped numbers.

Input: x is an array of size n+1 whose first n elements are positive
integers in the range l...k-1. p is array of size n+1, used for temporary
storage. q controls aggregation and compaction by defining the maximum intue
for any symbol during transformation: q must be at least k-l; if q<=n,
compaction is guaranteed; if k-l>n, compaction is never done; if q is
INT_MAX, the maximum number of symbols are aggregated into one.

Output: Returns an integer j in the range 1...q representing the size of the
new alphabet. If j<=n+1, the alphabet is compacted. The global variable r is
set to the number of old symbols grouped into one. Only x[n] is 0.*/
static tmap_bwt_sint_t 
QSufSortTransform(tmap_bwt_sint_t* __restrict V, tmap_bwt_sint_t* __restrict I, const tmap_bwt_sint_t numChar, const tmap_bwt_sint_t largestInputSymbol, 
                                 const tmap_bwt_sint_t smallestInputSymbol, const tmap_bwt_sint_t maxNewAlphabetSize, tmap_bwt_sint_t *numSymbolAggregated) 
{

    tmap_bwt_sint_t c, i, j;
    tmap_bwt_sint_t a;	// numSymbolAggregated
    tmap_bwt_sint_t mask;
    tmap_bwt_sint_t minSymbolInChunk = 0, maxSymbolInChunk = 0;
    tmap_bwt_sint_t newAlphabetSize;
    tmap_bwt_sint_t maxNumInputSymbol, maxNumBit, maxSymbol;

    maxNumInputSymbol = largestInputSymbol - smallestInputSymbol + 1;

    maxNumBit = BITS_IN_WORD - leadingZero(maxNumInputSymbol);
    maxSymbol = INT_MAX >> maxNumBit;

    c = maxNumInputSymbol;
    for (a = 0; a < numChar && maxSymbolInChunk <= maxSymbol && c <= maxNewAlphabetSize; a++) {
        minSymbolInChunk = (minSymbolInChunk << maxNumBit) | (V[a] - smallestInputSymbol + 1);
        maxSymbolInChunk = c;
        c = (maxSymbolInChunk << maxNumBit) | maxNumInputSymbol;
    }

    mask = (1 << (a-1) * maxNumBit) - 1;	/* mask masks off top old symbol from chunk.*/
    V[numChar] = smallestInputSymbol - 1;	/* emulate zero terminator.*/

#ifdef DEBUG
    // Section of code for maxSymbolInChunk > numChar removed!
    if (maxSymbolInChunk > numChar) {
        tmap_error("QSufSortTransform(): maxSymbolInChunk > numChar!\n", Exit, OutOfRange)
    }
#endif

    /* bucketing possible, compact alphabet.*/
    for (i=0; i<=maxSymbolInChunk; i++) {
        I[i] = 0;	/* zero transformation table.*/
    }
    c = minSymbolInChunk;
    for (i=a; i<=numChar; i++) {
        I[c] = 1;			/* mark used chunk symbol.*/
        c = ((c & mask) << maxNumBit) | (V[i] - smallestInputSymbol + 1);	/* shift in next old symbol in chunk.*/
    }
    for (i=1; i<a; i++) {	/* handle last r-1 positions.*/
        I[c] = 1;			/* mark used chunk symbol.*/
        c = (c & mask) << maxNumBit;	/* shift in next old symbol in chunk.*/
    }
    newAlphabetSize = 1;
    for (i=0; i<=maxSymbolInChunk; i++) {
        if (I[i]) {
            I[i] = newAlphabetSize;
            newAlphabetSize++;
        }
    }
    c = minSymbolInChunk;
    for (i=0, j=a; j<=numChar; i++, j++) {
        V[i] = I[c];						/* transform to new alphabet.*/
        c = ((c & mask) << maxNumBit) | (V[j] - smallestInputSymbol + 1);	/* shift in next old symbol in chunk.*/
    }
    for (; i<numChar; i++) {	/* handle last a-1 positions.*/
        V[i] = I[c];			/* transform to new alphabet.*/
        c = (c & mask) << maxNumBit;	/* shift right-end zero in chunk.*/
    }

    V[numChar] = 0;		/* end-of-string symbol is zero.*/

    *numSymbolAggregated = a;
    return newAlphabetSize;

}

/*
 * sais.c for sais-lite
 * Copyright (c) 2008 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#define chr(i) (cs == sizeof(int32_t) ? ((const int32_t *)T)[i]:((const uint8_t*)T)[i])

/* find the start or end of each bucket */
static void 
getCounts(const unsigned char *T, int32_t *C, int32_t n, int32_t k, int32_t cs)
{
  int32_t i;
  for (i = 0; i < k; ++i) C[i] = 0;
  for (i = 0; i < n; ++i) ++C[chr(i)];
}

static void 
getBuckets(const int32_t *C, int32_t *B, int32_t k, int32_t end)
{
  int32_t i, sum = 0;
  if (end) {
      for (i = 0; i < k; ++i) {
          sum += C[i];
          B[i] = sum;
      }
  } else {
      for (i = 0; i < k; ++i) {
          sum += C[i];
          B[i] = sum - C[i];
      }
  }
}

/* compute SA */
static void 
induceSA(const unsigned char *T, int32_t *SA, int32_t *C, int32_t *B, int32_t n, int32_t k, int32_t cs)
{
  int32_t *b, i, j;
  int32_t  c0, c1;
  /* compute SAl */
  if (C == B) getCounts(T, C, n, k, cs);
  getBuckets(C, B, k, 0);	/* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = chr(j)];
  *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
  for (i = 0; i < n; ++i) {
      j = SA[i], SA[i] = ~j;
      if (0 < j) {
          --j;
          if ((c0 = chr(j)) != c1) {
              B[c1] = b - SA;
              b = SA + B[c1 = c0];
          }
          *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
      }
  }
  /* compute SAs */
  if (C == B) getCounts(T, C, n, k, cs);
  getBuckets(C, B, k, 1);	/* find ends of buckets */
  for (i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
      if (0 < (j = SA[i])) {
          --j;
          if ((c0 = chr(j)) != c1) {
              B[c1] = b - SA;
              b = SA + B[c1 = c0];
          }
          *--b = ((j == 0) || (chr(j - 1) > c1)) ? ~j : j;
      } else SA[i] = ~j;
  }
}

/*
 * find the suffix array SA of T[0..n-1] in {0..k-1}^n use a working
 * space (excluding T and SA) of at most 2n+O(1) for a constant alphabet
 */
static int32_t 
tmap_sa_sais_main(const unsigned char *T, int32_t *SA, int32_t fs, int32_t n, int32_t k, int32_t cs)
{
  int32_t *C, *B, *RA;
  int32_t  i, j, c, m, p, q, plen, qlen, name;
  int32_t  c0, c1;
  int32_t  diff;

  /* stage 1: reduce the problem by at least 1/2 sort all the
   * S-substrings */
  if (k <= fs) {
      C = SA + n;
      B = (k <= (fs - k)) ? C + k : C;
  } else if ((C = B = tmap_malloc(k * sizeof(int32_t), "C = B =")) == NULL) return -2;
  getCounts(T, C, n, k, cs);
  getBuckets(C, B, k, 1);	/* find ends of buckets */
  for (i = 0; i < n; ++i) SA[i] = 0;
  for (i = n - 2, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
      if ((c0 = chr(i)) < (c1 + c)) c = 1;
      else if (c != 0) SA[--B[c1]] = i + 1, c = 0;
  }
  induceSA(T, SA, C, B, n, k, cs);
  if (fs < k) free(C);
  /* compact all the sorted substrings into the first m items of SA
   * 2*m must be not larger than n (proveable) */
  for (i = 0, m = 0; i < n; ++i) {
      p = SA[i];
      if ((0 < p) && (chr(p - 1) > (c0 = chr(p)))) {
          for (j = p + 1; (j < n) && (c0 == (c1 = chr(j))); ++j);
          if ((j < n) && (c0 < c1)) SA[m++] = p;
      }
  }
  for (i = m; i < n; ++i) SA[i] = 0;	/* init the name array buffer */
  /* store the length of all substrings */
  for (i = n - 2, j = n, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
      if ((c0 = chr(i)) < (c1 + c)) c = 1;
      else if (c != 0) {
          SA[m + ((i + 1) >> 1)] = j - i - 1;
          j = i + 1;
          c = 0;
      }
  }
  /* find the lexicographic names of all substrings */
  for (i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
      p = SA[i], plen = SA[m + (p >> 1)], diff = 1;
      if (plen == qlen) {
          for (j = 0; (j < plen) && (chr(p + j) == chr(q + j)); j++);
          if (j == plen) diff = 0;
      }
      if (diff != 0) ++name, q = p, qlen = plen;
      SA[m + (p >> 1)] = name;
  }

  /* stage 2: solve the reduced problem recurse if names are not yet
   * unique */
  if (name < m) {
      RA = SA + n + fs - m;
      for (i = n - 1, j = m - 1; m <= i; --i) {
          if (SA[i] != 0) RA[j--] = SA[i] - 1;
      }
      if (tmap_sa_sais_main((unsigned char *) RA, SA, fs + n - m * 2, m, name, sizeof(int)) != 0) return -2;
      for (i = n - 2, j = m - 1, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
          if ((c0 = chr(i)) < (c1 + c)) c = 1;
          else if (c != 0) RA[j--] = i + 1, c = 0; /* get p1 */
      }
      for (i = 0; i < m; ++i) SA[i] = RA[SA[i]]; /* get index */
  }
  /* stage 3: induce the result for the original problem */
  if (k <= fs) {
      C = SA + n;
      B = (k <= (fs - k)) ? C + k : C;
  } else if ((C = B = (int32_t *) tmap_malloc(k * sizeof(int32_t), "C = B =")) == NULL) return -2;
  /* put all left-most S characters into their buckets */
  getCounts(T, C, n, k, cs);
  getBuckets(C, B, k, 1);	/* find ends of buckets */
  for (i = m; i < n; ++i) SA[i] = 0; /* init SA[m..n-1] */
  for (i = m - 1; 0 <= i; --i) {
      j = SA[i], SA[i] = 0;
      SA[--B[chr(j)]] = j;
  }
  induceSA(T, SA, C, B, n, k, cs);
  if (fs < k) free(C);
  return 0;
}

uint32_t 
tmap_sa_gen_short(const uint8_t *T, int32_t *SA, uint32_t n)
{
  if ((T == NULL) || (SA == NULL) || (n < 0)) return -1;
  SA[0] = n;
  if (n <= 1) {
      if (n == 1) SA[1] = 0;
      return 0;
  }
  return tmap_sa_sais_main(T, SA+1, 0, n, 256, 1);
}

int
tmap_sa_bwt2sa_main(int argc, char *argv[])
{
  int c, intv = TMAP_SA_INTERVAL, help=0;

  while((c = getopt(argc, argv, "i:vh")) >= 0) {
      switch(c) {
        case 'i': intv = atoi(optarg); break;
        case 'v': tmap_progress_set_verbosity(1); break;
        case 'h': help = 1; break;
        default: return 1;
      }
  }
  if(1 != argc - optind || 1 == help) {
      tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [-i INT -vh] <in.fasta>\n", PACKAGE, argv[0]);
      return 1;
  }
  if(intv <= 0 || (1 < intv && 0 != (intv % 2))) {
      tmap_error("option -i out of range", Exit, CommandLineArgument);
  }

  tmap_sa_bwt2sa(argv[optind], intv);

  return 0;
}
