/*

   BWT-Index Construction

   This module constructs BWT and auxiliary data structures.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef FMAP_BWT_GEN_H
#define FMAP_BWT_GEN_H

/*! 
  BWT Generation Library
  */

/* START -- large BWT construction code */
#define ALPHABET_SIZE				4
#define BIT_PER_CHAR				2
#define CHAR_PER_WORD				16
#define CHAR_PER_BYTE				4

#define BITS_IN_WORD 32
#define BITS_IN_BYTE 8
#define BYTES_IN_WORD 4

#define ALL_ONE_MASK 0xFFFFFFFF
#define DNA_OCC_CNT_TABLE_SIZE_IN_WORD	65536

#define BITS_PER_OCC_VALUE			16
#define OCC_VALUE_PER_WORD			2
#define OCC_INTERVAL				256
#define OCC_INTERVAL_MAJOR			65536

#define TRUE    1
#define FALSE   0

#define BWTINC_INSERT_SORT_NUM_ITEM 7

#define average(value1, value2)					( ((value1) & (value2)) + ((value1) ^ (value2)) / 2 )
#define min(value1, value2)						( ((value1) < (value2)) ? (value1) : (value2) )
#define max(value1, value2)						( ((value1) > (value2)) ? (value1) : (value2) )
#define med3(a, b, c)							( a<b ? (b<c ? b : a<c ? c : a) : (b>c ? b : a>c ? c : a))
#define swap(a, b, t);							t = a; a = b; b = t;
#define truncateLeft(value, offset)				( (value) << (offset) >> (offset) )
#define truncateRight(value, offset)			( (value) >> (offset) << (offset) )
#define DNA_OCC_SUM_EXCEPTION(sum)			((sum & 0xfefefeff) == 0)

typedef struct SaIndexRange {
    uint32_t startSaIndex;
    uint32_t endSaIndex;
} SaIndexRange;

typedef struct BWT {
    uint32_t textLength;			// length of the text
    uint32_t saInterval;			// interval between two SA values stored explicitly
    uint32_t inverseSaInterval;		// interval between two inverse SA stored explicitly
    uint32_t inverseSa0;			// SA-1[0]
    uint32_t *cumulativeFreq;		// cumulative frequency
    uint32_t *bwtCode;				// BWT code
    uint32_t *occValue;				// Occurrence values stored explicitly
    uint32_t *occValueMajor;		// Occurrence values stored explicitly
    uint32_t *saValue;				// SA values stored explicitly
    uint32_t *inverseSa;			// Inverse SA stored explicitly
    SaIndexRange *saIndexRange;			// SA index range
    uint32_t saIndexRangeNumOfChar;			// Number of characters indexed in SA index range
    uint32_t *saValueOnBoundary;	// Pre-calculated frequently referred data
    uint32_t *decodeTable;			// For decoding BWT by table lookup
    uint32_t decodeTableGenerated;	// == TRUE if decode table is generated on load and will be freed
    uint32_t bwtSizeInWord;			// Temporary variable to hold the memory allocated
    uint32_t occSizeInWord;			// Temporary variable to hold the memory allocated
    uint32_t occMajorSizeInWord;	// Temporary variable to hold the memory allocated
    uint32_t saValueSize;			// Temporary variable to hold the memory allocated
    uint32_t inverseSaSize;			// Temporary variable to hold the memory allocated
    uint32_t saIndexRangeSize;		// Temporary variable to hold the memory allocated
} BWT;

typedef struct BWTInc {
    BWT *bwt;
    uint32_t numberOfIterationDone;
    uint32_t *cumulativeCountInCurrentBuild;
    uint32_t availableWord;
    uint32_t targetTextLength;
    float targetNBit;
    uint32_t buildSize;
    uint32_t initialMaxBuildSize;
    uint32_t incMaxBuildSize;
    uint32_t firstCharInLastIteration;
    uint32_t *workingMemory;
    uint32_t *packedText;
    uint8_t *textBuffer;
    uint32_t *packedShift;
} BWTInc;
/* END -- large BWT construction code */

/*! 
  create a bwt FASTA file from a packed FASTA file
  @param  fn_fasta      file name of the FASTA file
  @param  is_large      0 to use the short BWT construction algorith, 1 otherwise (large BWT construction algorithm) 
  @param  occ_interval  the desired occurrence interval
  @param  hash_width    the desired k-mer hash width
  */
void 
fmap_bwt_pac2bwt(const char *fn_fasta, uint32_t is_large, int32_t occ_interval, uint32_t hash_width);

/*! 
  @param  T  the input string
  @param  n  the length of the input string
  @return    the primary index if no error occurred, -1 or -2 otherwise
  */
uint32_t 
fmap_bwt_gen_short(uint8_t *T, uint32_t n);

#endif
