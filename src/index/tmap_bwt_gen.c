/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/*

   BWTConstruct.c		BWT-Index Construction

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_definitions.h"
#include "../io/tmap_file.h"
#include "tmap_refseq.h"
#include "tmap_bwt.h"
#include "tmap_sa.h"
#include "tmap_bwt_gen.h"

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

#define MIN_AVAILABLE_WORD 0x10000

#define average(value1, value2)					( ((value1) & (value2)) + ((value1) ^ (value2)) / 2 )
#define min(value1, value2)						( ((value1) < (value2)) ? (value1) : (value2) )
#define truncateLeft(value, offset)				( (value) << (offset) >> (offset) )
#define truncateRight(value, offset)			( (value) >> (offset) << (offset) )
#define DNA_OCC_SUM_EXCEPTION(sum)			((sum & 0xfefefeff) == 0)

/*!
 Structure used to generate large BWT strings
 */
typedef struct {
    tmap_bwt_int_t textLength;  /*!<  length of the text */
    tmap_bwt_int_t inverseSa0;  /*!<  SA-1[0] */
    tmap_bwt_int_t *cumulativeFreq;  /*!<  cumulative frequency */
    uint32_t *bwtCode;  /*!<  BWT code */
    uint32_t *occValue;  /*!<  Occurrence values stored explicitly */
    tmap_bwt_int_t *occValueMajor;  /*!<  Occurrence values stored explicitly */
    uint32_t *decodeTable;  /*!<  For decoding BWT by table lookup */
    tmap_bwt_int_t bwtSizeInWord;  /*!<  Temporary variable to hold the memory allocated */
    tmap_bwt_int_t occSizeInWord;  /*!<  Temporary variable to hold the memory allocated */
    tmap_bwt_int_t occMajorSizeInWord; /*!< Temporary variable to hold the memory allocated */
} tmap_bwt_gen_t;

/*!
 Structure used to generate large BWT strings; for each iteration
 */
typedef struct {
    tmap_bwt_gen_t *bwt;
    uint32_t numberOfIterationDone;
    tmap_bwt_int_t *cumulativeCountInCurrentBuild;
    tmap_bwt_int_t availableWord;
    tmap_bwt_int_t buildSize;
    tmap_bwt_int_t initialMaxBuildSize;
    tmap_bwt_int_t incMaxBuildSize;
    uint32_t firstCharInLastIteration;
    uint32_t *workingMemory;
    uint32_t *packedText;
    uint8_t *textBuffer;
    uint32_t *packedShift;
} tmap_bwt_gen_inc_t;

static tmap_bwt_int_t
TextLengthFromBytePacked(tmap_bwt_int_t bytePackedLength, uint32_t bitPerChar,
                         uint32_t lastByteLength)
{
  /*
  if (bytePackedLength > ALL_ONE_MASK / (BITS_IN_BYTE / bitPerChar)) {
      tmap_error("TextLengthFromBytePacked(): text length > 2^32", Exit, OutOfRange);
  }
  */
  return (bytePackedLength - 1) * (BITS_IN_BYTE / bitPerChar) + lastByteLength;
}

static void 
initializeVAL(uint32_t *startAddr, const tmap_bwt_int_t length, const uint32_t initValue)
{
  tmap_bwt_int_t i;
  for (i=0; i<length; i++) startAddr[i] = initValue;
}

static void 
initializeVALBIG(tmap_bwt_int_t *startAddr, const tmap_bwt_int_t length, const tmap_bwt_int_t initValue)
{
  tmap_bwt_int_t i;
  for (i=0; i<length; i++) startAddr[i] = initValue;
}

static void 
GenerateDNAOccCountTable(uint32_t *dnaDecodeTable)
{
  uint32_t i, j, c, t;

  for (i=0; i<DNA_OCC_CNT_TABLE_SIZE_IN_WORD; i++) {
      dnaDecodeTable[i] = 0;
      c = i;
      for (j=0; j<8; j++) {
          t = c & 0x00000003;
          dnaDecodeTable[i] += 1 << (t * 8);
          c >>= 2;
      }
  }

}

// for BWTIncCreate()
static tmap_bwt_int_t 
BWTOccValueMajorSizeInWord(const tmap_bwt_int_t numChar)
{
  tmap_bwt_int_t numOfOccValue;
  unsigned numOfOccIntervalPerMajor;
  numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1; // Value at both end for bi-directional encoding
  numOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
  return (numOfOccValue + numOfOccIntervalPerMajor - 1) / numOfOccIntervalPerMajor * ALPHABET_SIZE;
}

// for BWTIncCreate()
static tmap_bwt_int_t 
BWTOccValueMinorSizeInWord(const tmap_bwt_int_t numChar)
{
  tmap_bwt_int_t numOfOccValue;
  numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;		// Value at both end for bi-directional encoding
  return (numOfOccValue + OCC_VALUE_PER_WORD - 1) / OCC_VALUE_PER_WORD * ALPHABET_SIZE;
}

// for BWTIncCreate()
static tmap_bwt_int_t 
BWTResidentSizeInWord(const tmap_bwt_int_t numChar) {

    tmap_bwt_int_t numCharRoundUpToOccInterval;

    // The $ in BWT at the position of inverseSa0 is not encoded
    numCharRoundUpToOccInterval = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL;

    return (numCharRoundUpToOccInterval + CHAR_PER_WORD - 1) / CHAR_PER_WORD;

}

static void 
BWTIncSetBuildSizeAndTextAddr(tmap_bwt_gen_inc_t *bwtInc)
{
  tmap_bwt_int_t maxBuildSize;

  if (bwtInc->bwt->textLength == 0) {
      // initial build
      // Minus 2 because n+1 entries of seq and rank needed for n char
      maxBuildSize = (bwtInc->availableWord - (2 + OCC_INTERVAL / CHAR_PER_WORD) * (sizeof(tmap_bwt_int_t) / 4))
        / (2 * CHAR_PER_WORD + 1) * CHAR_PER_WORD / (sizeof(tmap_bwt_int_t) / 4);
      if (bwtInc->initialMaxBuildSize > 0) {
          bwtInc->buildSize = min(bwtInc->initialMaxBuildSize, maxBuildSize);
      } else {
          bwtInc->buildSize = maxBuildSize;
      }
  } else {
      // Minus 3 because n+1 entries of sorted rank, seq and rank needed for n char
      // Minus numberOfIterationDone because bwt slightly shift to left in each iteration
      maxBuildSize = (bwtInc->availableWord - bwtInc->bwt->bwtSizeInWord - bwtInc->bwt->occSizeInWord
                      - (3 + bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR) * (sizeof(tmap_bwt_int_t) / 4))
        / 3 / (sizeof(tmap_bwt_int_t) / 4);

      if (maxBuildSize < CHAR_PER_WORD) {
          tmap_error("BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction", Exit, OutOfRange);
      }
      if (bwtInc->incMaxBuildSize > 0) {
          bwtInc->buildSize = min(bwtInc->incMaxBuildSize, maxBuildSize);
      } else {
          bwtInc->buildSize = maxBuildSize;
      }
      if (bwtInc->buildSize < CHAR_PER_WORD) {
          bwtInc->buildSize = CHAR_PER_WORD;
      }
  }

  if (bwtInc->buildSize < CHAR_PER_WORD) {
      tmap_error("BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction", Exit, OutOfRange);
  }

  bwtInc->buildSize = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;

  bwtInc->packedText = bwtInc->workingMemory + 2 * (bwtInc->buildSize + 1) * (sizeof(tmap_bwt_int_t) / 4);
  bwtInc->textBuffer = (uint8_t*)(bwtInc->workingMemory + (bwtInc->buildSize + 1) * (sizeof(tmap_bwt_int_t) / 4));
}

// for ceilLog2()
uint32_t 
leadingZero(const uint32_t input)
{
  uint32_t l;
  static const uint32_t leadingZero8bit[256] = {8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
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

// for BitPerBytePackedChar()
static uint32_t 
ceilLog2(const uint32_t input)
{
  if (input <= 1) return 0;
  return BITS_IN_WORD - leadingZero(input - 1);

}

// for ConvertBytePackedToWordPacked()
static uint32_t 
BitPerBytePackedChar(const uint32_t alphabetSize)
{
  uint32_t bitPerChar;
  bitPerChar = ceilLog2(alphabetSize);
  // Return the largest number of bit that does not affect packing efficiency
  if (BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar) > bitPerChar)
    bitPerChar = BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar);
  return bitPerChar;
}

// for ConvertBytePackedToWordPacked()
static uint32_t 
BitPerWordPackedChar(const uint32_t alphabetSize)
{
  return ceilLog2(alphabetSize);
}

static void 
ConvertBytePackedToWordPacked(const  uint8_t *input, uint32_t *output, const uint32_t alphabetSize,
                              const tmap_bwt_int_t textLength)
{
  tmap_bwt_int_t i;
  uint32_t j, k, c;
  uint32_t bitPerBytePackedChar;
  uint32_t bitPerWordPackedChar;
  uint32_t charPerWord;
  uint32_t charPerByte;
  uint32_t bytePerIteration;
  tmap_bwt_int_t byteProcessed = 0;
  tmap_bwt_int_t wordProcessed = 0;
  uint32_t mask, shift;

  uint32_t buffer[BITS_IN_WORD];

  bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
  bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
  charPerByte = BITS_IN_BYTE / bitPerBytePackedChar;
  charPerWord = BITS_IN_WORD / bitPerWordPackedChar;

  bytePerIteration = charPerWord / charPerByte;
  mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
  shift = BITS_IN_WORD - BITS_IN_BYTE + bitPerBytePackedChar - bitPerWordPackedChar;

  while ((wordProcessed + 1) * charPerWord < textLength) {

      k = 0;
      for (i=0; i<bytePerIteration; i++) {
          c = (uint32_t)input[byteProcessed] << shift;
          for (j=0; j<charPerByte; j++) {
              buffer[k] = c & mask;
              c <<= bitPerBytePackedChar;
              k++;
          }
          byteProcessed++;
      }

      c = 0;
      for (i=0; i<charPerWord; i++) {
          c |= buffer[i] >> bitPerWordPackedChar * i;
      }
      output[wordProcessed] = c;
      wordProcessed++;

  }

  k = 0;
  for (i=0; i < (textLength - wordProcessed * charPerWord - 1) / charPerByte + 1; i++) {
      c = (uint32_t)input[byteProcessed] << shift;
      for (j=0; j<charPerByte; j++) {
          buffer[k] = c & mask;
          c <<= bitPerBytePackedChar;
          k++;
      }
      byteProcessed++;
  }

  c = 0;
  for (i=0; i<textLength - wordProcessed * charPerWord; i++) {
      c |= buffer[i] >> bitPerWordPackedChar * i;
  }
  output[wordProcessed] = c;
}

tmap_bwt_gen_t *
BWTCreate(const tmap_bwt_int_t textLength, uint32_t *decodeTable)
{
  tmap_bwt_gen_t *bwt;

  bwt = tmap_calloc(1, sizeof(tmap_bwt_gen_t), "bwt");

  bwt->textLength = 0;

  bwt->cumulativeFreq = tmap_calloc((ALPHABET_SIZE + 1), sizeof(tmap_bwt_int_t), "bwt->cumulativeFreq");
  initializeVALBIG(bwt->cumulativeFreq, ALPHABET_SIZE + 1, 0);

  bwt->bwtSizeInWord = 0;

  // Generate decode tables
  if (decodeTable == NULL) {
      bwt->decodeTable = tmap_calloc(DNA_OCC_CNT_TABLE_SIZE_IN_WORD, sizeof(uint32_t), "bwt->decodeTable");
      GenerateDNAOccCountTable(bwt->decodeTable);
  } else {
      bwt->decodeTable = decodeTable;
  }

  bwt->occValueMajor = tmap_calloc(BWTOccValueMajorSizeInWord(textLength), sizeof(tmap_bwt_int_t), "bwt->occValueMajor");

  bwt->occSizeInWord = 0;
  bwt->occValue = NULL;


  return bwt;
}

tmap_bwt_gen_inc_t *
BWTIncCreate(const tmap_bwt_int_t textLength, 
             const uint32_t initialMaxBuildSize, const uint32_t incMaxBuildSize)
{
  tmap_bwt_gen_inc_t *bwtInc;
  uint32_t i, n_iter;

  bwtInc = tmap_calloc(1, sizeof(tmap_bwt_gen_inc_t), "bwtInc");
  bwtInc->numberOfIterationDone = 0;
  bwtInc->bwt = BWTCreate(textLength, NULL);
  bwtInc->initialMaxBuildSize = initialMaxBuildSize;
  bwtInc->incMaxBuildSize = incMaxBuildSize;
  bwtInc->cumulativeCountInCurrentBuild = tmap_calloc((ALPHABET_SIZE + 1), sizeof(tmap_bwt_int_t), "bwtInc->cumumlativeCountInCurrentBuild");
  initializeVALBIG(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

  // Build frequently accessed data
  bwtInc->packedShift = tmap_calloc(CHAR_PER_WORD, sizeof(uint32_t), "bwtInc->packedShift");
  for (i=0; i<CHAR_PER_WORD; i++) {
      bwtInc->packedShift[i] = BITS_IN_WORD - (i+1) * BIT_PER_CHAR;
  }

  n_iter = (textLength - initialMaxBuildSize) / incMaxBuildSize + 1;
  bwtInc->availableWord = BWTResidentSizeInWord(textLength) + BWTOccValueMinorSizeInWord(textLength) // minimal memory requirement
    + OCC_INTERVAL / BIT_PER_CHAR * n_iter * 2 * (sizeof(tmap_bwt_int_t) / 4) // buffer at the end of occ array 
    + incMaxBuildSize/5 * 3 * (sizeof(tmap_bwt_int_t) / 4); // space for the 3 temporary arrays in each iteration
  if (bwtInc->availableWord < MIN_AVAILABLE_WORD) bwtInc->availableWord = MIN_AVAILABLE_WORD; // lh3: otherwise segfault when availableWord is too small
  //fprintf(stderr, "[%s] textLength=%ld, availableWord=%ld\n", __func__, (long)textLength, (long)bwtInc->availableWord);
  bwtInc->workingMemory = (unsigned*)calloc(bwtInc->availableWord, BYTES_IN_WORD);

  return bwtInc;
}

// for BWTIncConstruct()
static void 
BWTIncPutPackedTextToRank(const uint32_t *packedText, tmap_bwt_int_t* __restrict rank,
                          tmap_bwt_int_t* __restrict cumulativeCount, const tmap_bwt_int_t numChar)
{
  tmap_bwt_int_t i;
  uint32_t j;
  uint32_t c, t;
  uint32_t packedMask;
  tmap_bwt_int_t rankIndex, lastWord;
  uint32_t numCharInLastWord;

  lastWord = (numChar - 1) / CHAR_PER_WORD;
  numCharInLastWord = numChar - lastWord * CHAR_PER_WORD;

  packedMask = ALL_ONE_MASK >> (BITS_IN_WORD - BIT_PER_CHAR);
  rankIndex = numChar - 1;

  t = packedText[lastWord] >> (BITS_IN_WORD - numCharInLastWord * BIT_PER_CHAR);
  for (i=0; i<numCharInLastWord; i++) {
      c = t & packedMask;
      cumulativeCount[c+1]++;
      rank[rankIndex] = c;
      rankIndex--;
      t >>= BIT_PER_CHAR;
  }

  for (i=lastWord; i--;) {	// loop from lastWord - 1 to 0
      t = packedText[i];
      for (j=0; j<CHAR_PER_WORD; j++) {
          c = t & packedMask;
          cumulativeCount[c+1]++;
          rank[rankIndex] = c;
          rankIndex--;
          t >>= BIT_PER_CHAR;
      }
  }

  // Convert occurrence to cumulativeCount
  cumulativeCount[2] += cumulativeCount[1];
  cumulativeCount[3] += cumulativeCount[2];
  cumulativeCount[4] += cumulativeCount[3];
}

static void 
ForwardDNAAllOccCountNoLimit(const uint32_t*  dna, const tmap_bwt_int_t index,
                             tmap_bwt_int_t* __restrict occCount, const uint32_t*  dnaDecodeTable)
{
  static const uint32_t truncateRightMask[16] = { 0x00000000, 0xC0000000, 0xF0000000, 0xFC000000,
      0xFF000000, 0xFFC00000, 0xFFF00000, 0xFFFC0000,
      0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00,
      0xFFFFFF00, 0xFFFFFFC0, 0xFFFFFFF0, 0xFFFFFFFC };

  tmap_bwt_int_t iteration, i;
  uint32_t wordToCount, charToCount;
  uint32_t j, c;
  uint32_t sum;

  occCount[0] = 0;
  occCount[1] = 0;
  occCount[2] = 0;
  occCount[3] = 0;

  iteration = index / 256;
  wordToCount = (index - iteration * 256) / 16;
  charToCount = index - iteration * 256 - wordToCount * 16;

  for (i=0; i<iteration; i++) {

      sum = 0;
      for (j=0; j<16; j++) {
          sum += dnaDecodeTable[*dna >> 16];
          sum += dnaDecodeTable[*dna & 0x0000FFFF];
          dna++;
      }
      if (!DNA_OCC_SUM_EXCEPTION(sum)) {
          occCount[0] += sum & 0x000000FF;	sum >>= 8;
          occCount[1] += sum & 0x000000FF;	sum >>= 8;
          occCount[2] += sum & 0x000000FF;	sum >>= 8;
          occCount[3] += sum;
      } else {
          // only some or all of the 3 bits are on
          // in reality, only one of the four cases are possible
          if (sum == 0x00000100) {
              occCount[0] += 256;
          } else if (sum == 0x00010000) {
              occCount[1] += 256;
          } else if (sum == 0x01000000) {
              occCount[2] += 256;
          } else if (sum == 0x00000000) {
              occCount[3] += 256;
          } else {
              tmap_error("ForwardDNAAllOccCountNoLimit(): DNA occ sum exception", Exit, OutOfRange);
          }
      }

  }

  sum = 0;
  for (j=0; j<wordToCount; j++) {
      sum += dnaDecodeTable[*dna >> 16];
      sum += dnaDecodeTable[*dna & 0x0000FFFF];
      dna++;
  }

  if (charToCount > 0) {
      c = *dna & truncateRightMask[charToCount];	// increase count of 'a' by 16 - c;
      sum += dnaDecodeTable[c >> 16];
      sum += dnaDecodeTable[c & 0xFFFF];
      sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
  }

  occCount[0] += sum & 0x000000FF;	sum >>= 8;
  occCount[1] += sum & 0x000000FF;	sum >>= 8;
  occCount[2] += sum & 0x000000FF;	sum >>= 8;
  occCount[3] += sum;
}

static void 
BWTIncBuildPackedBwt(const tmap_bwt_int_t *relativeRank, uint32_t* __restrict bwt, const tmap_bwt_int_t numChar,
                     const tmap_bwt_int_t *cumulativeCount, const uint32_t *packedShift) {

    tmap_bwt_int_t i, r;
    uint32_t c;
    tmap_bwt_int_t previousRank, currentRank;
    tmap_bwt_int_t wordIndex, charIndex;
    tmap_bwt_int_t inverseSa0;

    inverseSa0 = previousRank = relativeRank[0];

    for (i=1; i<=numChar; i++) {
        currentRank = relativeRank[i];
        // previousRank > cumulativeCount[c] because $ is one of the char
        c = (previousRank > cumulativeCount[1]) + (previousRank > cumulativeCount[2]) 
          + (previousRank > cumulativeCount[3]);
        // set bwt for currentRank
        if (c > 0) {
            // c <> 'a'
            r = currentRank;
            if (r > inverseSa0) {
                // - 1 because $ at inverseSa0 is not encoded			
                r--;
            }
            wordIndex = r / CHAR_PER_WORD;
            charIndex = r - wordIndex * CHAR_PER_WORD;
            bwt[wordIndex] |= c << packedShift[charIndex];
        }
        previousRank = currentRank;
    }
}

static inline tmap_bwt_int_t
BWTOccValueExplicit(const tmap_bwt_gen_t *bwt, const tmap_bwt_int_t occIndexExplicit,
                    const uint32_t character)
{
  tmap_bwt_int_t occIndexMajor;

  occIndexMajor = occIndexExplicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR;

  if (occIndexExplicit % OCC_VALUE_PER_WORD == 0) {
      return bwt->occValueMajor[occIndexMajor * ALPHABET_SIZE + character] +
        (bwt->occValue[occIndexExplicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE + character] >> 16);

  } else {
      return bwt->occValueMajor[occIndexMajor * ALPHABET_SIZE + character] +
        (bwt->occValue[occIndexExplicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE + character] & 0x0000FFFF);
  }
}

static uint32_t 
ForwardDNAOccCount(const uint32_t*  dna, const uint32_t index, const uint32_t character,
                   const uint32_t*  dnaDecodeTable)
{
  static const uint32_t truncateRightMask[16] = { 0x00000000, 0xC0000000, 0xF0000000, 0xFC000000,
      0xFF000000, 0xFFC00000, 0xFFF00000, 0xFFFC0000,
      0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00,
      0xFFFFFF00, 0xFFFFFFC0, 0xFFFFFFF0, 0xFFFFFFFC };

  uint32_t wordToCount, charToCount;
  uint32_t i, c;
  uint32_t sum = 0;

  wordToCount = index / 16;
  charToCount = index - wordToCount * 16;

  for (i=0; i<wordToCount; i++) {
      sum += dnaDecodeTable[dna[i] >> 16];
      sum += dnaDecodeTable[dna[i] & 0x0000FFFF];
  }

  if (charToCount > 0) {
      c = dna[i] & truncateRightMask[charToCount];	// increase count of 'a' by 16 - c;
      sum += dnaDecodeTable[c >> 16];
      sum += dnaDecodeTable[c & 0xFFFF];
      sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
  }

  return (sum >> (character * 8)) & 0x000000FF;

}

static uint32_t 
BackwardDNAOccCount(const uint32_t*  dna, const uint32_t index, const uint32_t character,
                    const uint32_t*  dnaDecodeTable)
{
  static const uint32_t truncateLeftMask[16] =  { 0x00000000, 0x00000003, 0x0000000F, 0x0000003F,
      0x000000FF, 0x000003FF, 0x00000FFF, 0x00003FFF,
      0x0000FFFF, 0x0003FFFF, 0x000FFFFF, 0x003FFFFF,
      0x00FFFFFF, 0x03FFFFFF, 0x0FFFFFFF, 0x3FFFFFFF };

  uint32_t wordToCount, charToCount;
  uint32_t i, c;
  uint32_t sum = 0;

  wordToCount = index / 16;
  charToCount = index - wordToCount * 16;

  dna -= wordToCount + 1;

  if (charToCount > 0) {
      c = *dna & truncateLeftMask[charToCount];	// increase count of 'a' by 16 - c;
      sum += dnaDecodeTable[c >> 16];
      sum += dnaDecodeTable[c & 0xFFFF];
      sum += charToCount - 16;	// decrease count of 'a' by 16 - positionToProcess
  }

  for (i=0; i<wordToCount; i++) {
      dna++;
      sum += dnaDecodeTable[*dna >> 16];
      sum += dnaDecodeTable[*dna & 0x0000FFFF];
  }

  return (sum >> (character * 8)) & 0x000000FF;

}

tmap_bwt_int_t
BWTOccValue(const tmap_bwt_gen_t *bwt, tmap_bwt_int_t index, const uint32_t character) {

    tmap_bwt_int_t occValue;
    tmap_bwt_int_t occExplicitIndex, occIndex;

    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is subtracted by 1 for adjustment
    if (index > bwt->inverseSa0) {
        index--;
    }

    occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;	// Bidirectional encoding
    occIndex = occExplicitIndex * OCC_INTERVAL;
    occValue = BWTOccValueExplicit(bwt, occExplicitIndex, character);

    if (occIndex == index) {
        return occValue;
    }

    if (occIndex < index) {
        return occValue + ForwardDNAOccCount(bwt->bwtCode + occIndex / CHAR_PER_WORD, index - occIndex, character, bwt->decodeTable);
    } else {
        return occValue - BackwardDNAOccCount(bwt->bwtCode + occIndex / CHAR_PER_WORD, occIndex - index, character, bwt->decodeTable);
    }

}

static tmap_bwt_int_t 
BWTIncGetAbsoluteRank(tmap_bwt_gen_t *bwt, tmap_bwt_int_t* __restrict absoluteRank, tmap_bwt_int_t* __restrict seq,
                      const uint32_t *packedText, const tmap_bwt_int_t numChar,
                      const tmap_bwt_int_t* cumulativeCount, const uint32_t firstCharInLastIteration)
{
  tmap_bwt_int_t saIndex;
  tmap_bwt_int_t lastWord;
  uint32_t packedMask;
  tmap_bwt_int_t i;
  uint32_t c, t, j;
  tmap_bwt_int_t rankIndex;
  uint32_t shift;
  tmap_bwt_int_t seqIndexFromStart[ALPHABET_SIZE];
  tmap_bwt_int_t seqIndexFromEnd[ALPHABET_SIZE];

  for (i=0; i<ALPHABET_SIZE; i++) {
      seqIndexFromStart[i] = cumulativeCount[i];
      seqIndexFromEnd[i] = cumulativeCount[i+1] - 1;
  }

  shift = BITS_IN_WORD - BIT_PER_CHAR;
  packedMask = ALL_ONE_MASK >> shift;
  saIndex = bwt->inverseSa0;
  rankIndex = numChar - 1;

  lastWord = numChar / CHAR_PER_WORD;
  for (i=lastWord; i--;) {	// loop from lastWord - 1 to 0
      t = packedText[i];
      for (j=0; j<CHAR_PER_WORD; j++) {
          c = t & packedMask;
          saIndex = bwt->cumulativeFreq[c] + BWTOccValue(bwt, saIndex, c) + 1;
          // A counting sort using the first character of suffix is done here
          // If rank > inverseSa0 -> fill seq from end, otherwise fill seq from start -> to leave the right entry for inverseSa0
          if (saIndex > bwt->inverseSa0) {
              seq[seqIndexFromEnd[c]] = rankIndex;
              absoluteRank[seqIndexFromEnd[c]] = saIndex;
              seqIndexFromEnd[c]--;
          } else {
              seq[seqIndexFromStart[c]] = rankIndex;
              absoluteRank[seqIndexFromStart[c]] = saIndex;
              seqIndexFromStart[c]++;
          }
          rankIndex--;
          t >>= BIT_PER_CHAR;
      }
  }

  absoluteRank[seqIndexFromStart[firstCharInLastIteration]] = bwt->inverseSa0;	// representing the substring of all preceding characters
  seq[seqIndexFromStart[firstCharInLastIteration]] = numChar;

  return seqIndexFromStart[firstCharInLastIteration];
}

static void 
BWTIncSortKey(tmap_bwt_int_t* __restrict key, tmap_bwt_int_t* __restrict seq, const tmap_bwt_int_t numItem)
{
#define EQUAL_KEY_THRESHOLD	4	// Partition for equal key if data array size / the number of data with equal value with pivot < EQUAL_KEY_THRESHOLD

  int64_t lowIndex, highIndex, midIndex;
  int64_t lowPartitionIndex, highPartitionIndex;
  int64_t lowStack[32], highStack[32];
  uint32_t stackDepth;
  int64_t i, j;
  tmap_bwt_int_t tempSeq, tempKey;
  int64_t numberOfEqualKey;

  if (numItem < 2) return;

  stackDepth = 0;

  lowIndex = 0;
  highIndex = numItem - 1;

  for (;;) {

      for (;;) {

          // Sort small array of data
          if (highIndex - lowIndex < BWTINC_INSERT_SORT_NUM_ITEM) {	 // Insertion sort on smallest arrays
              for (i=lowIndex+1; i<=highIndex; i++) {
                  tempSeq = seq[i];
                  tempKey = key[i];
                  for (j = i; j > lowIndex && key[j-1] > tempKey; j--) {
                      seq[j] = seq[j-1];
                      key[j] = key[j-1];
                  }
                  if (j != i) {
                      seq[j] = tempSeq;
                      key[j] = tempKey;
                  }
              }
              break;
          }

          // Choose pivot as median of the lowest, middle, and highest data; sort the three data

          midIndex = average(lowIndex, highIndex);
          if (key[lowIndex] > key[midIndex]) {
              tempSeq = seq[lowIndex];
              tempKey = key[lowIndex];
              seq[lowIndex] = seq[midIndex];
              key[lowIndex] = key[midIndex];
              seq[midIndex] = tempSeq;
              key[midIndex] = tempKey;
          }
          if (key[lowIndex] > key[highIndex]) {
              tempSeq = seq[lowIndex];
              tempKey = key[lowIndex];
              seq[lowIndex] = seq[highIndex];
              key[lowIndex] = key[highIndex];
              seq[highIndex] = tempSeq;
              key[highIndex] = tempKey;
          }
          if (key[midIndex] > key[highIndex]) {
              tempSeq = seq[midIndex];
              tempKey = key[midIndex];
              seq[midIndex] = seq[highIndex];
              key[midIndex] = key[highIndex];
              seq[highIndex] = tempSeq;
              key[highIndex] = tempKey;
          }

          // Partition data

          numberOfEqualKey = 0;

          lowPartitionIndex = lowIndex + 1;
          highPartitionIndex = highIndex - 1;

          for (;;) {
              while (lowPartitionIndex <= highPartitionIndex && key[lowPartitionIndex] <= key[midIndex]) {
                  numberOfEqualKey += (key[lowPartitionIndex] == key[midIndex]);
                  lowPartitionIndex++;
              }
              while (lowPartitionIndex < highPartitionIndex) {
                  if (key[midIndex] >= key[highPartitionIndex]) {
                      numberOfEqualKey += (key[midIndex] == key[highPartitionIndex]);
                      break;
                  }
                  highPartitionIndex--;
              }
              if (lowPartitionIndex >= highPartitionIndex) {
                  break;
              }
              tempSeq = seq[lowPartitionIndex];
              tempKey = key[lowPartitionIndex];
              seq[lowPartitionIndex] = seq[highPartitionIndex];
              key[lowPartitionIndex] = key[highPartitionIndex];
              seq[highPartitionIndex] = tempSeq;
              key[highPartitionIndex] = tempKey;
              if (highPartitionIndex == midIndex) {
                  // partition key has been moved
                  midIndex = lowPartitionIndex;
              }
              lowPartitionIndex++;
              highPartitionIndex--;
          }

          // Adjust the partition index
          highPartitionIndex = lowPartitionIndex;
          lowPartitionIndex--;

          // move the partition key to end of low partition
          tempSeq = seq[midIndex];
          tempKey = key[midIndex];
          seq[midIndex] = seq[lowPartitionIndex];
          key[midIndex] = key[lowPartitionIndex];
          seq[lowPartitionIndex] = tempSeq;
          key[lowPartitionIndex] = tempKey;

          if (highIndex - lowIndex + BWTINC_INSERT_SORT_NUM_ITEM <= EQUAL_KEY_THRESHOLD * numberOfEqualKey) {

              // Many keys = partition key; separate the equal key data from the lower partition

              midIndex = lowIndex;

              for (;;) {
                  while (midIndex < lowPartitionIndex && key[midIndex] < key[lowPartitionIndex]) {
                      midIndex++;
                  }
                  while (midIndex < lowPartitionIndex && key[lowPartitionIndex] == key[lowPartitionIndex - 1]) {
                      lowPartitionIndex--;
                  }
                  if (midIndex >= lowPartitionIndex) {
                      break;
                  }
                  tempSeq = seq[midIndex];
                  tempKey = key[midIndex];
                  seq[midIndex] = seq[lowPartitionIndex - 1];
                  key[midIndex] = key[lowPartitionIndex - 1];
                  seq[lowPartitionIndex - 1] = tempSeq;
                  key[lowPartitionIndex - 1] = tempKey;
                  midIndex++;
                  lowPartitionIndex--;
              }

          }

          if (lowPartitionIndex - lowIndex > highIndex - highPartitionIndex) {
              // put the larger partition to stack
              lowStack[stackDepth] = lowIndex;
              highStack[stackDepth] = lowPartitionIndex - 1;
              stackDepth++;
              // sort the smaller partition first
              lowIndex = highPartitionIndex;
          } else {
              // put the larger partition to stack
              lowStack[stackDepth] = highPartitionIndex;
              highStack[stackDepth] = highIndex;
              stackDepth++;
              // sort the smaller partition first
              if (lowPartitionIndex > lowIndex) {
                  highIndex = lowPartitionIndex - 1;
              } else {
                  // all keys in the partition equals to the partition key
                  break;
              }
          }
          continue;
      }

      // Pop a range from stack
      if (stackDepth > 0) {
          stackDepth--;
          lowIndex = lowStack[stackDepth];
          highIndex = highStack[stackDepth];
          continue;
      } else return;
  }
}

static void 
BWTIncBuildRelativeRank(tmap_bwt_int_t* __restrict sortedRank, tmap_bwt_int_t* __restrict seq,
                        tmap_bwt_int_t* __restrict relativeRank, const tmap_bwt_int_t numItem,
                        tmap_bwt_int_t oldInverseSa0, const tmap_bwt_int_t *cumulativeCount)
{
  tmap_bwt_int_t i, c;
  tmap_bwt_int_t s, r;
  tmap_bwt_int_t lastRank, lastIndex;
  tmap_bwt_int_t oldInverseSa0RelativeRank = 0;
  tmap_bwt_int_t freq;

  lastIndex = numItem;
  lastRank = sortedRank[numItem];
  if (lastRank > oldInverseSa0) {
      sortedRank[numItem]--;	// to prepare for merging; $ is not encoded in bwt
  }
  s = seq[numItem];
  relativeRank[s] = numItem;
  if (lastRank == oldInverseSa0) {
      oldInverseSa0RelativeRank = numItem;
      oldInverseSa0++;	// so that this segment of code is not run again
      lastRank++;			// so that oldInverseSa0 become a sorted group with 1 item
  }

  c = ALPHABET_SIZE - 1;
  freq = cumulativeCount[c];

  for (i=numItem; i--;) {	// from numItem - 1 to 0
      r = sortedRank[i];
      if (r > oldInverseSa0) {
          sortedRank[i]--;	// to prepare for merging; $ is not encoded in bwt
      }
      s = seq[i];
      if (i < freq) {
          if (lastIndex >= freq) {
              lastRank++;	// to trigger the group across alphabet boundary to be split
          }
          c--;
          freq = cumulativeCount[c];
      }
      if (r == lastRank) {
          relativeRank[s] = lastIndex;
      } else {
          if (i == lastIndex - 1) {
              if (lastIndex < numItem && (tmap_bwt_sint_t)seq[lastIndex + 1] < 0) {
                  seq[lastIndex] = seq[lastIndex + 1] - 1;
              } else {
                  seq[lastIndex] = (tmap_bwt_int_t)-1;
              }
          }
          lastIndex = i;
          lastRank = r;
          relativeRank[s] = i;
          if (r == oldInverseSa0) {
              oldInverseSa0RelativeRank = i;
              oldInverseSa0++;	// so that this segment of code is not run again
              lastRank++;			// so that oldInverseSa0 become a sorted group with 1 item
          }
      }
  }

}

static void 
BWTIncBuildBwt(uint32_t*  seq, const tmap_bwt_int_t *relativeRank, const tmap_bwt_int_t numChar,
               const tmap_bwt_int_t *cumulativeCount)
{
  uint32_t c;
  tmap_bwt_int_t i, previousRank, currentRank;

  previousRank = relativeRank[0];

  for (i=1; i<=numChar; i++) {
      currentRank = relativeRank[i];
      c = (previousRank >= cumulativeCount[1]) + (previousRank >= cumulativeCount[2])
        + (previousRank >= cumulativeCount[3]);
      seq[currentRank] = c;
      previousRank = currentRank;
  }
}

static void 
BWTIncMergeBwt(const tmap_bwt_int_t* sortedRank, const uint32_t* oldBwt, const uint32_t *insertBwt,
               uint32_t* __restrict mergedBwt, const tmap_bwt_int_t numOldBwt, const tmap_bwt_int_t numInsertBwt)
{
  uint32_t bitsInWordMinusBitPerChar;
  tmap_bwt_int_t leftShift, rightShift;
  tmap_bwt_int_t o;
  tmap_bwt_int_t oIndex, iIndex, mIndex;
  tmap_bwt_int_t mWord, mChar, oWord, oChar;
  tmap_bwt_int_t numInsert;

  bitsInWordMinusBitPerChar = BITS_IN_WORD - BIT_PER_CHAR;

  oIndex = 0;
  iIndex = 0;
  mIndex = 0;

  mWord = 0;
  mChar = 0;

  mergedBwt[0] = 0;	// this can be cleared as merged Bwt slightly shift to the left in each iteration

  while (oIndex < numOldBwt) {

      // copy from insertBwt
      while (iIndex <= numInsertBwt && sortedRank[iIndex] <= oIndex) {
          if (sortedRank[iIndex] != 0) {	// special value to indicate that this is for new inverseSa0
              mergedBwt[mWord] |= insertBwt[iIndex] << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
              mIndex++;
              mChar++;
              if (mChar == CHAR_PER_WORD) {
                  mChar = 0;
                  mWord++;
                  mergedBwt[mWord] = 0;	// no need to worry about crossing mergedBwt boundary
              }
          }
          iIndex++;
      }

      // Copy from oldBwt to mergedBwt
      if (iIndex <= numInsertBwt) {
          o = sortedRank[iIndex];
      } else {
          o = numOldBwt;
      }
      numInsert = o - oIndex;

      oWord = oIndex / CHAR_PER_WORD;
      oChar = oIndex - oWord * CHAR_PER_WORD;
      if (oChar > mChar) {
          leftShift = (oChar - mChar) * BIT_PER_CHAR;
          rightShift = (CHAR_PER_WORD + mChar - oChar) * BIT_PER_CHAR;
          mergedBwt[mWord] = mergedBwt[mWord]
            | (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR))
            | (oldBwt[oWord+1] >> rightShift);
          oIndex += min(numInsert, CHAR_PER_WORD - mChar);
          while (o > oIndex) {
              oWord++;
              mWord++;
              mergedBwt[mWord] = (oldBwt[oWord] << leftShift) | (oldBwt[oWord+1] >> rightShift);
              oIndex += CHAR_PER_WORD;
          }
      } else if (oChar < mChar) {
          rightShift = (mChar - oChar) * BIT_PER_CHAR;
          leftShift = (CHAR_PER_WORD + oChar - mChar) * BIT_PER_CHAR;
          mergedBwt[mWord] = mergedBwt[mWord] 
            | (oldBwt[oWord] << (oChar * BIT_PER_CHAR) >> (mChar * BIT_PER_CHAR));
          oIndex += min(numInsert, CHAR_PER_WORD - mChar);
          while (o > oIndex) {
              oWord++;
              mWord++;
              mergedBwt[mWord] = (oldBwt[oWord-1] << leftShift) | (oldBwt[oWord] >> rightShift);
              oIndex += CHAR_PER_WORD;
          }
      } else { // oChar == mChar
          mergedBwt[mWord] = mergedBwt[mWord] | truncateLeft(oldBwt[oWord], mChar * BIT_PER_CHAR);
          oIndex += min(numInsert, CHAR_PER_WORD - mChar);
          while (o > oIndex) {
              oWord++;
              mWord++;
              mergedBwt[mWord] = oldBwt[oWord];
              oIndex += CHAR_PER_WORD;
          }
      }
      oIndex = o;
      mIndex += numInsert;

      // Clear the trailing garbage in mergedBwt
      mWord = mIndex / CHAR_PER_WORD;
      mChar = mIndex - mWord * CHAR_PER_WORD;
      if (mChar == 0) {
          mergedBwt[mWord] = 0;
      } else {
          mergedBwt[mWord] = truncateRight(mergedBwt[mWord], (BITS_IN_WORD - mChar * BIT_PER_CHAR));
      }

  }

  // copy from insertBwt
  while (iIndex <= numInsertBwt) {
      if (sortedRank[iIndex] != 0) {
          mergedBwt[mWord] |= insertBwt[iIndex] << (BITS_IN_WORD - (mChar + 1) * BIT_PER_CHAR);
          mIndex++;
          mChar++;
          if (mChar == CHAR_PER_WORD) {
              mChar = 0;
              mWord++;
              mergedBwt[mWord] = 0;	// no need to worry about crossing mergedBwt boundary
          }
      }
      iIndex++;
  }
}

void 
BWTClearTrailingBwtCode(tmap_bwt_gen_t *bwt)
{
  tmap_bwt_int_t bwtResidentSizeInWord;
  tmap_bwt_int_t wordIndex, offset;
  tmap_bwt_int_t i;

  bwtResidentSizeInWord = BWTResidentSizeInWord(bwt->textLength);

  wordIndex = bwt->textLength / CHAR_PER_WORD;
  offset = (bwt->textLength - wordIndex * CHAR_PER_WORD) * BIT_PER_CHAR;
  if (offset > 0) {
      bwt->bwtCode[wordIndex] = truncateRight(bwt->bwtCode[wordIndex], BITS_IN_WORD - offset);
  } else {
      if (wordIndex < bwtResidentSizeInWord) {
          bwt->bwtCode[wordIndex] = 0;
      }
  }

  for (i=wordIndex+1; i<bwtResidentSizeInWord; i++) {
      bwt->bwtCode[i] = 0;
  }
}


void 
BWTGenerateOccValueFromBwt(const uint32_t*  bwt, uint32_t* __restrict occValue,
                           tmap_bwt_int_t* __restrict occValueMajor,
                           const tmap_bwt_int_t textLength, const uint32_t*  decodeTable)
{
  tmap_bwt_int_t numberOfOccValueMajor, numberOfOccValue;
  uint32_t wordBetweenOccValue;
  tmap_bwt_int_t numberOfOccIntervalPerMajor;
  uint32_t c;
  tmap_bwt_int_t i, j;
  tmap_bwt_int_t occMajorIndex;
  tmap_bwt_int_t occIndex, bwtIndex;
  tmap_bwt_int_t sum;
  tmap_bwt_int_t tempOccValue0[ALPHABET_SIZE], tempOccValue1[ALPHABET_SIZE];

  wordBetweenOccValue = OCC_INTERVAL / CHAR_PER_WORD;

  // Calculate occValue
  // [lh3] by default: OCC_INTERVAL_MAJOR=65536, OCC_INTERVAL=256
  numberOfOccValue = (textLength + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;				// Value at both end for bi-directional encoding
  numberOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
  numberOfOccValueMajor = (numberOfOccValue + numberOfOccIntervalPerMajor - 1) / numberOfOccIntervalPerMajor;

  tempOccValue0[0] = 0;
  tempOccValue0[1] = 0;
  tempOccValue0[2] = 0;
  tempOccValue0[3] = 0;
  occValueMajor[0] = 0;
  occValueMajor[1] = 0;
  occValueMajor[2] = 0;
  occValueMajor[3] = 0;

  occIndex = 0;
  bwtIndex = 0;
  for (occMajorIndex=1; occMajorIndex<numberOfOccValueMajor; occMajorIndex++) {

      for (i=0; i<numberOfOccIntervalPerMajor/2; i++) {

          sum = 0;
          tempOccValue1[0] = tempOccValue0[0];
          tempOccValue1[1] = tempOccValue0[1];
          tempOccValue1[2] = tempOccValue0[2];
          tempOccValue1[3] = tempOccValue0[3];

          for (j=0; j<wordBetweenOccValue; j++) {
              c = bwt[bwtIndex];
              sum += decodeTable[c >> 16];
              sum += decodeTable[c & 0x0000FFFF];
              bwtIndex++;
          }
          if (!DNA_OCC_SUM_EXCEPTION(sum)) {
              tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
              tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
              tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
              tempOccValue1[3] += sum;
          } else {
              if (sum == 0x00000100) {
                  tempOccValue1[0] += 256;
              } else if (sum == 0x00010000) {
                  tempOccValue1[1] += 256;
              } else if (sum == 0x01000000) {
                  tempOccValue1[2] += 256;
              } else {
                  tempOccValue1[3] += 256;
              }
          }
          occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
          occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
          occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
          occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];
          tempOccValue0[0] = tempOccValue1[0];
          tempOccValue0[1] = tempOccValue1[1];
          tempOccValue0[2] = tempOccValue1[2];
          tempOccValue0[3] = tempOccValue1[3];
          sum = 0;

          occIndex++;

          for (j=0; j<wordBetweenOccValue; j++) {
              c = bwt[bwtIndex];
              sum += decodeTable[c >> 16];
              sum += decodeTable[c & 0x0000FFFF];
              bwtIndex++;
          }
          if (!DNA_OCC_SUM_EXCEPTION(sum)) {
              tempOccValue0[0] += (sum & 0x000000FF);	sum >>= 8;
              tempOccValue0[1] += (sum & 0x000000FF);	sum >>= 8;
              tempOccValue0[2] += (sum & 0x000000FF);	sum >>= 8;
              tempOccValue0[3] += sum;
          } else {
              if (sum == 0x00000100) {
                  tempOccValue0[0] += 256;
              } else if (sum == 0x00010000) {
                  tempOccValue0[1] += 256;
              } else if (sum == 0x01000000) {
                  tempOccValue0[2] += 256;
              } else {
                  tempOccValue0[3] += 256;
              }
          }
      }

      occValueMajor[occMajorIndex * 4 + 0] = occValueMajor[(occMajorIndex - 1) * 4 + 0] + tempOccValue0[0];
      occValueMajor[occMajorIndex * 4 + 1] = occValueMajor[(occMajorIndex - 1) * 4 + 1] + tempOccValue0[1];
      occValueMajor[occMajorIndex * 4 + 2] = occValueMajor[(occMajorIndex - 1) * 4 + 2] + tempOccValue0[2];
      occValueMajor[occMajorIndex * 4 + 3] = occValueMajor[(occMajorIndex - 1) * 4 + 3] + tempOccValue0[3];
      tempOccValue0[0] = 0;
      tempOccValue0[1] = 0;
      tempOccValue0[2] = 0;
      tempOccValue0[3] = 0;

  }

  while (occIndex < (numberOfOccValue-1)/2) {
      sum = 0;
      tempOccValue1[0] = tempOccValue0[0];
      tempOccValue1[1] = tempOccValue0[1];
      tempOccValue1[2] = tempOccValue0[2];
      tempOccValue1[3] = tempOccValue0[3];
      for (j=0; j<wordBetweenOccValue; j++) {
          c = bwt[bwtIndex];
          sum += decodeTable[c >> 16];
          sum += decodeTable[c & 0x0000FFFF];
          bwtIndex++;
      }
      if (!DNA_OCC_SUM_EXCEPTION(sum)) {
          tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
          tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
          tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
          tempOccValue1[3] += sum;
      } else {
          if (sum == 0x00000100) {
              tempOccValue1[0] += 256;
          } else if (sum == 0x00010000) {
              tempOccValue1[1] += 256;
          } else if (sum == 0x01000000) {
              tempOccValue1[2] += 256;
          } else {
              tempOccValue1[3] += 256;
          }
      }
      occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
      occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
      occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
      occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];
      tempOccValue0[0] = tempOccValue1[0];
      tempOccValue0[1] = tempOccValue1[1];
      tempOccValue0[2] = tempOccValue1[2];
      tempOccValue0[3] = tempOccValue1[3];
      sum = 0;
      occIndex++;

      for (j=0; j<wordBetweenOccValue; j++) {
          c = bwt[bwtIndex];
          sum += decodeTable[c >> 16];
          sum += decodeTable[c & 0x0000FFFF];
          bwtIndex++;
      }
      if (!DNA_OCC_SUM_EXCEPTION(sum)) {
          tempOccValue0[0] += (sum & 0x000000FF);	sum >>= 8;
          tempOccValue0[1] += (sum & 0x000000FF);	sum >>= 8;
          tempOccValue0[2] += (sum & 0x000000FF);	sum >>= 8;
          tempOccValue0[3] += sum;
      } else {
          if (sum == 0x00000100) {
              tempOccValue0[0] += 256;
          } else if (sum == 0x00010000) {
              tempOccValue0[1] += 256;
          } else if (sum == 0x01000000) {
              tempOccValue0[2] += 256;
          } else {
              tempOccValue0[3] += 256;
          }
      }
  }

  sum = 0;
  tempOccValue1[0] = tempOccValue0[0];
  tempOccValue1[1] = tempOccValue0[1];
  tempOccValue1[2] = tempOccValue0[2];
  tempOccValue1[3] = tempOccValue0[3];

  if (occIndex * 2 < numberOfOccValue - 1) {
      for (j=0; j<wordBetweenOccValue; j++) {
          c = bwt[bwtIndex];
          sum += decodeTable[c >> 16];
          sum += decodeTable[c & 0x0000FFFF];
          bwtIndex++;
      }
      if (!DNA_OCC_SUM_EXCEPTION(sum)) {
          tempOccValue1[0] += (sum & 0x000000FF);	sum >>= 8;
          tempOccValue1[1] += (sum & 0x000000FF);	sum >>= 8;
          tempOccValue1[2] += (sum & 0x000000FF);	sum >>= 8;
          tempOccValue1[3] += sum;
      } else {
          if (sum == 0x00000100) {
              tempOccValue1[0] += 256;
          } else if (sum == 0x00010000) {
              tempOccValue1[1] += 256;
          } else if (sum == 0x01000000) {
              tempOccValue1[2] += 256;
          } else {
              tempOccValue1[3] += 256;
          }
      }
  }

  occValue[occIndex * 4 + 0] = (tempOccValue0[0] << 16) | tempOccValue1[0];
  occValue[occIndex * 4 + 1] = (tempOccValue0[1] << 16) | tempOccValue1[1];
  occValue[occIndex * 4 + 2] = (tempOccValue0[2] << 16) | tempOccValue1[2];
  occValue[occIndex * 4 + 3] = (tempOccValue0[3] << 16) | tempOccValue1[3];

}

static void 
BWTIncConstruct(tmap_bwt_gen_inc_t *bwtInc, const tmap_bwt_int_t numChar)
{
  uint32_t i;
  tmap_bwt_int_t mergedBwtSizeInWord, mergedOccSizeInWord;
  uint32_t firstCharInThisIteration;

  tmap_bwt_int_t *relativeRank, *seq, *sortedRank;
  uint32_t *insertBwt, *mergedBwt;
  tmap_bwt_int_t newInverseSa0RelativeRank, oldInverseSa0RelativeRank, newInverseSa0;

#ifdef DEBUG
  if (numChar > bwtInc->buildSize) {
      tmap_error("BWTIncConstruct(): numChar > buildSize", Exit, OutOfRange);
  }
#endif

  mergedBwtSizeInWord = BWTResidentSizeInWord(bwtInc->bwt->textLength + numChar);
  mergedOccSizeInWord = BWTOccValueMinorSizeInWord(bwtInc->bwt->textLength + numChar);

  initializeVALBIG(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

  if (bwtInc->bwt->textLength == 0) {		// Initial build

      // Set address
      seq = (tmap_bwt_int_t*)bwtInc->workingMemory;
      relativeRank = seq + bwtInc->buildSize + 1;
      mergedBwt = insertBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord;	// build in place

      if((void*)(relativeRank + bwtInc->buildSize + 1) > (void*)bwtInc->packedText) {
          tmap_error("build error", Exit, OutOfRange);
      }
      if((void*)(relativeRank + bwtInc->buildSize + 1) > (void*)mergedBwt) {
          tmap_error("build error", Exit, OutOfRange);
      }

      BWTIncPutPackedTextToRank(bwtInc->packedText, relativeRank, bwtInc->cumulativeCountInCurrentBuild, numChar);

      firstCharInThisIteration = relativeRank[0];
      relativeRank[numChar] = 0;

      // Sort suffix
      QSufSortSuffixSort((tmap_bwt_sint_t*)relativeRank, (tmap_bwt_sint_t*)seq, (tmap_bwt_sint_t)numChar, (tmap_bwt_sint_t)ALPHABET_SIZE - 1, 0, FALSE);
      newInverseSa0 = relativeRank[0];

      // Clear BWT area
      initializeVAL(insertBwt, mergedBwtSizeInWord, 0);

      // Build BWT
      BWTIncBuildPackedBwt(relativeRank, insertBwt, numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->packedShift);

      // so that the cumulativeCount is not deducted
      bwtInc->firstCharInLastIteration = ALPHABET_SIZE;

  } else {		// Incremental build
      // Set address
      sortedRank = (tmap_bwt_int_t*)bwtInc->workingMemory;
      seq = sortedRank + bwtInc->buildSize + 1;
      insertBwt = (unsigned*)seq;
      relativeRank = seq + bwtInc->buildSize + 1;

      // Store the first character of this iteration
      firstCharInThisIteration = bwtInc->packedText[0] >> (BITS_IN_WORD - BIT_PER_CHAR);

      // Count occurrence of input text
      ForwardDNAAllOccCountNoLimit(bwtInc->packedText, numChar, bwtInc->cumulativeCountInCurrentBuild + 1, bwtInc->bwt->decodeTable);
      // Add the first character of the previous iteration to represent the inverseSa0 of the previous iteration
      bwtInc->cumulativeCountInCurrentBuild[bwtInc->firstCharInLastIteration + 1]++;
      bwtInc->cumulativeCountInCurrentBuild[2] += bwtInc->cumulativeCountInCurrentBuild[1];
      bwtInc->cumulativeCountInCurrentBuild[3] += bwtInc->cumulativeCountInCurrentBuild[2];
      bwtInc->cumulativeCountInCurrentBuild[4] += bwtInc->cumulativeCountInCurrentBuild[3];

      // Get rank of new suffix among processed suffix
      // The seq array is built into ALPHABET_SIZE + 2 groups; ALPHABET_SIZE groups + 1 group divided into 2 by inverseSa0 + inverseSa0 as 1 group
      oldInverseSa0RelativeRank = BWTIncGetAbsoluteRank(bwtInc->bwt, sortedRank, seq, bwtInc->packedText, 
                                                        numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->firstCharInLastIteration);

      // Sort rank by ALPHABET_SIZE + 2 groups (or ALPHABET_SIZE + 1 groups when inverseSa0 sit on the border of a group)
      for (i=0; i<ALPHABET_SIZE; i++) {
          if (bwtInc->cumulativeCountInCurrentBuild[i] > oldInverseSa0RelativeRank ||
              bwtInc->cumulativeCountInCurrentBuild[i+1] <= oldInverseSa0RelativeRank) {
              BWTIncSortKey(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], bwtInc->cumulativeCountInCurrentBuild[i+1] - bwtInc->cumulativeCountInCurrentBuild[i]);
          } else {
              if (bwtInc->cumulativeCountInCurrentBuild[i] < oldInverseSa0RelativeRank) {
                  BWTIncSortKey(sortedRank + bwtInc->cumulativeCountInCurrentBuild[i], seq + bwtInc->cumulativeCountInCurrentBuild[i], oldInverseSa0RelativeRank - bwtInc->cumulativeCountInCurrentBuild[i]);
              }
              if (bwtInc->cumulativeCountInCurrentBuild[i+1] > oldInverseSa0RelativeRank + 1) {
                  BWTIncSortKey(sortedRank + oldInverseSa0RelativeRank + 1, seq + oldInverseSa0RelativeRank + 1, bwtInc->cumulativeCountInCurrentBuild[i+1] - oldInverseSa0RelativeRank - 1);
              }
          }
      }

      // build relative rank; sortedRank is updated for merging to cater for the fact that $ is not encoded in bwt
      // the cumulative freq information is used to make sure that inverseSa0 and suffix beginning with different characters are kept in different unsorted groups)
      BWTIncBuildRelativeRank(sortedRank, seq, relativeRank, numChar, bwtInc->bwt->inverseSa0, bwtInc->cumulativeCountInCurrentBuild);
#ifdef DEBUG
      if (relativeRank[numChar] != oldInverseSa0RelativeRank) {
          tmap_error("BWTIncConstruct(): relativeRank[numChar] != oldInverseSa0RelativeRank", Exit, OutOfRange);
      }
#endif

      // Sort suffix
      QSufSortSuffixSort((tmap_bwt_sint_t*)relativeRank, (tmap_bwt_sint_t*)seq, (tmap_bwt_sint_t)numChar, (tmap_bwt_sint_t)numChar, 1, TRUE);

      newInverseSa0RelativeRank = relativeRank[0];
      newInverseSa0 = sortedRank[newInverseSa0RelativeRank] + newInverseSa0RelativeRank;

      sortedRank[newInverseSa0RelativeRank] = 0;	// a special value so that this is skipped in the merged bwt

      // Build BWT;  seq is overwritten by insertBwt
      BWTIncBuildBwt(insertBwt, relativeRank, numChar, bwtInc->cumulativeCountInCurrentBuild);

      // Merge BWT
      mergedBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord 
        - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR * sizeof(tmap_bwt_int_t) / 4;
      if(mergedBwt < insertBwt + numChar) {
          tmap_error("build error", Exit, OutOfRange);
      }
      // minus numberOfIteration * occInterval to create a buffer for merging
      BWTIncMergeBwt(sortedRank, bwtInc->bwt->bwtCode, insertBwt, mergedBwt, bwtInc->bwt->textLength, numChar);

  }

  // Build auxiliary structure and update info and pointers in BWT
  bwtInc->bwt->textLength += numChar;
  bwtInc->bwt->bwtCode = mergedBwt;
  bwtInc->bwt->bwtSizeInWord = mergedBwtSizeInWord;
  bwtInc->bwt->occSizeInWord = mergedOccSizeInWord;
  if (mergedBwt < bwtInc->workingMemory + mergedOccSizeInWord) {
      tmap_error("BWTIncConstruct() : Not enough memory allocated", Exit, OutOfRange);
  }

  bwtInc->bwt->occValue = mergedBwt - mergedOccSizeInWord;

  BWTClearTrailingBwtCode(bwtInc->bwt);
  BWTGenerateOccValueFromBwt(bwtInc->bwt->bwtCode, bwtInc->bwt->occValue, bwtInc->bwt->occValueMajor,
                             bwtInc->bwt->textLength, bwtInc->bwt->decodeTable);

  bwtInc->bwt->inverseSa0 = newInverseSa0;

  bwtInc->bwt->cumulativeFreq[1] += bwtInc->cumulativeCountInCurrentBuild[1] - (bwtInc->firstCharInLastIteration <= 0);
  bwtInc->bwt->cumulativeFreq[2] += bwtInc->cumulativeCountInCurrentBuild[2] - (bwtInc->firstCharInLastIteration <= 1);
  bwtInc->bwt->cumulativeFreq[3] += bwtInc->cumulativeCountInCurrentBuild[3] - (bwtInc->firstCharInLastIteration <= 2);
  bwtInc->bwt->cumulativeFreq[4] += bwtInc->cumulativeCountInCurrentBuild[4] - (bwtInc->firstCharInLastIteration <= 3);

  bwtInc->firstCharInLastIteration = firstCharInThisIteration;

  // Set build size and text address for the next build
  BWTIncSetBuildSizeAndTextAddr(bwtInc);
  bwtInc->numberOfIterationDone++;

}

tmap_bwt_gen_inc_t *
BWTIncConstructFromPacked(const char *inputFileName, 
                          const tmap_bwt_int_t initialMaxBuildSize, const tmap_bwt_int_t incMaxBuildSize)
{

  FILE *packedFile;
  tmap_bwt_int_t packedFileLen;
  tmap_bwt_int_t totalTextLength;
  tmap_bwt_int_t textToLoad, textSizeInByte;
  tmap_bwt_int_t processedTextLength;
  uint8_t lastByteLength;

  tmap_bwt_gen_inc_t *bwtInc;

  packedFile = (FILE*)fopen(inputFileName, "rb");

  if (packedFile == NULL) {
      tmap_error("BWTIncConstructFromPacked() : Cannot open inputFileName", Exit, OutOfRange);
  }

  fseek(packedFile, -1, SEEK_END);
  packedFileLen = ftell(packedFile);
  if ((long)packedFileLen < 0) {
      tmap_error("BWTIncConstructFromPacked: Cannot determine file length", Exit, OutOfRange);
  }
  if(1 != fread(&lastByteLength, sizeof(uint8_t), 1, packedFile)) {
      tmap_error(NULL, Exit, ReadFileError);
  }
  totalTextLength = TextLengthFromBytePacked(packedFileLen, BIT_PER_CHAR, lastByteLength);

  bwtInc = BWTIncCreate(totalTextLength, initialMaxBuildSize, incMaxBuildSize);

  BWTIncSetBuildSizeAndTextAddr(bwtInc);

  if (bwtInc->buildSize > totalTextLength) {
      textToLoad = totalTextLength;
  } else {
      textToLoad = totalTextLength - ((totalTextLength - bwtInc->buildSize + CHAR_PER_WORD - 1) / CHAR_PER_WORD * CHAR_PER_WORD);
  }
  textSizeInByte = textToLoad / CHAR_PER_BYTE;	// excluded the odd byte

  fseek(packedFile, -2, SEEK_CUR);
  fseek(packedFile, -((long)textSizeInByte), SEEK_CUR);
  if(textSizeInByte + 1 != fread(bwtInc->textBuffer, sizeof(uint8_t), textSizeInByte + 1, packedFile)) {
      tmap_error(NULL, Exit, ReadFileError);
  }
  fseek(packedFile, -((long)textSizeInByte + 1), SEEK_CUR);
  
  // base case
  ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
  BWTIncConstruct(bwtInc, textToLoad);

  // iterate
  processedTextLength = textToLoad;
  while (processedTextLength < totalTextLength) {
      textToLoad = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;
      if (textToLoad > totalTextLength - processedTextLength) {
          textToLoad = totalTextLength - processedTextLength;
      }
      textSizeInByte = textToLoad / CHAR_PER_BYTE;
      fseek(packedFile, -((long)textSizeInByte), SEEK_CUR);
      if(textSizeInByte != fread(bwtInc->textBuffer, sizeof(uint8_t), textSizeInByte, packedFile)) {
          tmap_error(NULL, Exit, ReadFileError);
      }
      fseek(packedFile, -((long)textSizeInByte), SEEK_CUR);
      ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
      BWTIncConstruct(bwtInc, textToLoad);
      processedTextLength += textToLoad;
      if (bwtInc->numberOfIterationDone % 10 == 0) {
          tmap_progress_print2("%.2lf%% complete with %llu iterations done and %llu bases processed", 
                               100.0 * processedTextLength / (double)totalTextLength,
                               (unsigned long long int)bwtInc->numberOfIterationDone, 
                               (unsigned long long int)processedTextLength);
      }
  }
  if (bwtInc->numberOfIterationDone % 10 != 0) {
      tmap_progress_print2("%.2lf%% complete with %llu iterations done and %llu bases processed", 
                           100.0 * processedTextLength / (double)totalTextLength,
                           (unsigned long long int)bwtInc->numberOfIterationDone, 
                           (unsigned long long int)processedTextLength);
  }
  return bwtInc;
}

void 
BWTFree(tmap_bwt_gen_t *bwt)
{
  if (bwt == 0) return;
  free(bwt->cumulativeFreq);
  //free(bwt->bwtCode); // not allocated?
  //free(bwt->occValue); // not allocated?
  free(bwt->occValueMajor);
  free(bwt->decodeTable);
  free(bwt);
}

void 
BWTIncFree(tmap_bwt_gen_inc_t *bwtInc)
{
  if (bwtInc == 0) return;
  BWTFree(bwtInc->bwt);
  free(bwtInc->workingMemory);
  free(bwtInc->cumulativeCountInCurrentBuild);
  free(bwtInc->packedShift);
  free(bwtInc);
}

static tmap_bwt_int_t 
BWTFileSizeInWord(const tmap_bwt_int_t numChar)
{
  // The $ in BWT at the position of inverseSa0 is not encoded
  return (numChar + CHAR_PER_WORD - 1) / CHAR_PER_WORD;
}

void 
BWTSaveBwtCodeAndOcc(tmap_bwt_t *bwt_out, const tmap_bwt_gen_t *bwt, const char *fn_fasta, int32_t occ_interval) 
{
  tmap_bwt_int_t i;
  tmap_bwt_t *bwt_tmp=NULL;

  // Move over to bwt data structure
  if(bwt_out->bwt_size != BWTFileSizeInWord(bwt->textLength)) {
      tmap_error(NULL, Exit, OutOfRange);
  }
  bwt_out->primary = bwt->inverseSa0;
  for(i=0;i<1+ALPHABET_SIZE;i++) {
      bwt_out->L2[i] = bwt->cumulativeFreq[i];
  }
  bwt_out->bwt = bwt->bwtCode; // shallow copy
  bwt_out->occ_interval = OCC_INTERVAL;

  // write
  bwt_out->hash_width = 0; // none yet
  tmap_bwt_write(fn_fasta, bwt_out);

  // free and nullify
  bwt_out->bwt = NULL;

  // update occurrence interval, if necessary
  if(occ_interval != bwt_out->occ_interval) {
      bwt_tmp = tmap_bwt_read(fn_fasta);
      tmap_bwt_update_occ_interval(bwt_tmp, occ_interval);
      tmap_bwt_write(fn_fasta, bwt_tmp);
      tmap_bwt_destroy(bwt_tmp);
  }
}

int32_t
tmap_bwt_tune_hash_width(uint64_t ref_len) 
{
  int32_t hash_width = 0;
  while(0 < ref_len) {
      ref_len >>= 2; // divide by four
      hash_width++;
  }
  if(hash_width < TMAP_BWT_HASH_WIDTH_AUTO_MIN) {
      hash_width = TMAP_BWT_HASH_WIDTH_AUTO_MIN;
  }
  else if(TMAP_BWT_HASH_WIDTH_AUTO_MAX < hash_width) {
      hash_width = TMAP_BWT_HASH_WIDTH_AUTO_MAX;
  }

  return hash_width;
}

void 
tmap_bwt_pac2bwt(const char *fn_fasta, uint32_t is_large, int32_t occ_interval, int32_t hash_width, int32_t check_hash)
{
  tmap_bwt_gen_inc_t *bwtInc=NULL;
  tmap_bwt_t *bwt=NULL;
  uint8_t *buf=NULL;
  tmap_bwt_int_t i;
  char *fn_pac=NULL;
  tmap_refseq_t *refseq=NULL;
  uint64_t ref_len;

  tmap_progress_print("constructing the BWT string from the packed FASTA");

  // read in packed FASTA
  refseq = tmap_refseq_read(fn_fasta);
  ref_len = refseq->len / 2; // NB: the reference is both forward and reverse compliment, so divide by two

  // initialization
  bwt = tmap_calloc(1, sizeof(tmap_bwt_t), "bwt");
  bwt->seq_len = refseq->len;
  bwt->bwt_size = (bwt->seq_len + 15) >> 4; // 2-bit packed
  bwt->version_id = TMAP_VERSION_ID;

  if(1 == is_large) {
      // destroy the reference sequence
      tmap_refseq_destroy(refseq);
      // create the bwt
      fn_pac = tmap_get_file_name(fn_fasta, TMAP_PAC_FILE);
      if(TMAP_PAC_COMPRESSION != TMAP_FILE_NO_COMPRESSION) { // the below uses fseek
          tmap_error("PAC compression not supported", Exit, OutOfRange);
      }
      bwtInc = BWTIncConstructFromPacked(fn_pac, 10000000, 10000000);
      BWTSaveBwtCodeAndOcc(bwt, bwtInc->bwt, fn_fasta, occ_interval);
      BWTIncFree(bwtInc);
      free(fn_pac);
  }
  else {
      // From bwtmisc.c at http://bio-bwa.sf.net

      // prepare sequence
      for(i=0;i<ALPHABET_SIZE+1;i++) {
          bwt->L2[i]=0;
      }
      buf = tmap_calloc(bwt->seq_len + 1, sizeof(uint8_t), "buf");
      for(i=0;i<bwt->seq_len;i++) {
          buf[i] = tmap_refseq_seq_i(refseq, i);
          ++bwt->L2[1+buf[i]];
      }
      for(i=2;i<ALPHABET_SIZE+1;i++) {
          bwt->L2[i] += bwt->L2[i-1];
      }
      // destroy the reference sequence
      tmap_refseq_destroy(refseq);

      // Burrows-Wheeler Transform
      bwt->primary = tmap_bwt_gen_short(buf, bwt->seq_len);
      bwt->bwt = tmap_calloc(bwt->bwt_size, sizeof(uint32_t), "bwt->bwt");
      for(i=0;i<bwt->seq_len;i++) {
          // 2-bit packing for DNA
          bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
      }
      free(buf);
      bwt->occ_interval = 1; 

      // update occurrence interval
      tmap_bwt_update_occ_interval(bwt, occ_interval);

      bwt->hash_width = 0; // none yet
      tmap_bwt_write(fn_fasta, bwt);
  }
  tmap_bwt_destroy(bwt);

  tmap_progress_print2("constructed the BWT string from the packed FASTA");

  if(INT32_MAX == hash_width) {
      hash_width = tmap_bwt_tune_hash_width(ref_len);
      tmap_progress_print2("setting the BWT hash width to %d", hash_width);
  }

  if(0 < hash_width) {
      bwt = tmap_bwt_read(fn_fasta); 
      tmap_bwt_gen_hash(bwt, hash_width, check_hash);
      tmap_bwt_write(fn_fasta, bwt);
      tmap_bwt_destroy(bwt);
  }
  else {
      tmap_progress_print("skipping occurrence hash creation");
  }
}

void 
tmap_bwt_update_hash(const char *fn_fasta, int32_t hash_width, int32_t check_hash)
{
  int32_t i;
  tmap_bwt_t *bwt;

  // read in the bwt
  bwt = tmap_bwt_read(fn_fasta); 

  // new hash width?
  if(hash_width != bwt->hash_width) {
      // free the previous hash
      if(NULL != bwt->hash_k) {
          for(i=0;i<bwt->hash_width;i++) {
              free(bwt->hash_k[i]);
          }
      }
      if(NULL != bwt->hash_l) {
          for(i=0;i<bwt->hash_width;i++) {
              free(bwt->hash_l[i]);
          }
      }
      free(bwt->hash_k);
      free(bwt->hash_l);
      bwt->hash_k = bwt->hash_l = NULL;

      // new hash
      tmap_bwt_gen_hash(bwt, hash_width, check_hash);

      // write
      tmap_bwt_write(fn_fasta, bwt);
  }
  
  // destroy
  tmap_bwt_destroy(bwt);
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

/**
 * Constructs the burrows-wheeler transformed string of a given string.
 * @param T[0..n-1] The input string.
 * @param n The length of the given string.
 * @return The primary index if no error occurred, -1 or -2 otherwise.
 */
uint32_t 
tmap_bwt_gen_short(uint8_t *T, uint32_t n)
{
  int32_t *SA=NULL;
  uint32_t i, primary = 0;
  SA = tmap_calloc(n+1, sizeof(int32_t), "SA");
  tmap_sa_gen_short(T, SA, n);

  for (i = 0; i <= n; ++i) {
      if (SA[i] == 0) primary = i;
      else SA[i] = T[SA[i] - 1];
  }
  for (i = 0; i < primary; ++i) T[i] = SA[i];
  for (; i < n; ++i) T[i] = SA[i + 1];
  free(SA);
  return primary;
}
