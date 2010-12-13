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
#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_progress.h"
#include "../util/fmap_definitions.h"
#include "../io/fmap_file.h"
#include "fmap_refseq.h"
#include "fmap_bwt.h"
#include "fmap_sa.h"
#include "fmap_bwt_gen.h"

static uint32_t 
TextLengthFromBytePacked(uint32_t bytePackedLength, uint32_t bitPerChar,
                         uint32_t lastByteLength)
{
  if (bytePackedLength > ALL_ONE_MASK / (BITS_IN_BYTE / bitPerChar)) {
      fmap_error("TextLengthFromBytePacked(): text length > 2^32", Exit, OutOfRange);
  }
  return (bytePackedLength - 1) * (BITS_IN_BYTE / bitPerChar) + lastByteLength;
}

static void 
initializeVAL(uint32_t *startAddr, const uint32_t length, const uint32_t initValue)
{
  uint32_t i;
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
static uint32_t 
BWTOccValueMajorSizeInWord(const uint32_t numChar)
{
  uint32_t numOfOccValue;
  uint32_t numOfOccIntervalPerMajor;
  numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1; // Value at both end for bi-directional encoding
  numOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
  return (numOfOccValue + numOfOccIntervalPerMajor - 1) / numOfOccIntervalPerMajor * ALPHABET_SIZE;
}

// for BWTIncCreate()
static uint32_t 
BWTOccValueMinorSizeInWord(const uint32_t numChar)
{
  uint32_t numOfOccValue;
  numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;		// Value at both end for bi-directional encoding
  return (numOfOccValue + OCC_VALUE_PER_WORD - 1) / OCC_VALUE_PER_WORD * ALPHABET_SIZE;
}

// for BWTIncCreate()
static uint32_t 
BWTResidentSizeInWord(const uint32_t numChar) {

    uint32_t numCharRoundUpToOccInterval;

    // The $ in BWT at the position of inverseSa0 is not encoded
    numCharRoundUpToOccInterval = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL;

    return (numCharRoundUpToOccInterval + CHAR_PER_WORD - 1) / CHAR_PER_WORD;

}

static void 
BWTIncSetBuildSizeAndTextAddr(fmap_bwt_gen_inc_t *bwtInc)
{
  uint32_t maxBuildSize;

  if (bwtInc->bwt->textLength == 0) {
      // initial build
      // Minus 2 because n+1 entries of seq and rank needed for n char
      maxBuildSize = (bwtInc->availableWord - 2 - OCC_INTERVAL / CHAR_PER_WORD)
        / (2 * CHAR_PER_WORD + 1) * CHAR_PER_WORD;
      if (bwtInc->initialMaxBuildSize > 0) {
          bwtInc->buildSize = min(bwtInc->initialMaxBuildSize, maxBuildSize);
      } else {
          bwtInc->buildSize = maxBuildSize;
      }
  } else {
      // Minus 3 because n+1 entries of sorted rank, seq and rank needed for n char
      // Minus numberOfIterationDone because bwt slightly shift to left in each iteration
      maxBuildSize = (bwtInc->availableWord - bwtInc->bwt->bwtSizeInWord - bwtInc->bwt->occSizeInWord - 3
                      - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR) 
        / 3;
      if (maxBuildSize < CHAR_PER_WORD) {
          fmap_error("BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction", Exit, OutOfRange);
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
      fmap_error("BWTIncSetBuildSizeAndTextAddr(): Not enough space allocated to continue construction", Exit, OutOfRange);
  }

  bwtInc->buildSize = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;

  bwtInc->packedText = bwtInc->workingMemory + 2 * (bwtInc->buildSize + 1);
  bwtInc->textBuffer = (uint8_t*)(bwtInc->workingMemory + bwtInc->buildSize + 1);

}

// for ceilLog2()
uint32_t 
leadingZero(const uint32_t input)
{
  uint32_t l;
  const static uint32_t leadingZero8bit[256] = {8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
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
                              const uint32_t textLength)
{
  uint32_t i, j, k;
  uint32_t c;
  uint32_t bitPerBytePackedChar;
  uint32_t bitPerWordPackedChar;
  uint32_t charPerWord;
  uint32_t charPerByte;
  uint32_t bytePerIteration;
  uint32_t byteProcessed = 0;
  uint32_t wordProcessed = 0;
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

fmap_bwt_gen_t *
BWTCreate(const uint32_t textLength, uint32_t *decodeTable)
{
  fmap_bwt_gen_t *bwt;

  bwt = fmap_calloc(1, sizeof(fmap_bwt_gen_t), "bwt");

  bwt->textLength = 0;
  bwt->inverseSa = 0;

  bwt->cumulativeFreq = fmap_calloc((ALPHABET_SIZE + 1), sizeof(uint32_t), "bwt->cumulativeFreq");
  initializeVAL(bwt->cumulativeFreq, ALPHABET_SIZE + 1, 0);

  bwt->bwtSizeInWord = 0;

  // Generate decode tables
  if (decodeTable == NULL) {
      bwt->decodeTable = fmap_calloc(DNA_OCC_CNT_TABLE_SIZE_IN_WORD, sizeof(uint32_t), "bwt->decodeTable");
      GenerateDNAOccCountTable(bwt->decodeTable);
  } else {
      bwt->decodeTable = decodeTable;
  }

  bwt->occValueMajor = fmap_calloc(BWTOccValueMajorSizeInWord(textLength), sizeof(uint32_t), "bwt->occValueMajor");

  bwt->occSizeInWord = 0;
  bwt->occValue = NULL;

  bwt->inverseSa = NULL;

  return bwt;
}

fmap_bwt_gen_inc_t *
BWTIncCreate(const uint32_t textLength, const float targetNBit,
             const uint32_t initialMaxBuildSize, const uint32_t incMaxBuildSize)
{
  fmap_bwt_gen_inc_t *bwtInc;
  uint32_t i;

  if (targetNBit == 0) {
      fmap_error("BWTIncCreate() : targetNBit = 0", Exit, OutOfRange);
  }

  bwtInc = fmap_calloc(1, sizeof(fmap_bwt_gen_inc_t), "bwtInc");
  bwtInc->numberOfIterationDone = 0;
  bwtInc->bwt = BWTCreate(textLength, NULL);
  bwtInc->initialMaxBuildSize = initialMaxBuildSize;
  bwtInc->incMaxBuildSize = incMaxBuildSize;
  bwtInc->targetNBit = targetNBit;
  bwtInc->cumulativeCountInCurrentBuild = fmap_calloc((ALPHABET_SIZE + 1), sizeof(uint32_t), "bwtInc->cumumlativeCountInCurrentBuild");
  initializeVAL(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

  // Build frequently accessed data
  bwtInc->packedShift = fmap_calloc(CHAR_PER_WORD, sizeof(uint32_t), "bwtInc->packedShift");
  for (i=0; i<CHAR_PER_WORD; i++) {
      bwtInc->packedShift[i] = BITS_IN_WORD - (i+1) * BIT_PER_CHAR;
  }

  bwtInc->targetTextLength = textLength;
  bwtInc->availableWord = (uint32_t)((textLength + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL / BITS_IN_WORD * bwtInc->targetNBit);
  if (bwtInc->availableWord < BWTResidentSizeInWord(textLength) + BWTOccValueMinorSizeInWord(textLength)) {
      fmap_error("BWTIncCreate() : targetNBit is too low", Exit, OutOfRange);
  }
  bwtInc->workingMemory = fmap_calloc(bwtInc->availableWord, BYTES_IN_WORD, "bwtInc->workingMemory");

  return bwtInc;

}

// for BWTIncConstruct()
static void 
BWTIncPutPackedTextToRank(const uint32_t *packedText, uint32_t* __restrict rank,
                          uint32_t* __restrict cumulativeCount, const uint32_t numChar)
{
  uint32_t i, j;
  uint32_t c, t;
  uint32_t packedMask;
  uint32_t rankIndex;
  uint32_t lastWord, numCharInLastWord;

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
ForwardDNAAllOccCountNoLimit(const uint32_t*  dna, const uint32_t index,
                             uint32_t* __restrict occCount, const uint32_t*  dnaDecodeTable)
{
  static const uint32_t truncateRightMask[16] = { 0x00000000, 0xC0000000, 0xF0000000, 0xFC000000,
      0xFF000000, 0xFFC00000, 0xFFF00000, 0xFFFC0000,
      0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00,
      0xFFFFFF00, 0xFFFFFFC0, 0xFFFFFFF0, 0xFFFFFFFC };

  uint32_t iteration, wordToCount, charToCount;
  uint32_t i, j, c;
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
              fmap_error("ForwardDNAAllOccCountNoLimit(): DNA occ sum exception", Exit, OutOfRange);
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
BWTIncBuildPackedBwt(const uint32_t *relativeRank, uint32_t* __restrict bwt, const uint32_t numChar,
                     const uint32_t *cumulativeCount, const uint32_t *packedShift) {

    uint32_t i, c, r;
    uint32_t previousRank, currentRank;
    uint32_t wordIndex, charIndex;
    uint32_t inverseSa0;

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

static inline uint32_t 
BWTOccValueExplicit(const fmap_bwt_gen_t *bwt, const uint32_t occIndexExplicit,
                    const uint32_t character)
{
  uint32_t occIndexMajor;

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

uint32_t 
BWTOccValue(const fmap_bwt_gen_t *bwt, uint32_t index, const uint32_t character) {

    uint32_t occValue;
    uint32_t occExplicitIndex, occIndex;

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

static uint32_t 
BWTIncGetAbsoluteRank(fmap_bwt_gen_t *bwt, uint32_t* __restrict absoluteRank, uint32_t* __restrict seq,
                      const uint32_t *packedText, const uint32_t numChar,
                      const uint32_t* cumulativeCount, const uint32_t firstCharInLastIteration)
{
  uint32_t saIndex;
  uint32_t lastWord;
  uint32_t packedMask;
  uint32_t i, j;
  uint32_t c, t;
  uint32_t rankIndex;
  uint32_t shift;
  uint32_t seqIndexFromStart[ALPHABET_SIZE];
  uint32_t seqIndexFromEnd[ALPHABET_SIZE];

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
BWTIncSortKey(uint32_t* __restrict key, uint32_t* __restrict seq, const uint32_t numItem)
{
#define EQUAL_KEY_THRESHOLD	4	// Partition for equal key if data array size / the number of data with equal value with pivot < EQUAL_KEY_THRESHOLD

  uint32_t lowIndex, highIndex, midIndex;
  uint32_t lowPartitionIndex, highPartitionIndex;
  uint32_t lowStack[32], highStack[32];
  uint32_t stackDepth;
  uint32_t i, j;
  uint32_t tempSeq, tempKey;
  uint32_t numberOfEqualKey;

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
BWTIncBuildRelativeRank(uint32_t* __restrict sortedRank, uint32_t* __restrict seq,
                        uint32_t* __restrict relativeRank, const uint32_t numItem,
                        uint32_t oldInverseSa0, const uint32_t *cumulativeCount)
{
  uint32_t i, c;
  uint32_t s, r;
  uint32_t lastRank, lastIndex;
  uint32_t oldInverseSa0RelativeRank = 0;
  uint32_t freq;

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
              if (lastIndex < numItem && (int)seq[lastIndex + 1] < 0) {
                  seq[lastIndex] = seq[lastIndex + 1] - 1;
              } else {
                  seq[lastIndex] = (uint32_t)-1;
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
BWTIncBuildBwt(uint32_t*  seq, const uint32_t *relativeRank, const uint32_t numChar,
               const uint32_t *cumulativeCount)
{
  uint32_t i, c;
  uint32_t previousRank, currentRank;

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
BWTIncMergeBwt(const uint32_t *sortedRank, const uint32_t* oldBwt, const uint32_t *insertBwt,
               uint32_t* __restrict mergedBwt, const uint32_t numOldBwt, const uint32_t numInsertBwt)
{
  uint32_t bitsInWordMinusBitPerChar;
  uint32_t leftShift, rightShift;
  uint32_t o;
  uint32_t oIndex, iIndex, mIndex;
  uint32_t mWord, mChar, oWord, oChar;
  uint32_t numInsert;

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
BWTClearTrailingBwtCode(fmap_bwt_gen_t *bwt)
{
  uint32_t bwtResidentSizeInWord;
  uint32_t wordIndex, offset;
  uint32_t i;

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
                           uint32_t* __restrict occValueMajor,
                           const uint32_t textLength, const uint32_t*  decodeTable)
{
  uint32_t numberOfOccValueMajor, numberOfOccValue;
  uint32_t wordBetweenOccValue;
  uint32_t numberOfOccIntervalPerMajor;
  uint32_t c;
  uint32_t i, j;
  uint32_t occMajorIndex;
  uint32_t occIndex, bwtIndex;
  uint32_t sum;
  uint32_t tempOccValue0[ALPHABET_SIZE], tempOccValue1[ALPHABET_SIZE];

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
BWTIncConstruct(fmap_bwt_gen_inc_t *bwtInc, const uint32_t numChar)
{
  uint32_t i;
  uint32_t mergedBwtSizeInWord, mergedOccSizeInWord;
  uint32_t firstCharInThisIteration;

  uint32_t *relativeRank, *seq, *sortedRank, *insertBwt, *mergedBwt;
  uint32_t newInverseSa0RelativeRank, oldInverseSa0RelativeRank, newInverseSa0;

#ifdef DEBUG
  if (numChar > bwtInc->buildSize) {
      fmap_error("BWTIncConstruct(): numChar > buildSize", Exit, OutOfRange);
  }
#endif

  mergedBwtSizeInWord = BWTResidentSizeInWord(bwtInc->bwt->textLength + numChar);
  mergedOccSizeInWord = BWTOccValueMinorSizeInWord(bwtInc->bwt->textLength + numChar);

  initializeVAL(bwtInc->cumulativeCountInCurrentBuild, ALPHABET_SIZE + 1, 0);

  if (bwtInc->bwt->textLength == 0) {		// Initial build

      // Set address
      seq = bwtInc->workingMemory;
      relativeRank = seq + bwtInc->buildSize + 1;
      mergedBwt = insertBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord;	// build in place

      BWTIncPutPackedTextToRank(bwtInc->packedText, relativeRank, bwtInc->cumulativeCountInCurrentBuild, numChar);

      firstCharInThisIteration = relativeRank[0];
      relativeRank[numChar] = 0;

      // Sort suffix
      QSufSortSuffixSort((int*)relativeRank, (int*)seq, (int)numChar, (int)ALPHABET_SIZE - 1, 0, FALSE);
      newInverseSa0 = relativeRank[0];

      // Clear BWT area
      initializeVAL(insertBwt, mergedBwtSizeInWord, 0);

      // Build BWT
      BWTIncBuildPackedBwt(relativeRank, insertBwt, numChar, bwtInc->cumulativeCountInCurrentBuild, bwtInc->packedShift);

      // so that the cumulativeCount is not deducted
      bwtInc->firstCharInLastIteration = ALPHABET_SIZE;

  } else {		// Incremental build
      // Set address
      sortedRank = bwtInc->workingMemory;
      seq = sortedRank + bwtInc->buildSize + 1;
      insertBwt = seq;
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
          fmap_error("BWTIncConstruct(): relativeRank[numChar] != oldInverseSa0RelativeRank", Exit, OutOfRange);
      }
#endif

      // Sort suffix
      QSufSortSuffixSort((int*)relativeRank, (int*)seq, (int)numChar, (int)numChar, 1, TRUE);

      newInverseSa0RelativeRank = relativeRank[0];
      newInverseSa0 = sortedRank[newInverseSa0RelativeRank] + newInverseSa0RelativeRank;

      sortedRank[newInverseSa0RelativeRank] = 0;	// a special value so that this is skipped in the merged bwt

      // Build BWT
      BWTIncBuildBwt(seq, relativeRank, numChar, bwtInc->cumulativeCountInCurrentBuild);

      // Merge BWT
      mergedBwt = bwtInc->workingMemory + bwtInc->availableWord - mergedBwtSizeInWord 
        - bwtInc->numberOfIterationDone * OCC_INTERVAL / BIT_PER_CHAR;
      // minus numberOfIteration * occInterval to create a buffer for merging
      BWTIncMergeBwt(sortedRank, bwtInc->bwt->bwtCode, insertBwt, mergedBwt, bwtInc->bwt->textLength, numChar);

  }

  // Build auxiliary structure and update info and pointers in BWT
  bwtInc->bwt->textLength += numChar;
  bwtInc->bwt->bwtCode = mergedBwt;
  bwtInc->bwt->bwtSizeInWord = mergedBwtSizeInWord;
  bwtInc->bwt->occSizeInWord = mergedOccSizeInWord;
  if (mergedBwt < bwtInc->workingMemory + mergedOccSizeInWord) {
      fmap_error("BWTIncConstruct() : Not enough memory allocated", Exit, OutOfRange);
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

fmap_bwt_gen_inc_t *
BWTIncConstructFromPacked(const char *inputFileName, const float targetNBit,
                          const uint32_t initialMaxBuildSize, const uint32_t incMaxBuildSize)
{

  FILE *packedFile;
  uint32_t packedFileLen;
  uint32_t totalTextLength;
  uint32_t textToLoad, textSizeInByte;
  uint32_t processedTextLength;
  char lastByteLength;

  fmap_bwt_gen_inc_t *bwtInc;

  packedFile = (FILE*)fopen(inputFileName, "rb");

  if (packedFile == NULL) {
      fmap_error("BWTIncConstructFromPacked() : Cannot open inputFileName", Exit, OutOfRange);
  }

  fseek(packedFile, -1, SEEK_END);
  packedFileLen = ftell(packedFile);
  if ((int)packedFileLen < 0) {
      fmap_error("BWTIncConstructFromPacked: Cannot determine file length", Exit, OutOfRange);
  }
  if(1 != fread(&lastByteLength, sizeof( char), 1, packedFile)) {
      fmap_error(NULL, Exit, ReadFileError);
  }
  totalTextLength = TextLengthFromBytePacked(packedFileLen, BIT_PER_CHAR, lastByteLength);

  bwtInc = BWTIncCreate(totalTextLength, targetNBit, initialMaxBuildSize, incMaxBuildSize);

  BWTIncSetBuildSizeAndTextAddr(bwtInc);

  if (bwtInc->buildSize > totalTextLength) {
      textToLoad = totalTextLength;
  } else {
      textToLoad = totalTextLength - ((totalTextLength - bwtInc->buildSize + CHAR_PER_WORD - 1) / CHAR_PER_WORD * CHAR_PER_WORD);
  }
  textSizeInByte = textToLoad / CHAR_PER_BYTE;	// excluded the odd byte

  fseek(packedFile, -2, SEEK_CUR);
  fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
  if(textSizeInByte + 1 != fread(bwtInc->textBuffer, sizeof( char), textSizeInByte + 1, packedFile)) {
      fmap_error(NULL, Exit, ReadFileError);
  }
  fseek(packedFile, -((int)textSizeInByte + 1), SEEK_CUR);

  ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
  BWTIncConstruct(bwtInc, textToLoad);

  processedTextLength = textToLoad;

  while (processedTextLength < totalTextLength) {
      textToLoad = bwtInc->buildSize / CHAR_PER_WORD * CHAR_PER_WORD;
      if (textToLoad > totalTextLength - processedTextLength) {
          textToLoad = totalTextLength - processedTextLength;
      }
      textSizeInByte = textToLoad / CHAR_PER_BYTE;
      fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
      if(textSizeInByte != fread(bwtInc->textBuffer, sizeof( char), textSizeInByte, packedFile)) {
          fmap_error(NULL, Exit, ReadFileError);
      }
      fseek(packedFile, -((int)textSizeInByte), SEEK_CUR);
      ConvertBytePackedToWordPacked(bwtInc->textBuffer, bwtInc->packedText, ALPHABET_SIZE, textToLoad);
      BWTIncConstruct(bwtInc, textToLoad);
      processedTextLength += textToLoad;
      if (bwtInc->numberOfIterationDone % 10 == 0) {
          fmap_progress_print2("%u iterations done with %u bases processed", bwtInc->numberOfIterationDone, processedTextLength);
      }
  }
  return bwtInc;
}

void 
BWTFree(fmap_bwt_gen_t *bwt)
{
  if (bwt == 0) return;
  free(bwt->cumulativeFreq);
  //free(bwt->bwtCode); // not allocated
  //free(bwt->occValue); // not allocated
  free(bwt->occValueMajor);
  free(bwt->inverseSa);
  free(bwt->decodeTable);
  free(bwt);
}

void 
BWTIncFree(fmap_bwt_gen_inc_t *bwtInc)
{
  if (bwtInc == 0) return;
  BWTFree(bwtInc->bwt);
  free(bwtInc->workingMemory);
  free(bwtInc->cumulativeCountInCurrentBuild);
  free(bwtInc->packedShift);
  free(bwtInc);
}

static uint32_t 
BWTFileSizeInWord(const uint32_t numChar)
{
  // The $ in BWT at the position of inverseSa0 is not encoded
  return (numChar + CHAR_PER_WORD - 1) / CHAR_PER_WORD;
}

void 
BWTSaveBwtCodeAndOcc(fmap_bwt_t *bwt_out, const fmap_bwt_gen_t *bwt, const char *fn_fasta, int32_t occ_interval, uint32_t is_rev) 
{
  uint32_t i;
  fmap_bwt_t *bwt_tmp=NULL;

  // Move over to bwt data structure
  if(bwt_out->bwt_size != BWTFileSizeInWord(bwt->textLength)) {
      fmap_error(NULL, Exit, OutOfRange);
  }
  bwt_out->primary = bwt->inverseSa0;
  for(i=0;i<1+ALPHABET_SIZE;i++) {
      bwt_out->L2[i] = bwt->cumulativeFreq[i];
  }
  bwt_out->bwt = bwt->bwtCode; // shallow copy
  bwt_out->occ_interval = OCC_INTERVAL;

  // write
  bwt_out->is_rev = is_rev;
  bwt_out->hash_width = 0; // none yet
  fmap_bwt_write(fn_fasta, bwt_out, is_rev);

  // nullify
  bwt_out->bwt = NULL;

  // update occurrence interval
  if(0 < occ_interval && bwt_out->occ_interval != occ_interval) {
      bwt_tmp = fmap_bwt_read(fn_fasta, is_rev);
      fmap_bwt_update_occ_interval(bwt_tmp, occ_interval);
      fmap_bwt_write(fn_fasta, bwt_tmp, is_rev);
      fmap_bwt_destroy(bwt_tmp);
  }
}

void 
fmap_bwt_pac2bwt(const char *fn_fasta, uint32_t is_large, int32_t occ_interval, uint32_t hash_width)
{
  // TODO: streamline occ_interval and hash creation
  fmap_bwt_gen_inc_t *bwtInc=NULL;
  fmap_bwt_t *bwt=NULL;
  uint8_t *buf=NULL;
  uint32_t i, is_rev;
  char *fn_pac=NULL;
  fmap_refseq_t *refseq=NULL;

  for(is_rev=0;is_rev<=1;is_rev++) { // forward/reverse

      if(0 == is_rev) {
          fmap_progress_print("constructing the BWT string from the packed FASTA");
      }
      else {
          fmap_progress_print("constructing the reverse BWT string from the reversed packed FASTA");
      }

      // read in packed FASTA
      refseq = fmap_refseq_read(fn_fasta, is_rev);

      // initialization
      bwt = fmap_calloc(1, sizeof(fmap_bwt_t), "bwt");
      bwt->seq_len = refseq->len;
      bwt->bwt_size = (bwt->seq_len + 15) >> 4;
      bwt->version_id = FMAP_VERSION_ID;

      if(1 == is_large) {
          fn_pac = fmap_get_file_name(fn_fasta, (0 == is_rev) ? FMAP_PAC_FILE : FMAP_REV_PAC_FILE);
          if(FMAP_PAC_COMPRESSION != FMAP_FILE_NO_COMPRESSION) { // the below uses fseek
              fmap_error("PAC compression not supported", Exit, OutOfRange);
          }
          bwtInc = BWTIncConstructFromPacked(fn_pac, 2.5, 10000000, 10000000);
          BWTSaveBwtCodeAndOcc(bwt, bwtInc->bwt, fn_fasta, occ_interval, is_rev);
          BWTIncFree(bwtInc);
          free(fn_pac);
      }
      else {
          // From bwtmisc.c at http://bio-bwa.sf.net

          // prepare sequence
          for(i=0;i<ALPHABET_SIZE+1;i++) {
              bwt->L2[i]=0;
          }
          buf = fmap_calloc(bwt->seq_len + 1, sizeof(uint8_t), "buf");
          for(i=0;i<bwt->seq_len;i++) {
              buf[i] = fmap_refseq_seq_i(refseq, i);
              ++bwt->L2[1+buf[i]];
          }
          for(i=2;i<ALPHABET_SIZE+1;i++) {
            bwt->L2[i] += bwt->L2[i-1];
          }

          // Burrows-Wheeler Transform
          bwt->primary = fmap_bwt_gen_short(buf, bwt->seq_len);
          bwt->bwt = fmap_calloc(bwt->bwt_size, sizeof(uint32_t), "bwt->bwt");
          for(i=0;i<bwt->seq_len;i++) {
              bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
          }
          free(buf);
          bwt->occ_interval = 1; 

          // update occurrence interval
          if(0 < occ_interval && bwt->occ_interval != occ_interval) {
              fmap_bwt_update_occ_interval(bwt, occ_interval);
          }

          bwt->is_rev = is_rev;
          bwt->hash_width = 0; // none yet
          fmap_bwt_write(fn_fasta, bwt, is_rev);
      }
      fmap_bwt_destroy(bwt);
      fmap_refseq_destroy(refseq);

      if(0 == is_rev) {
          fmap_progress_print2("constructed the BWT string from the packed FASTA");
      }
      else {
          fmap_progress_print2("constructed the reverse BWT string from the reversed packed FASTA");
      }
      
      bwt = fmap_bwt_read(fn_fasta, is_rev);
      fmap_bwt_gen_hash(bwt, hash_width);
      fmap_bwt_write(fn_fasta, bwt, is_rev);
      fmap_bwt_destroy(bwt);
  }
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
fmap_bwt_gen_short(uint8_t *T, uint32_t n)
{
  int32_t *SA=NULL;
  uint32_t i, primary = 0;
  SA = fmap_calloc(n+1, sizeof(int32_t), "SA");
  fmap_sa_gen_short(T, SA, n);

  for (i = 0; i <= n; ++i) {
      if (SA[i] == 0) primary = i;
      else SA[i] = T[SA[i] - 1];
  }
  for (i = 0; i < primary; ++i) T[i] = SA[i];
  for (; i < n; ++i) T[i] = SA[i + 1];
  free(SA);
  return primary;
}
