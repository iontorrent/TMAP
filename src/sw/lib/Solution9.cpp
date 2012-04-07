// Coder: Blazde
// Submission: 9
// URL http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=7212791&subnum=9
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <pmmintrin.h>
#include <cstdio>
#include "solution9.h"

using namespace std;

#define NDEBUG
#include <cassert>

int Solution9Reference::process(string b, string a, int qsc, int qec, int mm, int mi, int o, int e, int dir,
                                int *_opt, int *_te, int *_qe, int *_n_best) {
    int n = b.size(), m = a.size();

    int id = -1;
    if (qsc == 1 && qec == 1) id = 1;
    if (qsc == 0 && qec == 1) id = 2;
    if (qsc == 1 && qec == 0) id = 3;
    if (qsc == 0 && qec == 0) id = 4;

    if (id == 1 || id == 3) {
        for (int i=0; i <= m; i++) M[i][0] = 0;
        for (int j=0; j <= n; j++) M[0][j] = 0;
        for (int i=0; i <= m; i++) H[i][0] = V[i][0] = -INF;
        for (int j=0; j <= n; j++) H[0][j] = V[0][j] = -INF;

        for (int i=1; i <= m; i++)
          for (int j=1; j <= n; j++) {
              V[i][j] = max(M[i-1][j] + o + e, V[i-1][j] + e);
              H[i][j] = max(M[i][j-1] + o + e, H[i][j-1] + e);
              int mx = max(max(M[i-1][j-1], V[i-1][j-1]), H[i-1][j-1]);
              M[i][j] = max(0, mx + (a[i-1] == b[j-1] ? mm : mi));
          }
    } else {
        for (int j=0; j <= n; j++) M[0][j] = 0;
        for (int i=1; i <= m; i++) M[i][0] = -INF;
        for (int j=0; j <= n; j++) H[0][j] = V[0][j] = -INF;
        for (int i=1; i <= m; i++) {
            V[i][0] = -INF;
            H[i][0] = o + e * i;
        }

        for (int i=1; i <= m; i++)
          for (int j=1; j <= n; j++) {
              V[i][j] = max(M[i-1][j] + o + e, V[i-1][j] + e);
              H[i][j] = max(M[i][j-1] + o + e, H[i][j-1] + e);
              int mx = max(max(M[i-1][j-1], V[i-1][j-1]), H[i-1][j-1]);
              M[i][j] = mx + (a[i-1] == b[j-1] ? mm : mi);
          }
    }

    int minI = 1, maxI = m, minJ = 1, maxJ = n;
    if (id == 3 || id == 4) minI = m;

    int opt = -INF, query_end = -1, target_end = -1, n_best = 0;

    for (int i=minI; i <= maxI; i++)
      for (int j=minJ; j <= maxJ; j++) {
          opt = max(opt, M[i][j]);
          opt = max(opt, V[i][j]);
          opt = max(opt, H[i][j]);
      }

    n_best = 0;
    for (int i=minI; i <= maxI; i++)
      for (int j=minJ; j <= maxJ; j++)
        if (M[i][j] == opt || V[i][j] == opt || H[i][j] == opt)
          n_best++;

    if (dir == 0) {
        for (int i=minI; i <= maxI && query_end == -1; i++)
          for (int j=minJ; j <= maxJ && query_end == -1; j++)
            if (M[i][j] == opt || V[i][j] == opt || H[i][j] == opt) {
                query_end = i-1;
                target_end = j-1;
            }
    } else {
        for (int i=maxI; i >= minI && query_end == -1; i--)
          for (int j=maxJ; j >= minJ && query_end == -1; j--)
            if (M[i][j] == opt || V[i][j] == opt || H[i][j] == opt) {
                query_end = i-1;
                target_end = j-1;
            }
    }

    /*
    ostringstream oss;
    //oss << opt << " " << query_end << " " << target_end << " " << n_best;
    return oss.str();
    */
    (*_opt) = opt;
    (*_qe) = query_end;
    (*_te) = target_end;
    (*_n_best) = n_best;
    return opt;
}


#define MAXTARGETLENGTH 1024
#define MAXQUERYLENGTH 512

#define MINIMUMASSUMEDOPT 7

#define INFINITE16 15000
#define SIMDMULTIPLE16 8
#define SIMDLITTLESHIFT16 2
#define SIMDBIGSHIFT16 14

#define INFINITE8 255
#define SIMDMULTIPLE8 16
#define SIMDLITTLESHIFT8 1
#define SIMDBIGSHIFT8 15


const int CharMap[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 0 /* 'A' */, -1, 1 /* 'C' */, -1, -1, -1, 2 /* 'G' */, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, 3 /* 'T' */, -1,
};

#define SIMDTYPEX __m128i
#define SIMDTYPEMULTIPLE (sizeof(__m128i) / sizeof(SIMDTYPEX))
#define SIMDTYPE SIMDTYPEX __attribute__((aligned(128)))
#define SIMDARRAYWIDTH (SIMDTYPEMULTIPLE * (MAXTARGETLENGTH / SIMDMULTIPLE16 + 4))


SIMDTYPE MData[2][SIMDARRAYWIDTH];
SIMDTYPE HData[2][SIMDARRAYWIDTH];
SIMDTYPE VData[2][SIMDARRAYWIDTH];
SIMDTYPE TargetLookupData[4][SIMDARRAYWIDTH];

typedef __m128i (*MHVTRealType)[SIMDARRAYWIDTH / SIMDTYPEMULTIPLE];

static int QueryLookup[MAXQUERYLENGTH];
static int WorkQueue[SIMDARRAYWIDTH / SIMDTYPEMULTIPLE + 1];
//SIMDTYPE QueueNewMPlusOEData[SIMDARRAYWIDTH];


char __attribute__((aligned(128))) alignmentForce[128];    // Following consts don't get aligned to 16 without this hmm
const __m128i NegInfiniteSIMD16 = _mm_set1_epi16(-INFINITE16);
const __m128i TwiceInfinitesSIMD16 = _mm_set1_epi16(INFINITE16 * 2);
//const __m128i negInfiniteSIMDShiftedBig16 = _mm_srli_si128(negInfiniteSIMD16, SIMDBIGSHIFT16);    // Should be this but it triggers a bug in old GCC
const char __attribute__((aligned (16))) NegInfiniteSIMDShiftedBigData16[16] = {0x68, 0xC5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const short int* NegInfiniteSIMDShiftedBigDataShort16 = ((short int*)(NegInfiniteSIMDShiftedBigData16));
const __m128i NegInfiniteSIMDShiftedBig16 = *((__m128i*)(NegInfiniteSIMDShiftedBigData16));


Solution9::Solution9() {
    assert(NegInfiniteSIMDShiftedBigDataShort16[0] == -INFINITE16);
}

Solution9::~Solution9() {
}

int GetVariant(int queryStartClip, int queryEndClip) {
    return (4 - (queryEndClip << 1) - queryStartClip);
}
    
template<bool DIRECTION>
void inline Solution9::CheckOpt8(const unsigned char m, const int i, const int j) {
    assert(Opt >= 0 && Opt < 255);
    if(j > 0 && j <= TargetSize) {
        if(m > Opt) {
            Opt = m;
            NBest = 1;
            QueryEnd = i;
            TargetEnd = j;
        } else if(m == Opt) {
            NBest++;

            if(DIRECTION == true) {
                QueryEnd = i;
                TargetEnd = j;
            }
        }
    }
}


template<bool DIRECTION>
void inline Solution9::CheckOpt16(const short int m, const int i, const int j) {            
    if(j > 0 && j <= TargetSize) {
        if(m > Opt) {
            Opt = m;
            NBest = 1;
            QueryEnd = i;
            TargetEnd = j;
        } else if(m == Opt) {
            NBest++;

            if(DIRECTION == true) {
                QueryEnd = i;
                TargetEnd = j;
            }
        }
    }
}

template<bool DIRECTION>
void Solution9::ExtractOpts8(const __m128i newM, const int ltMask, const int i, const int jStart) {
    if((ltMask & 0x00FF) != 0x0000) {
        if(ltMask & 0x0001) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 0) & 0xFF, i, jStart + 0);
        if(ltMask & 0x0002) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 0) >> 8,   i, jStart + 1);
        if(ltMask & 0x0004) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 1) & 0xFF, i, jStart + 2);
        if(ltMask & 0x0008) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 1) >> 8,   i, jStart + 3);
        if(ltMask & 0x0010) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 2) & 0xFF, i, jStart + 4);
        if(ltMask & 0x0020) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 2) >> 8,   i, jStart + 5);
        if(ltMask & 0x0040) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 3) & 0xFF, i, jStart + 6);
        if(ltMask & 0x0080) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 3) >> 8,   i, jStart + 7);
    }
    if((ltMask & 0xFF00) != 0x0000) {
        if(ltMask & 0x0100) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 4) & 0xFF, i, jStart + 8);
        if(ltMask & 0x0200) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 4) >> 8,   i, jStart + 9);
        if(ltMask & 0x0400) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 5) & 0xFF, i, jStart + 10);
        if(ltMask & 0x0800) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 5) >> 8,   i, jStart + 11);
        if(ltMask & 0x1000) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 6) & 0xFF, i, jStart + 12);
        if(ltMask & 0x2000) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 6) >> 8,   i, jStart + 13);
        if(ltMask & 0x4000) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 7) & 0xFF, i, jStart + 14);
        if(ltMask & 0x8000) CheckOpt8<DIRECTION>(_mm_extract_epi16(newM, 7) >> 8,   i, jStart + 15);
    }    
}

template<bool DIRECTION>
void Solution9::ExtractOpts16(const __m128i newM, const int ltMask, const int i, const int jStart) {
    if((ltMask & 0x00FF) != 0x00FF) {
        if(!(ltMask & 0x0001)) CheckOpt16<DIRECTION>(_mm_extract_epi16(newM, 0), i, jStart + 0);
        if(!(ltMask & 0x0004)) CheckOpt16<DIRECTION>(_mm_extract_epi16(newM, 1), i, jStart + 1);
        if(!(ltMask & 0x0010)) CheckOpt16<DIRECTION>(_mm_extract_epi16(newM, 2), i, jStart + 2);
        if(!(ltMask & 0x0040)) CheckOpt16<DIRECTION>(_mm_extract_epi16(newM, 3), i, jStart + 3);
    }
    if((ltMask & 0xFF00) != 0xFF00) {
        if(!(ltMask & 0x0100)) CheckOpt16<DIRECTION>(_mm_extract_epi16(newM, 4), i, jStart + 4);
        if(!(ltMask & 0x0400)) CheckOpt16<DIRECTION>(_mm_extract_epi16(newM, 5), i, jStart + 5);
        if(!(ltMask & 0x1000)) CheckOpt16<DIRECTION>(_mm_extract_epi16(newM, 6), i, jStart + 6);
        if(!(ltMask & 0x4000)) CheckOpt16<DIRECTION>(_mm_extract_epi16(newM, 7), i, jStart + 7);
    }    
}









template<int VARIANT, bool DIRECTION, bool DETECTREDUNDANT> void Solution9::DoVariant8(int threshold) {
#ifndef NDEBUG
    //        char tag[128];
    //        sprintf(tag, "DoVariant8 v%d d%d r%d", VARIANT, (int)DIRECTION, (int)DETECTREDUNDANT);
#endif

    if(VARIANT == 2 || VARIANT == 4) {
        assert(false);
    }
    assert(TargetSizeInSIMDBlocks >= 1);
    assert(GapOpen <= 0);
    assert(GapExtension <= 0);


    MHVTRealType M = MData;
    MHVTRealType V = VData;
    MHVTRealType H = HData;
    MHVTRealType TargetLookup = TargetLookupData;
    //__m128i* QueueNewMPlusOE = (__m128i*)QueueNewMPlusOEData;
    __m128i optSIMD = _mm_set1_epi8(Opt & 0xFF);


    const __m128i negGapExtensionSIMD = _mm_set1_epi8(-GapExtension);
    const __m128i negGapOpenPlusgapExtensionSIMD = _mm_set1_epi8(-(GapOpen + GapExtension));
    const __m128i negGapExtensionTimes16SIMD = _mm_set1_epi8(-GapExtension * 16);
    const __m128i negGapExtensionTimes8SIMD = _mm_set1_epi8(-GapExtension * 8);
    const __m128i negGapExtensionTimes4SIMD = _mm_set1_epi8(-GapExtension * 4);
    const __m128i negGapExtensionTimes2SIMD = _mm_set1_epi8(-GapExtension * 2);
    const __m128i negMismatchSIMD = _mm_set1_epi8(-MismatchScore);
    const __m128i thresholdSIMD = _mm_set1_epi8((char)(threshold - 1));

    int jBlockSetup = 0;
    do {
        _mm_store_si128(M[0] + jBlockSetup, _mm_setzero_si128());
        _mm_store_si128(H[0] + jBlockSetup, _mm_setzero_si128());    // -Infinite = zero for 8bit version
        _mm_store_si128(V[0] + jBlockSetup, _mm_setzero_si128());
    } while(++jBlockSetup < (TargetSizeInSIMDBlocks + 2));


    for(int i = 1; i <= QuerySize; ++i) {
        int jSecondBlock = 0;
        do {
            __m128i newM = _mm_load_si128(M[(i - 1) & 0x1] + jSecondBlock);
            __m128i newM2 = _mm_load_si128(M[(i - 1) & 0x1] + jSecondBlock + 1);
            __m128i newMPlusOE = _mm_subs_epu8(newM, negGapOpenPlusgapExtensionSIMD);    // Unsigned saturared subtraction
            __m128i newMPlusOE2 = _mm_subs_epu8(newM2, negGapOpenPlusgapExtensionSIMD);

            __m128i oldV = _mm_load_si128(V[(i - 1) & 0x1] + jSecondBlock);
            __m128i oldV2 = _mm_load_si128(V[(i - 1) & 0x1] + jSecondBlock + 1);
            oldV = _mm_subs_epu8(oldV, negGapExtensionSIMD);
            oldV2 = _mm_subs_epu8(oldV2, negGapExtensionSIMD);
            __m128i newV = _mm_max_epu8(oldV, newMPlusOE);
            __m128i newV2 = _mm_max_epu8(oldV2, newMPlusOE2);
            _mm_store_si128(V[i & 0x1] + jSecondBlock, newV);
            _mm_store_si128(V[i & 0x1] + jSecondBlock + 1, newV2);

            if(DETECTREDUNDANT == true) {
                // XXX do only when needed
                _mm_store_si128(H[i & 0x1] + jSecondBlock, _mm_setzero_si128());
                _mm_store_si128(H[i & 0x1] + jSecondBlock + 1, _mm_setzero_si128());
            }

            jSecondBlock += 2;
        } while(jSecondBlock < TargetSizeInSIMDBlocks + 1);
        __m128i fixupV = _mm_load_si128(V[i & 0x1] + 0);
        int e = _mm_extract_epi16(fixupV, 0);    // No 8 bit extract...
        e = (e & 0xFF00);    // NEGATIVEINFINITE8 == 0
        fixupV = _mm_insert_epi16(fixupV, e, 0);
        _mm_store_si128(V[i & 0x1] + 0, fixupV);



        int queueLength;
        if(DETECTREDUNDANT == true) {
            queueLength = 0;
        }
        int currentQueryL = QueryLookup[i - 1];
        assert((currentQueryL & 0xFFFFFFFC) == 0);
        const __m128i* targetLookupCurrent = TargetLookup[currentQueryL];

        __m128i rollingNewM;
        if(VARIANT == 1 || VARIANT == 3) {
            rollingNewM = _mm_setzero_si128();
        }
        if(VARIANT == 2 || VARIANT == 4) {
            assert(false);
            //rollingNewM = NegInfiniteSIMDShiftedBig16;
        }


        int jInnerBlock = 0;
        do {
            __m128i thisV = _mm_load_si128(V[(i - 1) & 0x1] + jInnerBlock);
            __m128i thisH = _mm_load_si128(H[(i - 1) & 0x1] + jInnerBlock);
            __m128i thisM = _mm_load_si128(M[(i - 1) & 0x1] + jInnerBlock);
            __m128i mx = _mm_max_epu8(thisV, thisH);
            mx = _mm_max_epu8(mx, thisM);
            __m128i matchesSubMismatchScore = _mm_load_si128(targetLookupCurrent + jInnerBlock);    // Offset by the mismatch so it's all unsigned

            __m128i newM = _mm_adds_epu8(mx, matchesSubMismatchScore);
            newM = _mm_subs_epu8(newM, negMismatchSIMD);    // Unoffset again, there's a potential problem here if matches caused some saturation at the top end which is now subtracted from

            __m128i thisV2 = _mm_load_si128(V[(i - 1) & 0x1] + jInnerBlock + 1);
            __m128i thisH2 = _mm_load_si128(H[(i - 1) & 0x1] + jInnerBlock + 1);
            __m128i thisM2 = _mm_load_si128(M[(i - 1) & 0x1] + jInnerBlock  + 1);                
            __m128i mx2 = _mm_max_epu8(thisV2, thisH2);
            mx2 = _mm_max_epu8(mx2, thisM2);
            __m128i matchesSubMismatchScore2 = _mm_load_si128(targetLookupCurrent + jInnerBlock + 1);

            __m128i newM2 = _mm_adds_epu8(mx2, matchesSubMismatchScore2);
            newM2 = _mm_subs_epu8(newM2, negMismatchSIMD);

            if(VARIANT == 1 || VARIANT == 3) {
                newM = _mm_max_epu8(newM, _mm_setzero_si128());
                newM2 = _mm_max_epu8(newM2, _mm_setzero_si128());
            }

            __m128i toStoreM = _mm_slli_si128(newM, SIMDLITTLESHIFT8);
            toStoreM = _mm_or_si128(toStoreM, rollingNewM);
            rollingNewM = _mm_srli_si128(newM, SIMDBIGSHIFT8);
            _mm_store_si128(M[i & 0x1] + jInnerBlock, toStoreM);
            __m128i toStoreM2 = _mm_slli_si128(newM2, SIMDLITTLESHIFT8);
            toStoreM2 = _mm_or_si128(toStoreM2, rollingNewM);

            __m128i checkM;
            if(DETECTREDUNDANT == true) {
                __m128i checkMA = _mm_max_epu8(toStoreM, toStoreM2);
                // Signed comparison only (grr) so it can be...
                checkM = _mm_cmpgt_epi8(checkMA, negGapOpenPlusgapExtensionSIMD);    // Greater than a small +ve constant or...
                __m128i checkMSignedOtherWay = _mm_cmplt_epi8(checkMA, _mm_setzero_si128()); // >=128 and hence less than 0 in the signed world
                checkM = _mm_or_si128(checkM, checkMSignedOtherWay);
            }

            rollingNewM = _mm_srli_si128(newM2, SIMDBIGSHIFT8);
            _mm_store_si128(M[i & 0x1] + jInnerBlock +  1, toStoreM2);

            if(DETECTREDUNDANT == true) {
                int mask = _mm_movemask_epi8(checkM);
                //_mm_store_si128(QueueNewMPlusOE + queueLength, toStoreMPlusoe);
                WorkQueue[queueLength] = jInnerBlock;
                queueLength += (mask > 0);
            }
            jInnerBlock += 2;
        } while(jInnerBlock < TargetSizeInSIMDBlocks + 1);

        if(DETECTREDUNDANT == true) {
            WorkQueue[queueLength] = TargetSizeInSIMDBlocks + 1;    // So we can look one past the end below without checking length
        }



        __m128i rollingNewHZeroes;
        if(VARIANT == 1 || VARIANT == 3) {
            rollingNewHZeroes = _mm_setzero_si128();
        }


        bool continueThirdBlock;
        int jThirdBlock;
        int q;
        if(DETECTREDUNDANT == true) {
            q = 0;
            jThirdBlock = WorkQueue[q];
            continueThirdBlock = (q < queueLength);
        }
        if(DETECTREDUNDANT == false) {
            jThirdBlock = 0;
            continueThirdBlock = true;
        }

        while(continueThirdBlock) {
            //// Load dependencies

            __m128i newM = _mm_load_si128(M[i & 0x1] + jThirdBlock);
            __m128i newMPlusOE = _mm_subs_epu8(newM, negGapOpenPlusgapExtensionSIMD);
            __m128i newM2 = _mm_load_si128(M[i & 0x1] + jThirdBlock + 1);
            __m128i newMPlusOE2 = _mm_subs_epu8(newM2, negGapOpenPlusgapExtensionSIMD);

            int ltMask, ltMask2;
            if(VARIANT == 1 || VARIANT == 2) {
                // Roundabout signed less-than comparison
                __m128i optCheck = _mm_min_epu8(newM, optSIMD);
                __m128i optCheck2 = _mm_min_epu8(newM2, optSIMD);
                optCheck = _mm_cmpeq_epi8(optSIMD, optCheck);
                optCheck2 = _mm_cmpeq_epi8(optSIMD, optCheck2);

                ltMask = _mm_movemask_epi8(optCheck);
                ltMask2 = _mm_movemask_epi8(optCheck2);
            }
            if(VARIANT == 3 || VARIANT == 4) {
                // Still need to check for saturation
                __m128i optCheck = _mm_min_epu8(newM, thresholdSIMD);
                __m128i optCheck2 = _mm_min_epu8(newM2, thresholdSIMD);
                optCheck = _mm_cmpeq_epi8(optCheck, thresholdSIMD);
                optCheck2 = _mm_cmpeq_epi8(optCheck2, thresholdSIMD);

                ltMask = _mm_movemask_epi8(optCheck);
                ltMask2 = _mm_movemask_epi8(optCheck2);
            }

            __m128i rollingNewH = rollingNewHZeroes;

            __m128i oldH = _mm_subs_epu8(rollingNewH, negGapExtensionSIMD);
            __m128i newH = _mm_max_epu8(oldH, newMPlusOE);

            __m128i newH8Shift = _mm_subs_epu8(newH, negGapExtensionTimes16SIMD);
            __m128i newH2 = _mm_max_epu8(newMPlusOE2, newH8Shift);


            //// Propogate Hs

            __m128i newH4Shift = _mm_subs_epu8(newH, negGapExtensionTimes8SIMD);
            __m128i newH4Shift2 = _mm_subs_epu8(newH2, negGapExtensionTimes8SIMD);
            __m128i inbetweenH4Shift = _mm_srli_si128(newH4Shift, 8);
            newH4Shift = _mm_slli_si128(newH4Shift, 8);
            newH4Shift2 = _mm_slli_si128(newH4Shift2, 8);
            newH = _mm_max_epu8(newH, newH4Shift);
            newH2 = _mm_max_epu8(newH2, newH4Shift2);
            newH2 = _mm_max_epu8(newH2, inbetweenH4Shift);

            __m128i newH2Shift = _mm_subs_epu8(newH, negGapExtensionTimes4SIMD);
            __m128i newH2Shift2 = _mm_subs_epu8(newH2, negGapExtensionTimes4SIMD);
            __m128i inbetweenH2Shift = _mm_srli_si128(newH2Shift, 12);
            newH2Shift = _mm_slli_si128(newH2Shift, 4);
            newH2Shift2 = _mm_slli_si128(newH2Shift2, 4);
            newH = _mm_max_epu8(newH, newH2Shift);
            newH2 = _mm_max_epu8(newH2, newH2Shift2);
            newH2 = _mm_max_epu8(newH2, inbetweenH2Shift);

            __m128i newH1Shift = _mm_subs_epu8(newH, negGapExtensionTimes2SIMD);
            __m128i newH1Shift2 = _mm_subs_epu8(newH2, negGapExtensionTimes2SIMD);
            __m128i inbetweenH1Shift = _mm_srli_si128(newH1Shift, 14);
            newH1Shift = _mm_slli_si128(newH1Shift, 2);
            newH1Shift2 = _mm_slli_si128(newH1Shift2, 2);
            newH = _mm_max_epu8(newH, newH1Shift);
            newH2 = _mm_max_epu8(newH2, newH1Shift2);
            newH2 = _mm_max_epu8(newH2, inbetweenH1Shift);


            __m128i newH05Shift = _mm_subs_epu8(newH, negGapExtensionSIMD);
            __m128i newH05Shift2 = _mm_subs_epu8(newH2, negGapExtensionSIMD);
            __m128i inbetweenH05Shift = _mm_srli_si128(newH05Shift, 15);
            newH05Shift = _mm_slli_si128(newH05Shift, 1);
            newH05Shift2 = _mm_slli_si128(newH05Shift2, 1);
            newH = _mm_max_epu8(newH, newH05Shift);
            newH2 = _mm_max_epu8(newH2, newH05Shift2);
            newH2 = _mm_max_epu8(newH2, inbetweenH05Shift);


            //// Store new Hs

            __m128i toStoreH = _mm_slli_si128(newH, SIMDLITTLESHIFT8);
            __m128i toStoreH2 = _mm_slli_si128(newH2, SIMDLITTLESHIFT8);
            toStoreH = _mm_or_si128(toStoreH, rollingNewHZeroes);
            rollingNewHZeroes = _mm_srli_si128(newH, SIMDBIGSHIFT8);
            toStoreH2 = _mm_or_si128(toStoreH2, rollingNewHZeroes);
            rollingNewHZeroes = _mm_srli_si128(newH2, SIMDBIGSHIFT8);;


            _mm_store_si128(H[i & 0x1] + jThirdBlock, toStoreH);
            _mm_store_si128(H[i & 0x1] + jThirdBlock + 1, toStoreH2);


            if(VARIANT == 1 || VARIANT == 2) {
                if((ltMask | ltMask2) != 0x0000) {
                    // Order is important because of DIRECTION
                    if(ltMask != 0x0000) ExtractOpts8<DIRECTION>(newM, ltMask, i, (jThirdBlock * SIMDMULTIPLE8));
                    if(ltMask2 != 0x0000) ExtractOpts8<DIRECTION>(newM2, ltMask2, i, ((jThirdBlock + 1) * SIMDMULTIPLE8));
                    assert(Opt >= 0 && Opt < 256);
                    optSIMD = _mm_set1_epi8(Opt & 0xFF);
                }
            }
            if(VARIANT == 3 || VARIANT == 4) {
                if(ltMask != 0x0) {
                    int j = (jThirdBlock * SIMDMULTIPLE8);
                    while((ltMask & 0x1) == 0x0) {
                        j++;
                        ltMask = ltMask >> 1;
                    }
                    if(j <= TargetSize) {
                        assert(j > 0);
                        Opt = 255;
                        return;
                    }
                }
                if(ltMask2 != 0x0) {
                    int j = ((jThirdBlock + 1) * SIMDMULTIPLE8);
                    while((ltMask2 & 0x1) == 0x0) {
                        j++;
                        ltMask2 = ltMask2 >> 1;
                    }
                    if(j <= TargetSize) {
                        assert(j > 0);
                        Opt = 255;
                        return;
                    }
                }
            }

            if(DETECTREDUNDANT == true) {
                int nextH = _mm_extract_epi16(rollingNewHZeroes, 0) & 0xFF;
                q++;
                int nextJThirdBlock = jThirdBlock + 2;
                jThirdBlock = WorkQueue[q];
                continueThirdBlock = (q < queueLength);
                if(__builtin_expect(nextH > 0, 0)) {
                    if(jThirdBlock > nextJThirdBlock) {
                        jThirdBlock = nextJThirdBlock;
                        continueThirdBlock = (jThirdBlock < TargetSizeInSIMDBlocks + 1);
                        q--;
                    }
                }
            }
            if(DETECTREDUNDANT == false) {
                jThirdBlock += 2;
                continueThirdBlock = (jThirdBlock < TargetSizeInSIMDBlocks + 1);
            }
        }
    }

    if(VARIANT == 3 || VARIANT == 4) {
        int jOptCheckBlock = 0;
        do {
            __m128i finalM = _mm_load_si128(M[QuerySize & 0x1] + jOptCheckBlock);
            __m128i finalV = _mm_load_si128(V[QuerySize & 0x1] + jOptCheckBlock);
            __m128i finalH = _mm_load_si128(H[QuerySize & 0x1] + jOptCheckBlock);
            __m128i maxMVH = _mm_max_epu8(_mm_max_epu8(finalM, finalV), finalH);

            __m128i finalM2 = _mm_load_si128(M[QuerySize & 0x1] + jOptCheckBlock + 1);
            __m128i finalV2 = _mm_load_si128(V[QuerySize & 0x1] + jOptCheckBlock + 1);
            __m128i finalH2 = _mm_load_si128(H[QuerySize & 0x1] + jOptCheckBlock + 1);
            __m128i maxMVH2 = _mm_max_epu8(_mm_max_epu8(finalM2, finalV2), finalH2);


            __m128i optCheck = _mm_min_epu8(maxMVH, optSIMD);
            __m128i optCheck2 = _mm_min_epu8(maxMVH2, optSIMD);

            optCheck = _mm_cmpeq_epi8(optSIMD, optCheck);
            optCheck2 = _mm_cmpeq_epi8(optSIMD, optCheck2);

            int ltMask = _mm_movemask_epi8(optCheck);
            int ltMask2 = _mm_movemask_epi8(optCheck2);

            if((ltMask | ltMask2) != 0x0000) {
                // Order is important because of DIRECTION
                if(ltMask != 0x0000) ExtractOpts8<DIRECTION>(maxMVH, ltMask, QuerySize, (jOptCheckBlock * SIMDMULTIPLE8));
                if(ltMask2 != 0x0000) ExtractOpts8<DIRECTION>(maxMVH2, ltMask2, QuerySize, (jOptCheckBlock + 1) * SIMDMULTIPLE8);
                assert(Opt >= 0 && Opt < 256);
                optSIMD = _mm_set1_epi8(Opt & 0xFF);
            }

            jOptCheckBlock += 2;
        } while(jOptCheckBlock < TargetSizeInSIMDBlocks + 1);
    }

}

template<int VARIANT, bool DIRECTION, bool DETECTREDUNDANT> void Solution9::DoVariant16() {

#ifndef NDEBUG
    //        char tag[128];
    //        sprintf(tag, "DoVariant16 v%d d%d r%d", VARIANT, (int)DIRECTION, (int)DETECTREDUNDANT);
#endif

    assert(TargetSizeInSIMDBlocks >= 1);


    MHVTRealType M = MData;
    MHVTRealType V = VData;
    MHVTRealType H = HData;
    MHVTRealType TargetLookup = TargetLookupData;
    //__m128i* QueueNewMPlusOE = (__m128i*)QueueNewMPlusOEData;

    __m128i optSIMD = _mm_set1_epi16(Opt & 0xFF);

    const __m128i gapExtensionSIMD = _mm_set1_epi16(GapExtension);
    const __m128i gapOpenPlusgapExtensionSIMD = _mm_set1_epi16(GapOpen + GapExtension);
    const __m128i negGapOpenPlusgapExtensionSIMD = _mm_set1_epi16(-(GapOpen + GapExtension));
    const __m128i twiceInfinitesPlusOESIMD = _mm_add_epi16(gapOpenPlusgapExtensionSIMD, TwiceInfinitesSIMD16);
    const __m128i gapExtensionTimes8SIMD = _mm_slli_epi16(gapExtensionSIMD, 3);
    const __m128i gapExtensionTimes4SIMD = _mm_slli_epi16(gapExtensionSIMD, 2);
    const __m128i gapExtensionTimes2SIMD = _mm_slli_epi16(gapExtensionSIMD, 1);



    int jBlockSetup = 0;
    do {
        _mm_store_si128(M[0] + jBlockSetup, _mm_setzero_si128());
        _mm_store_si128(H[0] + jBlockSetup, NegInfiniteSIMD16);
        _mm_store_si128(V[0] + jBlockSetup, NegInfiniteSIMD16);
    } while(++jBlockSetup < (TargetSizeInSIMDBlocks + 2));



    for(int i = 1; i <= QuerySize; ++i) {
        int jSecondBlock = 0;
        do {
            __m128i newM = _mm_load_si128(M[(i - 1) & 0x1] + jSecondBlock);
            __m128i newM2 = _mm_load_si128(M[(i - 1) & 0x1] + jSecondBlock + 1);
            __m128i newMPlusOE = _mm_add_epi16(newM, gapOpenPlusgapExtensionSIMD);
            __m128i newMPlusOE2 = _mm_add_epi16(newM2, gapOpenPlusgapExtensionSIMD);

            __m128i oldV = _mm_load_si128(V[(i - 1) & 0x1] + jSecondBlock);
            __m128i oldV2 = _mm_load_si128(V[(i - 1) & 0x1] + jSecondBlock + 1);
            oldV = _mm_add_epi16(oldV, gapExtensionSIMD);
            oldV2 = _mm_add_epi16(oldV2, gapExtensionSIMD);
            __m128i newV = _mm_max_epi16(oldV, newMPlusOE);
            __m128i newV2 = _mm_max_epi16(oldV2, newMPlusOE2);
            _mm_store_si128(V[i & 0x1] + jSecondBlock, newV);
            _mm_store_si128(V[i & 0x1] + jSecondBlock + 1, newV2);

            if(DETECTREDUNDANT == true) {
                // XXX do only when needed
                _mm_store_si128(H[i & 0x1] + jSecondBlock, _mm_setzero_si128());
                _mm_store_si128(H[i & 0x1] + jSecondBlock + 1, _mm_setzero_si128());
            }

            jSecondBlock += 2;
        } while(jSecondBlock < TargetSizeInSIMDBlocks + 1);
        __m128i fixupV = _mm_load_si128(V[i & 0x1] + 0);
        fixupV = _mm_insert_epi16(fixupV, -INFINITE16, 0);
        _mm_store_si128(V[i & 0x1] + 0, fixupV);



        int queueLength;
        if(DETECTREDUNDANT == true) {
            queueLength = 0;
        }
        int currentQueryL = QueryLookup[i - 1];
        assert((currentQueryL & 0xFFFFFFFC) == 0);
        const __m128i* targetLookupCurrent = TargetLookup[currentQueryL];

        __m128i rollingNewM;
        if(VARIANT == 1 || VARIANT == 3) {
            rollingNewM = _mm_setzero_si128();
        }
        if(VARIANT == 2 || VARIANT == 4) {
            rollingNewM = NegInfiniteSIMDShiftedBig16;
        }



        int jInnerBlock = 0;
        do {
            __m128i thisV = _mm_load_si128(V[(i - 1) & 0x1] + jInnerBlock);
            __m128i thisH = _mm_load_si128(H[(i - 1) & 0x1] + jInnerBlock);
            __m128i thisM = _mm_load_si128(M[(i - 1) & 0x1] + jInnerBlock);
            __m128i mx = _mm_max_epi16(thisV, thisH);
            mx = _mm_max_epi16(mx, thisM);
            __m128i matches = _mm_load_si128(targetLookupCurrent + jInnerBlock);
            __m128i newM = _mm_add_epi16(mx, matches);


            __m128i thisV2 = _mm_load_si128(V[(i - 1) & 0x1] + jInnerBlock + 1);
            __m128i thisH2 = _mm_load_si128(H[(i - 1) & 0x1] + jInnerBlock + 1);
            __m128i thisM2 = _mm_load_si128(M[(i - 1) & 0x1] + jInnerBlock  + 1);                
            __m128i mx2 = _mm_max_epi16(thisV2, thisH2);
            mx2 = _mm_max_epi16(mx2, thisM2);
            __m128i matches2 = _mm_load_si128(targetLookupCurrent + jInnerBlock + 1);
            __m128i newM2 = _mm_add_epi16(mx2, matches2);

            if(VARIANT == 1 || VARIANT == 3) {
                newM = _mm_max_epi16(newM, _mm_setzero_si128());
                newM2 = _mm_max_epi16(newM2, _mm_setzero_si128());
            }

            __m128i toStoreM = _mm_slli_si128(newM, SIMDLITTLESHIFT16);
            toStoreM = _mm_or_si128(toStoreM, rollingNewM);
            rollingNewM = _mm_srli_si128(newM, SIMDBIGSHIFT16);
            _mm_store_si128(M[i & 0x1] + jInnerBlock, toStoreM);
            __m128i toStoreM2 = _mm_slli_si128(newM2, SIMDLITTLESHIFT16);
            toStoreM2 = _mm_or_si128(toStoreM2, rollingNewM);

            __m128i checkM;
            if(DETECTREDUNDANT == true) {
                checkM = _mm_max_epi16(toStoreM, toStoreM2);
                checkM = _mm_cmpgt_epi16(checkM, negGapOpenPlusgapExtensionSIMD);
            }

            rollingNewM = _mm_srli_si128(newM2, SIMDBIGSHIFT16);
            _mm_store_si128(M[i & 0x1] + jInnerBlock +  1, toStoreM2);

            if(DETECTREDUNDANT == true) {
                int mask = _mm_movemask_epi8(checkM);
                //_mm_store_si128(QueueNewMPlusOE + queueLength, toStoreMPlusoe);
                WorkQueue[queueLength] = jInnerBlock;
                queueLength += (mask > 0);
            }
            jInnerBlock += 2;
        } while(jInnerBlock < TargetSizeInSIMDBlocks + 1);

        if(DETECTREDUNDANT == true) {
            WorkQueue[queueLength] = TargetSizeInSIMDBlocks + 1;    // So we can look one past the end below without checking length
        }


        __m128i rollingNewHZeroes;
        if(VARIANT == 1 || VARIANT == 3) {
            rollingNewHZeroes = _mm_insert_epi16(_mm_setzero_si128(), -INFINITE16 + INFINITE16 * 2, 0);
        }
        if(VARIANT == 2 || VARIANT == 4) {
            rollingNewHZeroes = _mm_insert_epi16(_mm_setzero_si128(), GapOpen + (GapExtension * i) + INFINITE16 * 2, 0);
        }

        bool continueThirdBlock;
        int jThirdBlock;
        int q;
        if(DETECTREDUNDANT == true) {
            q = 0;
            jThirdBlock = WorkQueue[q];
            continueThirdBlock = (q < queueLength);
        }
        if(DETECTREDUNDANT == false) {
            jThirdBlock = 0;
            continueThirdBlock = true;
        }
        while(continueThirdBlock) {
            //// Load dependencies, Bias everything 2 * INFINITE16 so shifted in zeroes will always vanish after max()
            // XXX could get rid of the bias for variants 1 and 3

            __m128i newM = _mm_load_si128(M[i & 0x1] + jThirdBlock);
            __m128i newMPlusOE = _mm_add_epi16(newM, twiceInfinitesPlusOESIMD);
            __m128i newM2 = _mm_load_si128(M[i & 0x1] + jThirdBlock + 1);
            __m128i newMPlusOE2 = _mm_add_epi16(newM2, twiceInfinitesPlusOESIMD);    

            int ltMask, ltMask2;
            if(VARIANT == 1 || VARIANT == 2) {
                __m128i optCheck = _mm_cmplt_epi16(newM, optSIMD);
                __m128i optCheck2 = _mm_cmplt_epi16(newM2, optSIMD);
                ltMask = _mm_movemask_epi8(optCheck);
                ltMask2 = _mm_movemask_epi8(optCheck2);
            }

            __m128i rollingNewH = rollingNewHZeroes;

            __m128i oldH = _mm_add_epi16(rollingNewH, gapExtensionSIMD);
            __m128i newH = _mm_max_epi16(oldH, newMPlusOE);

            __m128i newH8Shift = _mm_add_epi16(newH, gapExtensionTimes8SIMD);
            __m128i newH2 = _mm_max_epi16(newMPlusOE2, newH8Shift);


            //// Propogate Hs

            __m128i newH4Shift = _mm_add_epi16(newH, gapExtensionTimes4SIMD);
            __m128i newH4Shift2 = _mm_add_epi16(newH2, gapExtensionTimes4SIMD);
            __m128i inbetweenH4Shift = _mm_srli_si128(newH4Shift, 8);
            newH4Shift = _mm_slli_si128(newH4Shift, 8);
            newH4Shift2 = _mm_slli_si128(newH4Shift2, 8);
            newH = _mm_max_epi16(newH, newH4Shift);
            newH2 = _mm_max_epi16(newH2, newH4Shift2);
            newH2 = _mm_max_epi16(newH2, inbetweenH4Shift);

            __m128i newH2Shift = _mm_add_epi16(newH, gapExtensionTimes2SIMD);
            __m128i newH2Shift2 = _mm_add_epi16(newH2, gapExtensionTimes2SIMD);
            __m128i inbetweenH2Shift = _mm_srli_si128(newH2Shift, 12);
            newH2Shift = _mm_slli_si128(newH2Shift, 4);
            newH2Shift2 = _mm_slli_si128(newH2Shift2, 4);
            newH = _mm_max_epi16(newH, newH2Shift);
            newH2 = _mm_max_epi16(newH2, newH2Shift2);
            newH2 = _mm_max_epi16(newH2, inbetweenH2Shift);

            __m128i newH1Shift = _mm_add_epi16(newH, gapExtensionSIMD);
            __m128i newH1Shift2 = _mm_add_epi16(newH2, gapExtensionSIMD);
            __m128i inbetweenH1Shift = _mm_srli_si128(newH1Shift, 14);
            newH1Shift = _mm_slli_si128(newH1Shift, 2);
            newH1Shift2 = _mm_slli_si128(newH1Shift2, 2);
            newH = _mm_max_epi16(newH, newH1Shift);
            newH2 = _mm_max_epi16(newH2, newH1Shift2);
            newH2 = _mm_max_epi16(newH2, inbetweenH1Shift);


            //// Store new Hs, unbias by 2 * INFINITE16

            __m128i toStoreH = _mm_slli_si128(newH, SIMDLITTLESHIFT16);
            __m128i toStoreH2 = _mm_slli_si128(newH2, SIMDLITTLESHIFT16);
            toStoreH = _mm_or_si128(toStoreH, rollingNewHZeroes);
            rollingNewHZeroes = _mm_srli_si128(newH, SIMDBIGSHIFT16);
            toStoreH2 = _mm_or_si128(toStoreH2, rollingNewHZeroes);
            rollingNewHZeroes = _mm_srli_si128(newH2, SIMDBIGSHIFT16);;

            toStoreH = _mm_sub_epi16(toStoreH, TwiceInfinitesSIMD16);
            toStoreH2 = _mm_sub_epi16(toStoreH2, TwiceInfinitesSIMD16);


            _mm_store_si128(H[i & 0x1] + jThirdBlock, toStoreH);
            _mm_store_si128(H[i & 0x1] + jThirdBlock + 1, toStoreH2);


            if(VARIANT == 1 || VARIANT == 2) {
                if((ltMask & ltMask2) != 0xFFFF) {
                    // Order is important because of DIRECTION
                    if(ltMask != 0xFFFF) ExtractOpts16<DIRECTION>(newM, ltMask, i, (jThirdBlock * SIMDMULTIPLE16));
                    if(ltMask2 != 0xFFFF) ExtractOpts16<DIRECTION>(newM2, ltMask2, i, ((jThirdBlock + 1) * SIMDMULTIPLE16));
                    optSIMD = _mm_set1_epi16(Opt);
                }
            }


            if(DETECTREDUNDANT == true) {
                int nextH = _mm_extract_epi16(rollingNewHZeroes, 0);
                q++;
                int nextJThirdBlock = jThirdBlock + 2;
                jThirdBlock = WorkQueue[q];
                continueThirdBlock = (q < queueLength);
                if(__builtin_expect(nextH > 2 * INFINITE16, 0)) {
                    if(jThirdBlock > nextJThirdBlock) {
                        jThirdBlock = nextJThirdBlock;
                        continueThirdBlock = (jThirdBlock < TargetSizeInSIMDBlocks + 1);
                        q--;
                    }
                }
            }
            if(DETECTREDUNDANT == false) {
                jThirdBlock += 2;
                continueThirdBlock = (jThirdBlock < TargetSizeInSIMDBlocks + 1);
            }
        }

    }

    if(VARIANT == 3 || VARIANT == 4) {
        int jOptCheckBlock = 0;
        do {
            __m128i finalM = _mm_load_si128(M[QuerySize & 0x1] + jOptCheckBlock);
            __m128i finalV = _mm_load_si128(V[QuerySize & 0x1] + jOptCheckBlock);
            __m128i finalH = _mm_load_si128(H[QuerySize & 0x1] + jOptCheckBlock);
            __m128i maxMVH = _mm_max_epi16(_mm_max_epi16(finalM, finalV), finalH);

            __m128i finalM2 = _mm_load_si128(M[QuerySize & 0x1] + jOptCheckBlock + 1);
            __m128i finalV2 = _mm_load_si128(V[QuerySize & 0x1] + jOptCheckBlock + 1);
            __m128i finalH2 = _mm_load_si128(H[QuerySize & 0x1] + jOptCheckBlock + 1);
            __m128i maxMVH2 = _mm_max_epi16(_mm_max_epi16(finalM2, finalV2), finalH2);

            __m128i optCheck = _mm_cmplt_epi16(maxMVH, optSIMD);
            __m128i optCheck2 = _mm_cmplt_epi16(maxMVH2, optSIMD);
            int ltMask = _mm_movemask_epi8(optCheck);
            int ltMask2 = _mm_movemask_epi8(optCheck2);

            if((ltMask & ltMask2) != 0xFFFF) {
                // Order is important because of DIRECTION
                if(ltMask != 0xFFFF) ExtractOpts16<DIRECTION>(maxMVH, ltMask, QuerySize, (jOptCheckBlock * SIMDMULTIPLE16));
                if(ltMask2 != 0xFFFF) ExtractOpts16<DIRECTION>(maxMVH2, ltMask2, QuerySize, (jOptCheckBlock + 1) * SIMDMULTIPLE16);
                optSIMD = _mm_set1_epi16(Opt);
            }

            jOptCheckBlock += 2;
        } while(jOptCheckBlock < TargetSizeInSIMDBlocks + 1);
    }
}

void Solution9::DoSetup8(const string& target, const string& query) {
    assert(MismatchScore <= 0);
    assert(MatchScore >= 0);

    MHVTRealType TargetLookup = TargetLookupData;

    int targetSizeRoundedUpForSIMD = ((1 + TargetSize) + SIMDMULTIPLE8 - 1) & ~(SIMDMULTIPLE8 - 1);
    TargetSizeInSIMDBlocks = targetSizeRoundedUpForSIMD / SIMDMULTIPLE8;

    for(int jBlock = 0; jBlock < TargetSizeInSIMDBlocks + 1; jBlock += 8) {
        _mm_prefetch(TargetLookup[0] + jBlock, _MM_HINT_T0);
        _mm_prefetch(TargetLookup[1] + jBlock, _MM_HINT_T0);
        _mm_prefetch(TargetLookup[2] + jBlock, _MM_HINT_T0);
        _mm_prefetch(TargetLookup[3] + jBlock, _MM_HINT_T0);
    }

    // Mismatch in this table will equal zero, match will be the match score plus the magnitude of the mismatch, so all values are offset by mismtch and postitive

    for(int i = 0; i < QuerySize; i++) QueryLookup[i] = CharMap[int(query[i])];

    for(int j = 0; j < TargetSizeInSIMDBlocks; j++) {
        _mm_store_si128(TargetLookup[0] + j, _mm_setzero_si128());
        _mm_store_si128(TargetLookup[1] + j, _mm_setzero_si128());
        _mm_store_si128(TargetLookup[2] + j, _mm_setzero_si128());
        _mm_store_si128(TargetLookup[3] + j, _mm_setzero_si128());
    }

    for(int j = 0; j < TargetSize; j++) ((unsigned char*)(TargetLookup[CharMap[int(target[j])]]))[j] = (unsigned char)(MatchScore - MismatchScore);


    for(int jBlock = 0; jBlock < TargetSizeInSIMDBlocks + 1; jBlock += 8) {
        _mm_prefetch(((__m128i*)MData[0]) + jBlock, _MM_HINT_T0);
        _mm_prefetch(((__m128i*)HData[0]) + jBlock, _MM_HINT_T0);
        _mm_prefetch(((__m128i*)VData[0]) + jBlock, _MM_HINT_T0);
        _mm_prefetch(((__m128i*)MData[1]) + jBlock, _MM_HINT_T0);
        _mm_prefetch(((__m128i*)HData[1]) + jBlock, _MM_HINT_T0);
        _mm_prefetch(((__m128i*)VData[1]) + jBlock, _MM_HINT_T0);
    }

    Opt = MINIMUMASSUMEDOPT;
    //OptSIMD = _mm_set1_epi8(Opt);
    QueryEnd = TargetEnd = NBest = -1;
}

void Solution9::DoSetup16(const string& target, const string& query) {
    MHVTRealType TargetLookup = TargetLookupData;

    int targetSizeRoundedUpForSIMD = ((1 + TargetSize) + SIMDMULTIPLE16 - 1) & ~(SIMDMULTIPLE16 - 1);
    TargetSizeInSIMDBlocks = targetSizeRoundedUpForSIMD / SIMDMULTIPLE16;

    for(int jBlock = 0; jBlock < TargetSizeInSIMDBlocks + 1; jBlock += 4) {
        _mm_prefetch(TargetLookup[0] + jBlock, _MM_HINT_T0);
        _mm_prefetch(TargetLookup[1] + jBlock, _MM_HINT_T0);
        _mm_prefetch(TargetLookup[2] + jBlock, _MM_HINT_T0);
        _mm_prefetch(TargetLookup[3] + jBlock, _MM_HINT_T0);
    }

    for(int i = 0; i < QuerySize; i++) QueryLookup[i] = CharMap[int(query[i])];
    __m128i targetLookupBuilder = _mm_set1_epi16(MismatchScore);
    for(int j = 0; j < TargetSizeInSIMDBlocks; j++) {
        _mm_store_si128(TargetLookup[0] + j, targetLookupBuilder);
        _mm_store_si128(TargetLookup[1] + j, targetLookupBuilder);
        _mm_store_si128(TargetLookup[2] + j, targetLookupBuilder);
        _mm_store_si128(TargetLookup[3] + j, targetLookupBuilder);
    }

    for(int j = 0; j < TargetSize; j++) ((short int*)(TargetLookup[CharMap[int(target[j])]]))[j] = MatchScore;

    for(int jBlock = 0; jBlock < TargetSizeInSIMDBlocks + 1; jBlock += 4) {
        _mm_prefetch(((__m128i*)MData[0]) + jBlock, _MM_HINT_T0);
        _mm_prefetch(((__m128i*)HData[0]) + jBlock, _MM_HINT_T0);
        _mm_prefetch(((__m128i*)VData[0]) + jBlock, _MM_HINT_T0);
        _mm_prefetch(((__m128i*)MData[1]) + jBlock, _MM_HINT_T0);
        _mm_prefetch(((__m128i*)HData[1]) + jBlock, _MM_HINT_T0);
        _mm_prefetch(((__m128i*)VData[1]) + jBlock, _MM_HINT_T0);
    }

    Opt = MINIMUMASSUMEDOPT;
    //OptSIMD = _mm_set1_epi16(Opt);
    QueryEnd = TargetEnd = NBest = -1;

}

int Solution9::process(string& target, string& query, int queryStartClip, int queryEndClip,
                          int matchScore, int mismatchScore, int gapOpen, int gapExtension, int direction,
                          int *_opt, int *_te, int *_qe, int *_n_best) {

    TargetSize = target.size();
    QuerySize = query.size();
    MatchScore = matchScore;
    MismatchScore = mismatchScore;
    GapOpen = gapOpen;
    GapExtension = gapExtension;
    int variant = GetVariant(queryStartClip, queryEndClip);


    int use8 = ((variant == 1 || variant == 3) & (matchScore == 1));
    int threshold8 = 255 + mismatchScore;

    bool detectRedundant = ((variant == 1) | (variant == 3)) & (matchScore == 1);
    switch(variant + (direction << 2)) {
      case 1 + (0 << 2):    // Variant 1, direction 0
        if(use8) {
            DoSetup8(target, query);
            if(detectRedundant) DoVariant8<1, false, true>(threshold8);
            else DoVariant8<1, false, false>(threshold8);
        }
        if(!use8 || Opt >= threshold8) {
            DoSetup16(target, query);
            if(detectRedundant) DoVariant16<1, false, true>();
            else DoVariant16<1, false, false>();
        }
        break;
      case 2 + (0 << 2):
        DoSetup16(target, query);
        DoVariant16<2, false, false>();
        break;
      case 3 + (0 << 2):
        if(use8) {
            DoSetup8(target, query);
            if(detectRedundant) DoVariant8<3, false, true>(threshold8);
            else DoVariant8<3, false, false>(threshold8);
        }
        if(!use8 || Opt >= threshold8) {
            DoSetup16(target, query);
            if(detectRedundant) DoVariant16<3, false, true>();
            else DoVariant16<3, false, false>();
        }
        break;
      case 4 + (0 << 2):
        DoSetup16(target, query);
        DoVariant16<4, false, false>();
        break;
      case 1 + (1 << 2):    // Variant 1, direction 1
        if(use8) {
            DoSetup8(target, query);
            if(detectRedundant) DoVariant8<1, true, true>(threshold8);
            else DoVariant8<1, true, false>(threshold8);
        }
        if(!use8 || Opt >= threshold8) {
            DoSetup16(target, query);
            if(detectRedundant) DoVariant16<1, true, true>();
            else DoVariant16<1, true, false>();
        }
        break;
      case 2 + (1 << 2):
        DoSetup16(target, query);
        DoVariant16<2, true, false>();
        break;
      case 3 + (1 << 2):
        if(use8) {
            DoSetup8(target, query);
            if(detectRedundant) DoVariant8<3, true, true>(threshold8);
            else DoVariant8<3, true, false>(threshold8);
        }
        if(!use8 || Opt >= threshold8) {
            DoSetup16(target, query);
            if(detectRedundant) DoVariant16<3, true, true>();
            else DoVariant16<3, true, false>();
        }
        break;
      case 4 + (1 << 2):
        DoSetup16(target, query);
        DoVariant16<4, true, false>();
        break;
      default:
        assert(false);
    }

    // "(query_end, target_end) is the location of the cell corresponding to the optimal alignment"
    // or not...
    QueryEnd -= 1;
    TargetEnd -= 1;

    if(Opt >= INFINITE16 || Opt <= -INFINITE16 || Opt <= MINIMUMASSUMEDOPT || (detectRedundant && Opt <= -(gapOpen + gapExtension))) {
        assert(false);
        return Reference.process(target, query, queryStartClip, queryEndClip, matchScore, mismatchScore, gapOpen, gapExtension, direction, _opt, _te, _qe, _n_best);
    }

    /*
    char resultBuffer[128];
    string resultString(resultBuffer, sprintf(resultBuffer, "%i %i %i %i", Opt, QueryEnd, TargetEnd, NBest));
    return resultString;
    */
    
    (*_opt) = Opt;
    (*_te) = TargetEnd;
    (*_qe) = QueryEnd;
    (*_n_best) = NBest;
}

