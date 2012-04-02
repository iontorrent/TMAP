// Coder: Psyho
// Submission: 40
// URL: http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=10597114&subnum=40
#define INLINE   __attribute__ ((always_inline))
#define NOINLINE __attribute__ ((noinline))

#define ALIGNED __attribute__ ((aligned(16)))

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

//#define PRINTF
#define CONVERT
#define FAST_UPDATE
#define USE_HASHING
#define FULL_HEURISTIC
#define HMSIZE 512

/*
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <emmintrin.h>
#include <map>
*/
#include <stdlib.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "tmap_sw.h"
#include "tmap_vsw_definitions.h"
#include "tmap_vsw_s1.h"

//using namespace std;

#define FOR(i,a,b)  for((i)=(a);i<(b);++(i))
#define REP(i,a)    FOR(i,0,a)

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)

#define INF (tmap_vsw16_max_value>>1) //(1<<14)

#define swap(_a, _b) do { \
    __m128i *T = _a; \
    _a = _b; \
    _b = T; \
} while(0)

// TODO
//const int MAX_DIMB = 1024 + 64;

typedef struct {
    int32_t opt;
    int32_t query_end;
    int32_t target_end;
    int32_t n_best;
} result_t;

static INLINE int16_t 
findMax16(const __m128i *mMax) 
{
  __m128i mshift = _mm_srli_si128((*mMax), 2);
  __m128i m = _mm_max_epi16(mshift, (*mMax));
  mshift = _mm_srli_si128(m, 4);
  m = _mm_max_epi16(mshift, m);
  int16_t* p = (int16_t*)&m;
  return max(p[0], p[4]);
}

INLINE int16_t 
findMax16Simple(const __m128i *mMax) 
{
  int16_t* p = (int16_t*)mMax;
  int16_t mx = p[0];
  mx = max(mx, p[1]);
  mx = max(mx, p[2]);
  mx = max(mx, p[3]);
  mx = max(mx, p[4]);
  mx = max(mx, p[5]);
  mx = max(mx, p[6]);
  mx = max(mx, p[7]);
  return mx;    
}

INLINE uint8_t 
findMax8(const __m128i *mMax) 
{
  __m128i mshift = _mm_srli_si128((*mMax), 1);
  __m128i m = _mm_max_epu8(mshift, (*mMax));
  mshift = _mm_srli_si128(m, 2);
  m = _mm_max_epu8(mshift, m);
  uint8_t* p = (uint8_t*)&m;
  uint8_t mx = p[0];
  mx = max(mx, p[4]);
  mx = max(mx, p[8]);
  mx = max(mx, p[12]);
  return mx;    
}

NOINLINE uint64_t 
hashDNA(const uint8_t *s, const int len) 
{
  uint64_t h = 1;
  int32_t j;
  if (len < 16) {
      REP(j, len) h = h * 1337 + s[j];
      return h;
  }

  __m128i mhash = _mm_set1_epi16(1);
  __m128i mmul = _mm_set1_epi16(13);
  __m128i mzero = _mm_setzero_si128();
  __m128i m0, m1a, m1b;
  int i = 0;
  while (i + 16 < len) {
      m0 = _mm_loadu_si128((__m128i*)&s[i]);
      m1a = _mm_unpacklo_epi8(m0, mzero);
      m1b = _mm_unpackhi_epi8(m0, mzero);
      mhash = _mm_mullo_epi16(mhash, mmul);
      mhash = _mm_add_epi16(mhash, m1a);
      mhash = _mm_mullo_epi16(mhash, mmul);
      mhash = _mm_add_epi16(mhash, m1b);
      i += 16;
  }
  m0 = _mm_loadu_si128((__m128i*)&s[len - 16]);
  m1a = _mm_unpacklo_epi8(m0, mzero);
  m1b = _mm_unpackhi_epi8(m0, mzero);
  mhash = _mm_mullo_epi16(mhash, mmul);
  mhash = _mm_add_epi16(mhash, m1a);
  mhash = _mm_mullo_epi16(mhash, mmul);
  mhash = _mm_add_epi16(mhash, m1b);

  uint32_t* p = (uint32_t*)&mhash;
  h = h * 1337 + p[0];
  h = h * 1337 + p[1];
  h = h * 1337 + p[2];
  h = h * 1337 + p[3];
  return h;
}

#define updateResult(_i, _curMax, _type) do { \
    int32_t _j; \
    if (lastMax >= _curMax && lastMax >= opt) { \
        _type* MM = (_type*)M1; \
        if (lastMax > opt) { \
            opt = lastMax; \
            result->n_best = 0; \
            res_min_pos = 1 << 30; \
            res_max_pos = 0; \
        } \
        if (1 != sizeof(_type)) { \
            REP(_j, len) { \
                if (MM[_j] == opt) { \
                    result->n_best++; \
                    int p = (_i << 10) + POS[_j]; \
                    res_min_pos = min(res_min_pos, p); \
                    res_max_pos = max(res_max_pos, p); \
                } \
            } \
        } else { \
            __m128i mopt = _mm_set1_epi8(opt); \
            int _k; \
            for (_k = 0; _k < len; _k += 16) if (_mm_movemask_epi8(_mm_cmpeq_epi8(*(__m128i*)&((uint8_t*)M1)[_k], mopt))) { \
                FOR(_j, _k, _k + 16) { \
                    if (MM[_j] == opt) { \
                        result->n_best++; \
                        int p = (_i << 10) + POS[_j]; \
                        res_min_pos = min(res_min_pos, p); \
                        res_max_pos = max(res_max_pos, p); \
                    } \
                } \
            } \
        } \
    } \
    lastMax = _curMax; \
} while(0)

#define updateResult8(_i, _curMax) updateResult(_i, _curMax, uint8_t)
#define updateResult16(_i, _curMax) updateResult(_i, _curMax, int16_t)

#define updateResultLast(_val, _type) do { \
    int32_t _j; \
    swap(M0, M1); \
    _type* MM = (_type*)M1; \
    REP(_j, INVALID_POS_NO) MM[INVALID_POS[_j]] = _val; \
    REP(_j, len) { \
        if (MM[_j] > opt) { \
            opt = MM[_j]; \
            result->n_best = 1; \
            res_min_pos = res_max_pos = (m << 10) + POS[_j]; \
        } else if (MM[_j] == opt) { \
            result->n_best++; \
            int p = (m << 10) + POS[_j]; \
            res_min_pos = min(res_min_pos, p); \
            res_max_pos = max(res_max_pos, p); \
        } \
    } \
} while(0)

#define updateResultLast8(_val) updateResultLast(_val, uint8_t)
#define updateResultLast16(_val) updateResultLast(_val, int16_t)

#define processFastVariantB16BitA(a, mm, mi, o, e, qec) do { \
    int32_t _i; \
    __m128i mo = _mm_set1_epi16(o + e); \
    __m128i me = _mm_set1_epi16(e); \
    __m128i minf = _mm_set1_epi16(-INF); \
 \
    REP(_i, segNo) M0[_i] = _mm_setzero_si128(); \
    REP(_i, segNo) V[_i] = minf; \
 \
    REP(_i, m) { \
        __m128i mmin;             \
        if (qec) { \
            mmin = _mm_set1_epi16(max(max(lastMax, opt) - (m - _i) * mm, o + e + e * _i)); \
        } else { \
            mmin = _mm_set1_epi16(o + e + e * _i); \
        } \
 \
        __m128i mM = _mm_load_si128(&M0[segNo - 1]); \
        mM = _mm_slli_si128(mM, 2); \
        if (_i) mM = _mm_insert_epi16(mM, o + e * _i    , 0); \
 \
        __m128i mH = _mm_set1_epi16(o + e + e * _i); \
 \
        __m128i mMax; \
        if (qec)  \
          mMax = _mm_setzero_si128(); \
 \
        __m128i *P = XP[a[_i]]; \
 \
        REP(j, segNo) { \
            __m128i mP = _mm_load_si128(&P[j]); \
            __m128i mV = _mm_load_si128(&V[j]); \
 \
            mM = _mm_add_epi16(mM, mP); \
 \
            if (qec) \
              mMax = _mm_max_epi16(mMax, mM); \
 \
            mM = _mm_max_epi16(mM, mV); \
            mM = _mm_max_epi16(mM, mH); \
            mM = _mm_max_epi16(mM, mmin); \
 \
            _mm_store_si128(&M1[j], mM); \
 \
            mM = _mm_add_epi16(mM, mo); \
            mV = _mm_add_epi16(mV, me); \
            mH = _mm_add_epi16(mH, me); \
            mV = _mm_max_epi16(mV, mM); \
            mH = _mm_max_epi16(mH, mM); \
 \
            mM  = _mm_load_si128(&M0[j]); \
 \
            _mm_store_si128(&V[j], mV); \
        } \
 \
        while (1) { \
            mH = _mm_slli_si128(mH, 2); \
            mH = _mm_insert_epi16(mH, -INF, 0); \
 \
            REP(j, segNo) { \
                mM = _mm_load_si128(&M1[j]); \
 \
                __m128i mtmp = _mm_add_epi16(mM, mo); \
                __m128i mcmp = _mm_cmpgt_epi16(mH, mtmp); \
                if (_mm_movemask_epi8(mcmp) == 0) goto outB11; \
 \
                mM = _mm_max_epi16(mM, mH); \
                _mm_store_si128(&M1[j], mM); \
 \
                mH = _mm_add_epi16(mH, me); \
            } \
        } \
outB11: \
 \
        swap(M0, M1); \
        if (qec) { \
            int curMax = findMax16(&mMax); \
            updateResult16(_i, curMax); \
            if (curMax + (m - _i) * mm < opt) goto finish; \
        } \
    } \
 \
    if (qec == 0 || lastMax >= opt) { \
        updateResultLast16(0); \
    } \
 \
} while(0)

#define processFastVariantB16BitB(a, mm, mi, o, e, qec) do { \
    int32_t i, j; \
 \
    __m128i mo = _mm_set1_epi16(o + e); \
    __m128i me = _mm_set1_epi16(e); \
    __m128i minf = _mm_set1_epi16(-INF); \
 \
    REP(i, segNo) M0[i] = _mm_setzero_si128(); \
    REP(i, segNo) V[i] = minf; \
 \
    int midstep = qec ? m : m * 4 / 5; \
    REP(i, midstep) { \
        __m128i mmin;             \
        if (qec) { \
            mmin = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * mm, o + e + e * i)); \
        } else { \
            mmin = _mm_set1_epi16(o + e + e * i); \
        } \
 \
        __m128i mM = _mm_load_si128(&M0[segNo - 1]); \
        mM = _mm_slli_si128(mM, 2); \
        if (i) mM = _mm_insert_epi16(mM, o + e * i    , 0); \
 \
        __m128i mH = _mm_set1_epi16(o + e + e * i); \
 \
        __m128i mMax; \
        if (qec)  \
          mMax = _mm_setzero_si128(); \
 \
        __m128i *P = XP[a[i]]; \
 \
        REP(j, segNo) { \
            __m128i mP = _mm_load_si128(&P[j]); \
            __m128i mV = _mm_load_si128(&V[j]); \
 \
            mM = _mm_add_epi16(mM, mP); \
 \
            if (qec) \
              mMax = _mm_max_epi16(mMax, mM); \
 \
            mM = _mm_max_epi16(mM, mV); \
            mM = _mm_max_epi16(mM, mH); \
            mM = _mm_max_epi16(mM, mmin); \
 \
            _mm_store_si128(&M1[j], mM); \
 \
            mM = _mm_add_epi16(mM, mo); \
            mV = _mm_add_epi16(mV, me); \
            mH = _mm_add_epi16(mH, me); \
            mV = _mm_max_epi16(mV, mM); \
            mH = _mm_max_epi16(mH, mM); \
 \
            mM  = _mm_load_si128(&M0[j]); \
 \
            _mm_store_si128(&V[j], mV); \
        } \
 \
        while (1) { \
            mH = _mm_slli_si128(mH, 2); \
            mH = _mm_insert_epi16(mH, -INF, 0); \
 \
            REP(j, segNo) { \
                mM = _mm_load_si128(&M1[j]); \
 \
                __m128i mtmp = _mm_add_epi16(mM, mo); \
                __m128i mcmp = _mm_cmpgt_epi16(mH, mtmp); \
                if (_mm_movemask_epi8(mcmp) == 0) goto outB12; \
 \
                mM = _mm_max_epi16(mM, mH); \
                _mm_store_si128(&M1[j], mM); \
 \
                mH = _mm_add_epi16(mH, me); \
            } \
        } \
outB12: \
 \
        swap(M0, M1); \
        if (qec) { \
            int curMax = findMax16(&mMax); \
            updateResult16(i, curMax); \
            if (curMax + (m - i) * mm < opt) goto finish; \
        } \
 \
    } \
 \
    if (!qec) { \
        FOR(i, midstep, m) { \
            __m128i mmin;             \
            if (qec) { \
                mmin = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * mm, o + e + e * i)); \
            } else { \
                mmin = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * (mm - e) + o, o + e + e * i)); \
            } \
 \
            __m128i mM = _mm_load_si128(&M0[segNo - 1]); \
            mM = _mm_slli_si128(mM, 2); \
            if (i) mM = _mm_insert_epi16(mM, o + e * i    , 0); \
 \
            __m128i mH = _mm_set1_epi16(o + e + e * i); \
 \
            __m128i mMax; \
            mMax = _mm_setzero_si128(); \
 \
            __m128i *P = XP[a[i]]; \
 \
            REP(j, segNo) { \
                __m128i mP = _mm_load_si128(&P[j]); \
                __m128i mV = _mm_load_si128(&V[j]); \
 \
                mM = _mm_add_epi16(mM, mP); \
 \
                mMax = _mm_max_epi16(mMax, mM); \
 \
                mM = _mm_max_epi16(mM, mV); \
                mM = _mm_max_epi16(mM, mH); \
                mM = _mm_max_epi16(mM, mmin); \
 \
                _mm_store_si128(&M1[j], mM); \
 \
                mM = _mm_add_epi16(mM, mo); \
                mV = _mm_add_epi16(mV, me); \
                mH = _mm_add_epi16(mH, me); \
                mV = _mm_max_epi16(mV, mM); \
                mH = _mm_max_epi16(mH, mM); \
 \
                mM  = _mm_load_si128(&M0[j]); \
 \
                _mm_store_si128(&V[j], mV); \
            } \
 \
            while (1) { \
                mH = _mm_slli_si128(mH, 2); \
                mH = _mm_insert_epi16(mH, -INF, 0); \
 \
                REP(j, segNo) { \
                    mM = _mm_load_si128(&M1[j]); \
 \
                    __m128i mtmp = _mm_add_epi16(mM, mo); \
                    __m128i mcmp = _mm_cmpgt_epi16(mH, mtmp); \
                    if (_mm_movemask_epi8(mcmp) == 0) goto outB22; \
 \
                    mM = _mm_max_epi16(mM, mH); \
                    _mm_store_si128(&M1[j], mM); \
 \
                    mH = _mm_add_epi16(mH, me); \
                } \
            } \
outB22: \
 \
            swap(M0, M1); \
            if (qec) { \
                int curMax = findMax16(&mMax); \
                updateResult16(i, curMax); \
                if (curMax + (m - i) * mm < opt) goto finish; \
            } else { \
                lastMax = findMax16(&mMax); \
            } \
 \
        } \
    } \
 \
    if (qec == 0 || lastMax >= opt) { \
        updateResultLast16(0); \
    } \
 \
} while(0)

#define processFastVariantA16Bit(a, mm, mi, o, e, iter, qec) do { \
    int32_t i, j; \
 \
    __m128i mo = _mm_set1_epi16(o + e); \
    __m128i me = _mm_set1_epi16(e); \
    __m128i mmin = _mm_set1_epi16(tmap_vsw16_min_value); \
    __m128i mmin0 = _mm_setr_epi16(tmap_vsw16_min_value, 0, 0, 0, 0, 0, 0, 0); \
 \
    if (iter == -1) { \
        REP(i, segNo) M0[i] = mmin; \
        REP(i, segNo) V[i] = mmin; \
    } \
 \
    FOR(i, iter + 1, m) { \
        __m128i mH = mmin; \
 \
        __m128i mco; \
        if (qec) { \
            mco = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * mm, tmap_vsw16_min_value)); \
        } \
 \
        __m128i mM = _mm_load_si128(&M0[segNo - 1]); \
        mM = _mm_slli_si128(mM, 2); \
        mM = _mm_or_si128(mM, mmin0); \
 \
        __m128i *P = XP[a[i]]; \
        __m128i mMax; \
        if (qec)  \
          mMax = mmin; \
 \
        REP(j, segNo) { \
            __m128i mV = _mm_load_si128(&V[j]); \
 \
            mM = _mm_adds_epi16(mM, P[j]); \
 \
            if (qec)  \
              mMax = _mm_max_epi16(mMax, mM); \
 \
            mM = _mm_max_epi16(mM, mH); \
            mM = _mm_max_epi16(mM, mV); \
 \
            _mm_store_si128(&M1[j], mM); \
 \
            mM = _mm_adds_epi16(mM, mo); \
            mV = _mm_adds_epi16(mV, me); \
            mV = _mm_max_epi16(mV, mM); \
            mH = _mm_adds_epi16(mH, me); \
            mH = _mm_max_epi16(mH, mM); \
 \
            mM = _mm_load_si128(&M0[j]); \
 \
            _mm_store_si128(&V[j], mV); \
        } \
 \
 \
        while (1) { \
            mH = _mm_slli_si128(mH, 2); \
            mH = _mm_or_si128(mH, mmin0); \
            REP(j, segNo) { \
                mM = _mm_load_si128(&M1[j]); \
 \
                __m128i mtmp = _mm_adds_epi16(mM, mo); \
                if (qec)  \
                  mtmp = _mm_max_epi16(mtmp, mco); \
                __m128i mcmp = _mm_cmpgt_epi16(mH, mtmp); \
                if (_mm_movemask_epi8(mcmp) == 0) goto outA16_##qec; \
 \
                mM = _mm_max_epi16(mM, mH); \
                _mm_store_si128(&M1[j], mM); \
 \
                mH = _mm_adds_epi16(mH, me); \
            } \
        } \
outA16_##qec: \
 \
 \
        swap(M0, M1); \
        if (qec) { \
            int curMax = findMax16Simple(&mMax); \
            updateResult16(i, curMax); \
            if (curMax + (m - i) * mm < opt) { \
                opt -= tmap_vsw16_min_value; \
                goto finish; \
            } \
        } \
 \
    } \
 \
    if (qec == 0 || lastMax >= opt) { \
        updateResultLast16(tmap_vsw16_min_value); \
    } \
 \
    opt -= tmap_vsw16_min_value; \
 \
} while(0) 

#define processFastVariantA8Bit(a, mm, mi, o, e, qec, iter) do { \
    const __m128i mo = _mm_set1_epi8(-(o + e)); \
    const __m128i me = _mm_set1_epi8(-e); \
    const __m128i m1 = _mm_set1_epi8(1); \
    const __m128i mmaxrange = _mm_set1_epi8(255 + mi - mm); \
    const __m128i mmi = _mm_set1_epi8(-mi); \
    const __m128i mzero = _mm_setzero_si128(); \
 \
    REP(i, segNo) M0[i] = mzero; \
    REP(i, segNo) V[i] = mzero; \
 \
    int midstep = qec ? m : m * 4 / 5; \
    REP(i, midstep) { \
        __m128i mco; \
        if (qec) mco = _mm_set1_epi8(max(max(lastMax, opt) - (m - i) * mm, 0)); \
 \
        __m128i mH = mzero; \
 \
        __m128i mM = _mm_load_si128(&M0[segNo - 1]); \
        mM = _mm_slli_si128(mM, 1); \
 \
        __m128i *P = XP[a[i]]; \
        __m128i mMax = mzero; \
 \
        REP(j, segNo) { \
            __m128i mV = _mm_load_si128(&V[j]); \
 \
            mM = _mm_adds_epu8(mM, P[j]); \
            mM = _mm_subs_epu8(mM, mmi); \
 \
            mMax = _mm_max_epu8(mMax, mM); \
 \
            mM = _mm_max_epu8(mM, mH); \
            mM = _mm_max_epu8(mM, mV); \
 \
            _mm_store_si128(&M1[j], mM); \
 \
            mM = _mm_subs_epu8(mM, mo); \
            mV = _mm_subs_epu8(mV, me); \
            mV = _mm_max_epu8(mV, mM); \
            mH = _mm_subs_epu8(mH, me); \
            mH = _mm_max_epu8(mH, mM); \
 \
            mM = _mm_load_si128(&M0[j]); \
 \
            _mm_store_si128(&V[j], mV); \
        } \
 \
 \
 \
        while (1) { \
            mH = _mm_slli_si128(mH, 1); \
            REP(j, segNo) { \
                mM = _mm_load_si128(&M1[j]); \
 \
                __m128i mtmp = _mm_subs_epu8(mM, mo); \
                mtmp = _mm_adds_epu8(mtmp, m1); \
                if (qec) mtmp = _mm_max_epu8(mtmp, mco); \
                mtmp = _mm_max_epu8(mtmp, mH); \
                __m128i mcmp = _mm_cmpeq_epi8(mH, mtmp); \
                if (_mm_movemask_epi8(mcmp) == 0) goto outA81_##qec; \
 \
                mM = _mm_max_epu8(mM, mH); \
                _mm_store_si128(&M1[j], mM); \
 \
                mH = _mm_subs_epu8(mH, me); \
            } \
        } \
outA81_##qec: \
 \
        swap(M0, M1); \
        if (qec) { \
            int curMax = findMax8(&mMax); \
            updateResult8(i, curMax); \
            if (curMax > 255 + mi - mm) { \
                iter = i; \
                goto processFastVariantA8BitBreak; \
            } \
 \
 \
            /* \
               if (curMax > 255 + mi - mm) { \
               int c0 = (m - i) * mm; \
               int c1 = (m - i) * -e + o; \
            c0 = min(24, c0); \
            __m128i mshift = _mm_set1_epi8(c0); \
            REP(j, segNo) { \
            M0[j] = _mm_subs_epu8(M0[j], mshift); \
            V[j] = _mm_subs_epu8(V[j], mshift); \
            } \
            shift += c0; \
            opt -= c0; \
            lastMax -= c0; \
            curMax -= c0; \
            } \
            */ \
 \
            if (curMax + (m - i) * mm < opt) { \
                iter = -1; \
                goto processFastVariantA8BitBreak; \
            } \
        } else { \
            mMax = _mm_max_epu8(mMax, mmaxrange); \
            __m128i mcmp = _mm_cmpgt_epi8(mMax, mmaxrange); \
            if (_mm_movemask_epi8(mcmp)) { \
                iter = i; \
                goto processFastVariantA8BitBreak; \
            } \
        } \
 \
    } \
 \
    if (!qec) { \
        FOR(i, midstep, m) { \
            __m128i mco; \
            if (qec)  \
              mco = _mm_set1_epi8(max(max(lastMax, opt) - (m - i) * mm, 0)); \
            else \
              mco = _mm_set1_epi8(max(max(lastMax, opt) - (m - i) * (mm - e) + o, 0)); \
 \
            __m128i mH = mzero; \
 \
            __m128i mM = _mm_load_si128(&M0[segNo - 1]); \
            mM = _mm_slli_si128(mM, 1); \
 \
            __m128i *P = XP[a[i]]; \
            __m128i mMax = mzero; \
 \
            REP(j, segNo) { \
                __m128i mV = _mm_load_si128(&V[j]); \
 \
                mM = _mm_adds_epu8(mM, P[j]); \
                mM = _mm_subs_epu8(mM, mmi); \
 \
                mMax = _mm_max_epu8(mMax, mM); \
 \
                mM = _mm_max_epu8(mM, mH); \
                mM = _mm_max_epu8(mM, mV); \
 \
                _mm_store_si128(&M1[j], mM); \
 \
                mM = _mm_subs_epu8(mM, mo); \
                mV = _mm_subs_epu8(mV, me); \
                mV = _mm_max_epu8(mV, mM); \
                mH = _mm_subs_epu8(mH, me); \
                mH = _mm_max_epu8(mH, mM); \
 \
                mM = _mm_load_si128(&M0[j]); \
 \
                _mm_store_si128(&V[j], mV); \
            } \
 \
 \
 \
            while (1) { \
                mH = _mm_slli_si128(mH, 1); \
                REP(j, segNo) { \
                    mM = _mm_load_si128(&M1[j]); \
 \
                    __m128i mtmp = _mm_subs_epu8(mM, mo); \
                    mtmp = _mm_adds_epu8(mtmp, m1); \
                    mtmp = _mm_max_epu8(mtmp, mco); \
                    mtmp = _mm_max_epu8(mtmp, mH); \
                    __m128i mcmp = _mm_cmpeq_epi8(mH, mtmp); \
                    if (_mm_movemask_epi8(mcmp) == 0) goto outA82_##qec; \
 \
                    mM = _mm_max_epu8(mM, mH); \
                    _mm_store_si128(&M1[j], mM); \
 \
                    mH = _mm_subs_epu8(mH, me); \
                } \
            } \
outA82_##qec: \
 \
            swap(M0, M1); \
            if (qec) { \
                int curMax = findMax8(&mMax); \
                updateResult8(i, curMax); \
                if (curMax > 255 + mi - mm) { \
                    iter = i; \
                    goto processFastVariantA8BitBreak; \
                } \
 \
                if (curMax + (m - i) * mm < opt) { \
                    iter = -1; \
                    goto processFastVariantA8BitBreak; \
                } \
            } else { \
                lastMax = findMax8(&mMax); \
                if (lastMax > 255 + mi - mm) { \
                    iter = i; \
                    goto processFastVariantA8BitBreak; \
                } \
            } \
 \
        } \
    } \
 \
    if (qec == 0 || lastMax >= opt) { \
        updateResultLast8(0); \
    } \
 \
    iter = -1; \
} while(0)

NOINLINE  void convertTable(__m128i *T0, __m128i *T1, int segNo) {
    int32_t i;
    const int segNoX = segNo / 2;
    const __m128i mlo = _mm_set1_epi16(0x00FF);
    const __m128i mminval = _mm_set1_epi16(tmap_vsw16_min_value);
    __m128i *T2 = &T1[segNoX];
    REP(i, segNoX) {
        __m128i m0 = _mm_load_si128(&T0[i]);
        __m128i m1 = _mm_srli_si128(m0, 1);
        m0 = _mm_and_si128(m0, mlo);
        m0 = _mm_add_epi16(m0, mminval);
        m1 = _mm_and_si128(m1, mlo);
        m1 = _mm_add_epi16(m1, mminval);
        _mm_store_si128(&T1[i], m0);
        _mm_store_si128(&T2[i], m1);
    }
}
#define calcInvalidPos(_value, _type, _size) do { \
    int32_t _i; \
    INVALID_POS_NO = len - n; \
    if (INVALID_POS_NO <= segNo) { \
        REP(_i, INVALID_POS_NO) { \
            INVALID_POS[_i] = len - 1 - _i * _size; \
            ((_type)XP[0])[INVALID_POS[_i]] = _value; \
            ((_type)XP[1])[INVALID_POS[_i]] = _value; \
            ((_type)XP[2])[INVALID_POS[_i]] = _value; \
            ((_type)XP[3])[INVALID_POS[_i]] = _value; \
        } \
    } else { \
        INVALID_POS_NO = 0; \
        REP(_i, len) if (POS[_i] >= n) { \
            INVALID_POS[INVALID_POS_NO] = _i; \
            ((_type)XP[0])[_i] = _value; \
            ((_type)XP[1])[_i] = _value; \
            ((_type)XP[2])[_i] = _value; \
            ((_type)XP[3])[_i] = _value; \
            INVALID_POS_NO++; \
        } \
    } \
} while(0)

#define convert16Bit(buffer, b, qsc, mm, mi) do { \
    segNo *= 2; \
    int len64 = len; \
    int32_t i, c, k; \
 \
    int ptr = 0; \
    __m128i *XM0 = (__m128i*)&buffer[ptr]; ptr += len64; \
    __m128i *XM1 = (__m128i*)&buffer[ptr]; ptr += len64; \
    __m128i *XV = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[0] = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[1] = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[2] = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[3] = (__m128i*)&buffer[ptr]; ptr += len64; \
    POS = &buffer[ptr]; ptr += len64; \
    int16_t* STARG = &buffer[ptr]; \
 \
    if (M0 > M1) swap(XM0, XM1); \
 \
    convertTable(V, XV, segNo); \
    convertTable(M0, XM0, segNo); \
 \
    M0 = XM0; \
    M1 = XM1; \
    V = XV; \
 \
    __m128i mPos = _mm_setr_epi16(0, segNo, segNo * 2, segNo * 3, segNo * 4, segNo * 5, segNo * 6, segNo * 7); \
    __m128i mOne = _mm_set1_epi16(1); \
    REP(k, segNo) { \
        _mm_store_si128((__m128i*)&POS[k * 8], mPos); \
        mPos = _mm_add_epi16(mPos, mOne); \
    } \
 \
    const int16_t SCHAR[4] = {'A', 'C', 'G', 'T'}; \
    __m128i mmul = _mm_set1_epi16(mm - mi); \
    __m128i madd = _mm_set1_epi16(mi); \
    REP(i, len) STARG[i] = b[POS[i]]; \
    REP(c, 4) { \
        __m128i mChar = _mm_set1_epi16(SCHAR[c]); \
        REP(k, segNo) { \
            __m128i mTarg = _mm_load_si128((__m128i*)&STARG[k * 8]); \
            __m128i mnew = _mm_cmpeq_epi16(mTarg, mChar); \
            mnew = _mm_srli_epi16(mnew, 1); \
            mnew = _mm_min_epi16(mnew, mmul); \
            mnew = _mm_add_epi16(mnew, madd); \
            _mm_store_si128(&XP[c][k], mnew); \
        } \
    } \
 \
    calcInvalidPos(mi, int16_t*, 8); \
 \
} while(0)

#define preprocess16Bit(buffer, b, qsc, mm, mi) do { \
    int32_t i, c, k; \
    segNo = (n + 7) / 8; \
    len = segNo * 8; \
    int len64 = len; \
 \
    int ptr = 0; \
    M0 = (__m128i*)&buffer[ptr]; ptr += len64; \
    M1 = (__m128i*)&buffer[ptr]; ptr += len64; \
    V = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[0] = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[1] = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[2] = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[3] = (__m128i*)&buffer[ptr]; ptr += len64; \
    POS = &buffer[ptr]; ptr += len64; \
    int16_t* STARG = &buffer[ptr]; \
 \
      { \
        __m128i mPos = _mm_setr_epi16(0, segNo, segNo * 2, segNo * 3, segNo * 4, segNo * 5, segNo * 6, segNo * 7); \
        __m128i mOne = _mm_set1_epi16(1); \
        REP(k, segNo) { \
            _mm_store_si128((__m128i*)&POS[k * 8], mPos); \
            mPos = _mm_add_epi16(mPos, mOne); \
        } \
 \
        const int16_t SCHAR[4] = {'A', 'C', 'G', 'T'}; \
        __m128i mmul = _mm_set1_epi16(mm - mi); \
        __m128i madd = _mm_set1_epi16(mi); \
        REP(i, len) STARG[i] = b[POS[i]]; \
        REP(c, 4) { \
            __m128i mChar = _mm_set1_epi16(SCHAR[c]); \
            REP(k, segNo) { \
                __m128i mTarg = _mm_load_si128((__m128i*)&STARG[k * 8]); \
                __m128i mnew = _mm_cmpeq_epi16(mTarg, mChar); \
                mnew = _mm_srli_epi16(mnew, 1); \
                mnew = _mm_min_epi16(mnew, mmul); \
                mnew = _mm_add_epi16(mnew, madd); \
                _mm_store_si128(&XP[c][k], mnew); \
            } \
        } \
      } \
 \
    calcInvalidPos(mi, int16_t*, 8); \
} while(0)

#define preprocess8Bit(buffer, b, qsc, mm, mi) do { \
    int32_t i, c, k; \
    segNo = (n + 15) / 16; \
    len = segNo * 16; \
    int len64 = len / 2; \
 \
    int ptr = 0; \
    M0 = (__m128i*)&buffer[ptr]; ptr += len64; \
    M1 = (__m128i*)&buffer[ptr]; ptr += len64; \
    V = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[0] = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[1] = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[2] = (__m128i*)&buffer[ptr]; ptr += len64; \
    XP[3] = (__m128i*)&buffer[ptr]; ptr += len64; \
    POS = &buffer[ptr]; ptr += len64 * 2; \
    int8_t* STARG = (int8_t*)&buffer[ptr]; \
 \
      { \
        __m128i mPos0 = _mm_setr_epi16(0, segNo, segNo * 2, segNo * 3, segNo * 4, segNo * 5, segNo * 6, segNo * 7); \
        __m128i mPos1 = _mm_setr_epi16(segNo * 8, segNo * 9, segNo * 10, segNo * 11, segNo * 12, segNo * 13, segNo * 14, segNo * 15); \
        __m128i mOne = _mm_set1_epi16(1); \
        REP(k, segNo) { \
            _mm_store_si128((__m128i*)&POS[k * 16 + 0], mPos0); \
            _mm_store_si128((__m128i*)&POS[k * 16 + 8], mPos1); \
            mPos0 = _mm_add_epi16(mPos0, mOne); \
            mPos1 = _mm_add_epi16(mPos1, mOne); \
        } \
 \
        const int8_t SCHAR[4] = {'A', 'C', 'G', 'T'}; \
        __m128i mmul = _mm_set1_epi8(mm - mi); \
        REP(i, len) STARG[i] = b[POS[i]]; \
        REP(c, 4) { \
            __m128i mChar = _mm_set1_epi8(SCHAR[c]); \
            REP(k, segNo) { \
                __m128i mTarg = _mm_load_si128((__m128i*)&STARG[k * 16]); \
                __m128i mnew = _mm_cmpeq_epi8(mTarg, mChar); \
                mnew = _mm_min_epu8(mnew, mmul); \
                _mm_store_si128(&XP[c][k], mnew); \
            } \
        } \
      } \
 \
    calcInvalidPos(0, int8_t*, 16); \
} while(0)

static void
process(tmap_vsw_data_s1_t *vsw, const uint8_t *bs, const int32_t n, const uint8_t *as, int32_t m,
        int32_t qsc, int32_t qec,
        int32_t mm, int32_t mi, int32_t o, int32_t e, int32_t dir,
        result_t *result)

{
  const uint8_t *a = NULL, *b = NULL;
  tmap_vsw_s1_qres_t *htdata = vsw->htdata;
  int32_t res_min_pos = 0, res_max_pos = 0, opt = 0, pos;
  int32_t lastMax, segNo, len;
  int16_t *buffer = vsw->buffer;
  __m128i* XP[4];
  __m128i* M0;
  __m128i* M1;
  __m128i* V;
  int16_t* POS;
  int16_t INVALID_POS[16];
  int32_t INVALID_POS_NO;

  int32_t iter = -1;
  int32_t i, j;
  
  a = as;
  b = bs;

#ifdef USE_HASHING
  uint64_t curHashA = hashDNA(a, m);
  uint64_t curHashB = hashDNA(b, n);

  uint64_t hash = curHashA ^ curHashB;
  if (qsc) hash ^= 0x5555555555555555ULL;
  if (qec) hash ^= 0xCCCCCCCCCCCCCCCCULL;
  uint32_t hh = hash >> 8;
  int32_t htpos = hash & (HMSIZE - 1);

  if (htdata[htpos].hash == hh) {
      tmap_vsw_s1_qres_t res = htdata[htpos];
      opt = res.opt;
      result->n_best = res.n_best;
      res_min_pos = res.res_min_pos;
      res_max_pos = res.res_max_pos;
      goto similar;
  }
#endif

#ifdef FULL_HEURISTIC
  if (qsc == 0) {
      int corr = 0;
      int c0 = 0;
      int c1 = 0;
      int offs = n - m + 1;
      //uint32_t h = 0;
      REP(i, offs) {
          REP(j, m) if (a[j] != b[i+j]) {
              goto next;
          }
          if (corr == 0) c0 = i;
          c1 = i;
          corr++;
next: ;
      }
      if (corr) {
          result->n_best = corr;
          opt = mm * m;
          res_min_pos = (m << 10) + c0 + m - 1 - (1 << 10);
          res_max_pos = (m << 10) + c1 + m - 1 - (1 << 10);
          goto similar;
      }
  }
#endif

  //REP(i, m) a[i] = DNA_CONV[a[i]];

  opt = tmap_vsw16_min_value;
  lastMax = tmap_vsw16_min_value;

  if (qsc == 0) goto process16bit; 

  preprocess8Bit(buffer, b, qsc, mm, mi);


  if (qsc) {
      if (qec) {
          processFastVariantA8Bit(a, mm, mi, o, e, 1, iter);
      } else {
          processFastVariantA8Bit(a, mm, mi, o, e, 0, iter);
      }
  }
processFastVariantA8BitBreak:

  if (iter < 0) goto finish;
#ifdef CONVERT
  convert16Bit(buffer, b, qsc, mm, mi);
  opt += tmap_vsw16_min_value;
  lastMax += tmap_vsw16_min_value;
  goto go16bit;
#else
  opt = tmap_vsw16_min_value;
  lastMax = tmap_vsw16_min_value;
  iter = -1;
#endif

process16bit:        
  preprocess16Bit(buffer, b, qsc, mm, mi);

go16bit:


  if (qsc) {
      if (qec) {
          processFastVariantA16Bit(a, mm, mi, o, e, iter, 1);
      } else {
          processFastVariantA16Bit(a, mm, mi, o, e, iter, 0);
      }
  } else {
      if (qec) {
          processFastVariantB16BitA(a, mm, mi, o, e, 1);
      } else {
          processFastVariantB16BitB(a, mm, mi, o, e, 0);
      }
  }

finish:
  res_min_pos -= 1 << 10;
  res_max_pos -= 1 << 10;

#ifdef USE_HASHING
    {
      tmap_vsw_s1_qres_t res = htdata[htpos];
      res.hash = hh;
      res.opt = opt;
      res.n_best = result->n_best;
      res.res_min_pos = res_min_pos;
      res.res_max_pos = res_max_pos;
    }
#endif

similar:
  // store results
  pos = (dir == 0) ? res_min_pos : res_max_pos;
  result->target_end = pos & 1023;
  result->query_end = pos >> 10;
  result->opt = opt;
}

/*
int32_t
tmap_vsw_process_s1(tmap_vsw_data_s1_t *vsw,
                    const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen,
                    int32_t qsc, int32_t qec, tmap_vsw_opt_t *opt,
                    int32_t direction, int32_t score_thr,
                    int32_t *query_end, int32_t *target_end,
                    int32_t *n_best, int32_t *overflow)
                    */

static tmap_vsw_data_s1_t*
tmap_vsw_data_init_s1_helper(const uint8_t *query, int32_t qlen, int32_t tlen, int32_t qsc, int32_t qec, tmap_vsw_opt_t *opt)
{
  tmap_vsw_data_s1_t *vsw = NULL;
  vsw = tmap_calloc(1, sizeof(tmap_vsw_data_s1_t), "vsw");

  vsw->mem_qlen = qlen;
  vsw->mem_tlen = tlen;
  tmap_roundup32(vsw->mem_qlen);
  tmap_roundup32(vsw->mem_tlen);

  vsw->htdata = tmap_calloc(vsw->mem_qlen, sizeof(tmap_vsw_s1_qres_t), "htdata");
  vsw->buffer = tmap_calloc(vsw->mem_tlen * 9, sizeof(int16_t), "buffer"); // NB: need nine of these...

  vsw->query_start_clip = qsc;
  vsw->query_end_clip = qec;
  vsw->opt = opt;
  vsw->max_qlen = INT32_MAX;
  vsw->max_tlen = INT32_MAX;
  return vsw;
}

tmap_vsw_data_s1_t*
tmap_vsw_data_init_s1(const uint8_t *query, int32_t qlen, int32_t query_start_clip, int32_t query_end_clip, tmap_vsw_opt_t *opt)
{
  return tmap_vsw_data_init_s1_helper(query, qlen, qlen, query_start_clip, query_end_clip, opt);
}

tmap_vsw_data_s1_t*
tmap_vsw_data_update_s1(tmap_vsw_data_s1_t *vsw, const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen)
{
  int32_t mem_qlen, mem_tlen, qsc, qec;
  tmap_vsw_opt_t *opt = NULL;

  if(vsw->mem_qlen < qlen || vsw->mem_tlen < tlen) { // update
      // NB: do not re-use memory, since we want to try to get it all in one
      // block
      mem_qlen = (vsw->mem_qlen < qlen) ? qlen : vsw->mem_qlen;
      mem_tlen = (vsw->mem_tlen < tlen) ? tlen : vsw->mem_tlen;
      qsc = vsw->query_start_clip;
      qec = vsw->query_end_clip;
      opt = vsw->opt;
      tmap_vsw_data_destroy_s1(vsw); // destroy
      vsw = tmap_vsw_data_init_s1_helper(query, mem_qlen, mem_tlen, qsc, qec, opt);
  }
  return vsw;
}

void
tmap_vsw_data_destroy_s1(tmap_vsw_data_s1_t *vsw)
{
  if(NULL == vsw) return;
  free(vsw->htdata);
  free(vsw->buffer);
  free(vsw);
}

int32_t
tmap_vsw_process_s1(tmap_vsw_data_s1_t *vsw,
                    const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen,
                    int32_t query_start_clip, int32_t query_end_clip, tmap_vsw_opt_t *opt,
                    int32_t direction, int32_t score_thr,
                    int32_t *query_end, int32_t *target_end,
                    int32_t *n_best, int32_t *overflow)
{
  result_t result;
  process(vsw, target, tlen, query, qlen,
         query_start_clip, query_end_clip,
         opt->score_match, -opt->pen_mm, -opt->pen_gapo, -opt->pen_gape, direction, &result);
  if(result.opt < score_thr) {
      (*query_end) = (*target_end) = -1;
      (*n_best) = 0;
      return INT16_MIN;
  }
  else {
      (*query_end) = result.query_end;
      (*target_end) = result.target_end;
      (*n_best) = result.n_best;
      return result.opt;
  }
}
