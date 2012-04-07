// Coder: Psyho
// Submission: 40
// URL: http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=10597114&subnum=40
#define INLINE   __attribute__ ((always_inline))
#define NOINLINE __attribute__ ((noinline))

#define ALIGNED __attribute__ ((aligned(16)))

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

//#define LIMIT 20000
//#define PRINTF
#define CONVERT
//#define USE_COUNTERS
#define FAST_UPDATE
#define USE_HASHING
#define FULL_HEURISTIC

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
#include <stdint.h>
#include "../../util/tmap_definitions.h"
#include "../../util/tmap_alloc.h"
#include "Solution4.h"

using namespace std;

#define FOR(i,a,b)  for(int i=(a);i<(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ZERO(m)     memset(m,0,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define S           size()
#define LL          int64_t
#define LD          long double
#define MP          make_pair
#define X           first
#define Y           second
#define VC          vector
#define PII         pair <int, int>
#define VI          VC < int >
#define VVI         VC < VI >
#define VD          VC < double >
#define VS          VC < string >
#define DB(a)       cout << #a << ": " << (a) << endl;


#ifdef USE_COUNTERS
#define ADD(x) x++
#else
#define ADD(x) ;
#endif

/*
void print(VI v) {cout << "[";if (v.S) cout << v[0];FOR(i, 1, v.S) cout << ", " << v[i];cout << "]\n";}
void print(VC < LL > v) {cout << "[";if (v.S) cout << v[0];FOR(i, 1, v.S) cout << ", " << v[i];cout << "]\n";}
void print(VD v) {cout << "[";if (v.S) cout << v[0];FOR(i, 1, v.S) cout << ", " << v[i];cout << "]\n";}
void print(VS v) {cout << "[";if (v.S) cout << v[0];FOR(i, 1, v.S) cout << ", " << v[i];cout << "]\n";}
template<class T> string i2s(T x) {ostringstream o; o << x; return o.str(); }
VS splt(string s, char c = ' ') {VS rv; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) rv.PB(s.substr(p, np - p)); p = np + 1;} if (p < s.S) rv.PB(s.substr(p)); return rv;}
*/

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)

const int INF = (1<<14);
const int MIN_VAL = -(1<<15);

INLINE int16_t Solution4::findMax16(const __m128i &mMax) {
    __m128i mshift = _mm_srli_si128(mMax, 2);
    __m128i m = _mm_max_epi16(mshift, mMax);
    mshift = _mm_srli_si128(m, 4);
    m = _mm_max_epi16(mshift, m);
    int16_t* p = (int16_t*)&m;
    return max(p[0], p[4]);
}

INLINE int16_t Solution4::findMax16Simple(const __m128i &mMax) {
    int16_t* p = (int16_t*)&mMax;
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

/*
   INLINE uint8_t findMax8(const __m128i &mMax) {
   __m128i mshift = _mm_srli_si128(mMax, 1);
   __m128i m = _mm_max_epu8(mshift, mMax);
   uint8_t* p = (uint8_t*)&m;
   return max(max(max(p[0], p[2]), max(p[4], p[6])), max(max(p[8], p[10]), max(p[12], p[14])));    
   }
   */

INLINE uint8_t Solution4::findMax8(const __m128i &mMax) {
    __m128i mshift = _mm_srli_si128(mMax, 1);
    __m128i m = _mm_max_epu8(mshift, mMax);
    mshift = _mm_srli_si128(m, 2);
    m = _mm_max_epu8(mshift, m);
    uint8_t* p = (uint8_t*)&m;
    uint8_t mx = p[0];
    mx = max(mx, p[4]);
    mx = max(mx, p[8]);
    mx = max(mx, p[12]);
    return mx;    
}

/*
void show8(const __m128i &m) {
    uint8_t* p = (uint8_t*)&m;
    print(VI(p, p + 16));
}

void show16(const __m128i &m) {
    int16_t* p = (int16_t*)&m;
    print(VI(p, p + 8));
}
*/

NOINLINE uint64_t Solution4::hashDNA(const string &s, const int len) {
    uint64_t h = 1;
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

int Solution4::resize(int a, int b) {
    int i, MAX_DIMA_OLD=MAX_DIMA;
    if(MAX_DIMB < b) {
        MAX_DIMB = b;
        tmap_roundup32(MAX_DIMB);
        BUFFER = (int16_t*)tmap_realloc(BUFFER, (MAX_DIMB + 64) * 9 * sizeof(int16_t), "BUFFER");
    }
    if(MAX_DIMA < a) {
        MAX_DIMA = a;
        tmap_roundup32(MAX_DIMA);
        HTDATA = (qres_t*)tmap_realloc(HTDATA, (MAX_DIMA + 64) * sizeof(qres_t), "HTDATA"); 
        for(i=MAX_DIMA_OLD;i<MAX_DIMA;i++) {
            HTDATA[i].hash = -1;
        }
    }
    return 1;
}

Solution4::Solution4() {
    int i;
    DNA_CONV['A'] = 0;
    DNA_CONV['C'] = 1;
    DNA_CONV['G'] = 2;
    DNA_CONV['T'] = 3;
    // Initial dimentions 
    MAX_DIMA = 512;
    MAX_DIMB = 1024;
    BUFFER = (int16_t*)tmap_malloc((MAX_DIMB + 64) * 9 * sizeof(int16_t), "BUFFER");
    HTDATA = (qres_t*)tmap_malloc((MAX_DIMA + 64) * sizeof(qres_t), "HTDATA"); 
    for(i=0;i<MAX_DIMA;i++) {
        HTDATA[i].hash = -1;
    }
}

Solution4::~Solution4() {
    free(BUFFER);
    free(HTDATA);
}

template <class T, bool BYTE> INLINE void Solution4::updateResult(int i, int curMax) {
    if (lastMax >= curMax && lastMax >= opt) {
        T* MM = (T*)M1;
        if (lastMax > opt) {
            opt = lastMax;
            n_best = 0;
            res_min_pos = 1 << 30;
            res_max_pos = 0;
        }

#ifdef FAST_UPDATE
        if (!BYTE) {
#endif
            REP(j, len) {
                if (MM[j] == opt) {
                    n_best++;
                    int p = (i << 10) + POS[j];
                    res_min_pos = min(res_min_pos, p);
                    res_max_pos = max(res_max_pos, p);
                }
            }
#ifdef FAST_UPDATE
        } else {
            __m128i mopt = _mm_set1_epi8(opt);
            for (int k = 0; k < len; k += 16) if (_mm_movemask_epi8(_mm_cmpeq_epi8(*(__m128i*)&((uint8_t*)M1)[k], mopt))) {
                FOR(j, k, k + 16) {
                    if (MM[j] == opt) {
                        n_best++;
                        int p = (i << 10) + POS[j];
                        res_min_pos = min(res_min_pos, p);
                        res_max_pos = max(res_max_pos, p);
                    }
                }
            }
        }
#endif
    }
    lastMax = curMax;
}

template <class T, int DEFAULT_VALUE> INLINE void Solution4::updateResultLast() {
    swap(M0, M1);
    T* MM = (T*)M1;
    REP(j, INVALID_POS_NO) MM[INVALID_POS[j]] = DEFAULT_VALUE;
    REP(j, len) {
        if (MM[j] > opt) {
            opt = MM[j];
            n_best = 1;
            res_min_pos = res_max_pos = (m << 10) + POS[j];
        } else if (MM[j] == opt) {
            n_best++;
            int p = (m << 10) + POS[j];
            res_min_pos = min(res_min_pos, p);
            res_max_pos = max(res_max_pos, p);
        }
    }
}

template <int qec> NOINLINE void Solution4::processFastVariantB16BitA(string &a, int mm, int mi, int o, int e) {

    __m128i mo = _mm_set1_epi16(o + e);
    __m128i me = _mm_set1_epi16(e);
    __m128i minf = _mm_set1_epi16(-INF);

    REP(i, segNo) M0[i] = _mm_setzero_si128();
    REP(i, segNo) V[i] = minf;

    REP(i, m) {
        __m128i mmin;            
        if (qec) {
            mmin = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * mm, o + e + e * i));
        } else {
            //mmin = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * (mm - e) + o, o + e + e * i));
            mmin = _mm_set1_epi16(o + e + e * i);
        }

        __m128i mM = _mm_load_si128(&M0[segNo - 1]);
        mM = _mm_slli_si128(mM, 2);
        if (i) mM = _mm_insert_epi16(mM, o + e * i    , 0);

        __m128i mH = _mm_set1_epi16(o + e + e * i);

        __m128i mMax;
        if (qec) 
          mMax = _mm_setzero_si128();

        __m128i *P = XP[(int)a[i]];

        REP(j, segNo) {
            __m128i mP = _mm_load_si128(&P[j]);
            __m128i mV = _mm_load_si128(&V[j]);

            mM = _mm_add_epi16(mM, mP);

            if (qec)
              mMax = _mm_max_epi16(mMax, mM);

            mM = _mm_max_epi16(mM, mV);
            mM = _mm_max_epi16(mM, mH);
            mM = _mm_max_epi16(mM, mmin);

            _mm_store_si128(&M1[j], mM);

            mM = _mm_add_epi16(mM, mo);
            mV = _mm_add_epi16(mV, me);
            mH = _mm_add_epi16(mH, me);
            mV = _mm_max_epi16(mV, mM);
            mH = _mm_max_epi16(mH, mM);

            mM  = _mm_load_si128(&M0[j]);

            _mm_store_si128(&V[j], mV);
        }

        while (true) {
            mH = _mm_slli_si128(mH, 2);
            mH = _mm_insert_epi16(mH, -INF, 0);

            REP(j, segNo) {
                mM = _mm_load_si128(&M1[j]);

                __m128i mtmp = _mm_add_epi16(mM, mo);
                __m128i mcmp = _mm_cmpgt_epi16(mH, mtmp);
                if (_mm_movemask_epi8(mcmp) == 0) goto outB11;

                mM = _mm_max_epi16(mM, mH);
                _mm_store_si128(&M1[j], mM);

                mH = _mm_add_epi16(mH, me);
            }
        }
outB11:

        swap(M0, M1);
        if (qec) {
            int curMax = findMax16(mMax);
            updateResult<int16_t, false>(i, curMax);
            if (curMax + (m - i) * mm < opt) return;
        } else {
            //lastMax = findMax16(mMax);
        }

    }

    if (qec == 0 || lastMax >= opt) {
        updateResultLast<int16_t, 0>();
    }

}

template <int qec> NOINLINE void Solution4::processFastVariantB16BitB(string &a, int mm, int mi, int o, int e) {

    __m128i mo = _mm_set1_epi16(o + e);
    __m128i me = _mm_set1_epi16(e);
    __m128i minf = _mm_set1_epi16(-INF);

    REP(i, segNo) M0[i] = _mm_setzero_si128();
    REP(i, segNo) V[i] = minf;

    int midstep = qec ? m : m * 4 / 5;
    REP(i, midstep) {
        __m128i mmin;            
        if (qec) {
            mmin = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * mm, o + e + e * i));
        } else {
            //mmin = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * (mm - e) + o, o + e + e * i));
            mmin = _mm_set1_epi16(o + e + e * i);
        }

        __m128i mM = _mm_load_si128(&M0[segNo - 1]);
        mM = _mm_slli_si128(mM, 2);
        if (i) mM = _mm_insert_epi16(mM, o + e * i    , 0);

        __m128i mH = _mm_set1_epi16(o + e + e * i);

        __m128i mMax;
        if (qec) 
          mMax = _mm_setzero_si128();

        __m128i *P = XP[(int)a[i]];

        REP(j, segNo) {
            __m128i mP = _mm_load_si128(&P[j]);
            __m128i mV = _mm_load_si128(&V[j]);

            mM = _mm_add_epi16(mM, mP);

            if (qec)
              mMax = _mm_max_epi16(mMax, mM);

            mM = _mm_max_epi16(mM, mV);
            mM = _mm_max_epi16(mM, mH);
            mM = _mm_max_epi16(mM, mmin);

            _mm_store_si128(&M1[j], mM);

            mM = _mm_add_epi16(mM, mo);
            mV = _mm_add_epi16(mV, me);
            mH = _mm_add_epi16(mH, me);
            mV = _mm_max_epi16(mV, mM);
            mH = _mm_max_epi16(mH, mM);

            mM  = _mm_load_si128(&M0[j]);

            _mm_store_si128(&V[j], mV);
        }

        while (true) {
            mH = _mm_slli_si128(mH, 2);
            mH = _mm_insert_epi16(mH, -INF, 0);

            REP(j, segNo) {
                mM = _mm_load_si128(&M1[j]);

                __m128i mtmp = _mm_add_epi16(mM, mo);
                __m128i mcmp = _mm_cmpgt_epi16(mH, mtmp);
                if (_mm_movemask_epi8(mcmp) == 0) goto outB12;

                mM = _mm_max_epi16(mM, mH);
                _mm_store_si128(&M1[j], mM);

                mH = _mm_add_epi16(mH, me);
            }
        }
outB12:

        swap(M0, M1);
        if (qec) {
            int curMax = findMax16(mMax);
            updateResult<int16_t, false>(i, curMax);
            if (curMax + (m - i) * mm < opt) return;
        } else {
            //lastMax = findMax16(mMax);
        }

    }

    if (!qec) {
        FOR(i, midstep, m) {
            __m128i mmin;            
            if (qec) {
                mmin = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * mm, o + e + e * i));
            } else {
                mmin = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * (mm - e) + o, o + e + e * i));
                //mmin = _mm_set1_epi16(o + e + e * i);
            }

            __m128i mM = _mm_load_si128(&M0[segNo - 1]);
            mM = _mm_slli_si128(mM, 2);
            if (i) mM = _mm_insert_epi16(mM, o + e * i    , 0);

            __m128i mH = _mm_set1_epi16(o + e + e * i);

            __m128i mMax;
            //if (qec) 
            mMax = _mm_setzero_si128();

            __m128i *P = XP[(int)a[i]];

            REP(j, segNo) {
                __m128i mP = _mm_load_si128(&P[j]);
                __m128i mV = _mm_load_si128(&V[j]);

                mM = _mm_add_epi16(mM, mP);

                //if (qec)
                mMax = _mm_max_epi16(mMax, mM);

                mM = _mm_max_epi16(mM, mV);
                mM = _mm_max_epi16(mM, mH);
                mM = _mm_max_epi16(mM, mmin);

                _mm_store_si128(&M1[j], mM);

                mM = _mm_add_epi16(mM, mo);
                mV = _mm_add_epi16(mV, me);
                mH = _mm_add_epi16(mH, me);
                mV = _mm_max_epi16(mV, mM);
                mH = _mm_max_epi16(mH, mM);

                mM  = _mm_load_si128(&M0[j]);

                _mm_store_si128(&V[j], mV);
            }

            while (true) {
                mH = _mm_slli_si128(mH, 2);
                mH = _mm_insert_epi16(mH, -INF, 0);

                REP(j, segNo) {
                    mM = _mm_load_si128(&M1[j]);

                    __m128i mtmp = _mm_add_epi16(mM, mo);
                    __m128i mcmp = _mm_cmpgt_epi16(mH, mtmp);
                    if (_mm_movemask_epi8(mcmp) == 0) goto outB22;

                    mM = _mm_max_epi16(mM, mH);
                    _mm_store_si128(&M1[j], mM);

                    mH = _mm_add_epi16(mH, me);
                }
            }
outB22:

            swap(M0, M1);
            if (qec) {
                int curMax = findMax16(mMax);
                updateResult<int16_t, false>(i, curMax);
                if (curMax + (m - i) * mm < opt) return;
            } else {
                lastMax = findMax16(mMax);
            }

        }
    }

    if (qec == 0 || lastMax >= opt) {
        updateResultLast<int16_t, 0>();
    }

}

template <int qec> NOINLINE void Solution4::processFastVariantA16Bit(string &a, int mm, int mi, int o, int e, int iter) {
    __m128i mo = _mm_set1_epi16(o + e);
    __m128i me = _mm_set1_epi16(e);
    __m128i mmin = _mm_set1_epi16(MIN_VAL);
    __m128i mmin0 = _mm_setr_epi16(MIN_VAL, 0, 0, 0, 0, 0, 0, 0);

    if (iter == -1) {
        REP(i, segNo) M0[i] = mmin;
        REP(i, segNo) V[i] = mmin;
    }

    FOR(i, iter + 1, m) {
        __m128i mH = mmin;

        __m128i mco;
        if (qec) {
            mco = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * mm, MIN_VAL));
        } else {
            //mco = _mm_set1_epi16(max(max(lastMax, opt) - (m - i) * (mm - e) + o, MIN_VAL));
        }

        __m128i mM = _mm_load_si128(&M0[segNo - 1]);
        mM = _mm_slli_si128(mM, 2);
        mM = _mm_or_si128(mM, mmin0);
        //mM = _mm_insert_epi16(mM, MIN_VAL, 0);

        __m128i *P = XP[(int)a[i]];
        __m128i mMax;
        if (qec) 
          mMax = mmin;

        REP(j, segNo) {
            __m128i mV = _mm_load_si128(&V[j]);

            mM = _mm_adds_epi16(mM, P[j]);

            if (qec) 
              mMax = _mm_max_epi16(mMax, mM);

            mM = _mm_max_epi16(mM, mH);
            mM = _mm_max_epi16(mM, mV);

            _mm_store_si128(&M1[j], mM);

            mM = _mm_adds_epi16(mM, mo);
            mV = _mm_adds_epi16(mV, me);
            mV = _mm_max_epi16(mV, mM);
            mH = _mm_adds_epi16(mH, me);
            mH = _mm_max_epi16(mH, mM);

            mM = _mm_load_si128(&M0[j]);

            _mm_store_si128(&V[j], mV);
            ADD(count0);
        }


        while (true) {
            mH = _mm_slli_si128(mH, 2);
            mH = _mm_or_si128(mH, mmin0);
            //mH = _mm_insert_epi16(mH, MIN_VAL, 0);
            REP(j, segNo) {
                ADD(count1);
                mM = _mm_load_si128(&M1[j]);

                __m128i mtmp = _mm_adds_epi16(mM, mo);
                if (qec) 
                  mtmp = _mm_max_epi16(mtmp, mco);
                __m128i mcmp = _mm_cmpgt_epi16(mH, mtmp);
                if (_mm_movemask_epi8(mcmp) == 0) goto outA16;

                mM = _mm_max_epi16(mM, mH);
                _mm_store_si128(&M1[j], mM);

                mH = _mm_adds_epi16(mH, me);
            }
        }
outA16:


        swap(M0, M1);
        if (qec) {
            int curMax = findMax16Simple(mMax);
            updateResult<int16_t, false>(i, curMax);
            if (curMax + (m - i) * mm < opt) {
                opt -= MIN_VAL;
                return;
            }
        } else {
            //lastMax = findMax16Simple(mMax);
        }

    }

    if (qec == 0 || lastMax >= opt) {
        updateResultLast<int16_t, MIN_VAL>();
    }

    opt -= MIN_VAL;

}

template <int qec> NOINLINE int Solution4::processFastVariantA8Bit(const string &a, const int mm, const int mi, const int o, const int e) {
    const __m128i mo = _mm_set1_epi8(-(o + e));
    const __m128i me = _mm_set1_epi8(-e);
    const __m128i m1 = _mm_set1_epi8(1);
    //const __m128i m128 = _mm_set1_epi8(128);
    const __m128i mmaxrange = _mm_set1_epi8(255 + mi - mm);
    const __m128i mmi = _mm_set1_epi8(-mi);
    const __m128i mzero = _mm_setzero_si128();

    REP(i, segNo) M0[i] = mzero;
    REP(i, segNo) V[i] = mzero;

    int midstep = qec ? m : m * 4 / 5;
    REP(i, midstep) {
        __m128i mco;
        if (qec) mco = _mm_set1_epi8(max(max(lastMax, opt) - (m - i) * mm, 0));

        __m128i mH = mzero;

        __m128i mM = _mm_load_si128(&M0[segNo - 1]);
        mM = _mm_slli_si128(mM, 1);

        __m128i *P = XP[(int)a[i]];
        __m128i mMax = mzero;

        REP(j, segNo) {
            __m128i mV = _mm_load_si128(&V[j]);

            mM = _mm_adds_epu8(mM, P[j]);
            mM = _mm_subs_epu8(mM, mmi);

            mMax = _mm_max_epu8(mMax, mM);

            mM = _mm_max_epu8(mM, mH);
            mM = _mm_max_epu8(mM, mV);

            _mm_store_si128(&M1[j], mM);

            mM = _mm_subs_epu8(mM, mo);
            mV = _mm_subs_epu8(mV, me);
            mV = _mm_max_epu8(mV, mM);
            mH = _mm_subs_epu8(mH, me);
            mH = _mm_max_epu8(mH, mM);

            mM = _mm_load_si128(&M0[j]);

            _mm_store_si128(&V[j], mV);
            ADD(count0);
        }



        while (true) {
            mH = _mm_slli_si128(mH, 1);
            REP(j, segNo) {
                ADD(count1);
                mM = _mm_load_si128(&M1[j]);

                __m128i mtmp = _mm_subs_epu8(mM, mo);
                mtmp = _mm_adds_epu8(mtmp, m1);
                if (qec) mtmp = _mm_max_epu8(mtmp, mco);
                mtmp = _mm_max_epu8(mtmp, mH);
                __m128i mcmp = _mm_cmpeq_epi8(mH, mtmp);
                if (_mm_movemask_epi8(mcmp) == 0) goto outA81;

                mM = _mm_max_epu8(mM, mH);
                _mm_store_si128(&M1[j], mM);

                mH = _mm_subs_epu8(mH, me);
            }
        }
outA81:

        swap(M0, M1);
        if (qec) {
            int curMax = findMax8(mMax);
            updateResult<uint8_t, true>(i, curMax);
            if (curMax > 255 + mi - mm) return i;


            /*
               if (curMax > 255 + mi - mm) {
               int c0 = (m - i) * mm;
               int c1 = (m - i) * -e + o;
            //if (c0 + c1 > curMax - 10) return i;
            c0 = min(24, c0);
            __m128i mshift = _mm_set1_epi8(c0);
            REP(j, segNo) {
            M0[j] = _mm_subs_epu8(M0[j], mshift);
            V[j] = _mm_subs_epu8(V[j], mshift);
            }
            shift += c0;
            opt -= c0;
            lastMax -= c0;
            curMax -= c0;
            }
            */

            if (curMax + (m - i) * mm < opt) {
                return -1;
            }
        } else {
            mMax = _mm_max_epu8(mMax, mmaxrange);
            __m128i mcmp = _mm_cmpgt_epi8(mMax, mmaxrange);
            if (_mm_movemask_epi8(mcmp)) return i;
        }

    }

    if (!qec) {
        FOR(i, midstep, m) {
            __m128i mco;
            if (qec) 
              mco = _mm_set1_epi8(max(max(lastMax, opt) - (m - i) * mm, 0));
            else
              mco = _mm_set1_epi8(max(max(lastMax, opt) - (m - i) * (mm - e) + o, 0));

            __m128i mH = mzero;

            __m128i mM = _mm_load_si128(&M0[segNo - 1]);
            mM = _mm_slli_si128(mM, 1);

            __m128i *P = XP[(int)a[i]];
            __m128i mMax = mzero;

            REP(j, segNo) {
                __m128i mV = _mm_load_si128(&V[j]);

                mM = _mm_adds_epu8(mM, P[j]);
                mM = _mm_subs_epu8(mM, mmi);

                mMax = _mm_max_epu8(mMax, mM);

                mM = _mm_max_epu8(mM, mH);
                mM = _mm_max_epu8(mM, mV);

                _mm_store_si128(&M1[j], mM);

                mM = _mm_subs_epu8(mM, mo);
                mV = _mm_subs_epu8(mV, me);
                mV = _mm_max_epu8(mV, mM);
                mH = _mm_subs_epu8(mH, me);
                mH = _mm_max_epu8(mH, mM);

                mM = _mm_load_si128(&M0[j]);

                _mm_store_si128(&V[j], mV);
                ADD(count0);
            }



            while (true) {
                mH = _mm_slli_si128(mH, 1);
                REP(j, segNo) {
                    ADD(count1);
                    mM = _mm_load_si128(&M1[j]);

                    __m128i mtmp = _mm_subs_epu8(mM, mo);
                    mtmp = _mm_adds_epu8(mtmp, m1);
                    //if (qec) 
                    mtmp = _mm_max_epu8(mtmp, mco);
                    mtmp = _mm_max_epu8(mtmp, mH);
                    __m128i mcmp = _mm_cmpeq_epi8(mH, mtmp);
                    if (_mm_movemask_epi8(mcmp) == 0) goto outA82;

                    mM = _mm_max_epu8(mM, mH);
                    _mm_store_si128(&M1[j], mM);

                    mH = _mm_subs_epu8(mH, me);
                }
            }
outA82:

            swap(M0, M1);
            if (qec) {
                int curMax = findMax8(mMax);
                updateResult<uint8_t, true>(i, curMax);
                if (curMax > 255 + mi - mm) return i;


                /*
                   if (curMax > 255 + mi - mm) {
                   int c0 = (m - i) * mm;
                   int c1 = (m - i) * -e + o;
                //if (c0 + c1 > curMax - 10) return i;
                c0 = min(24, c0);
                __m128i mshift = _mm_set1_epi8(c0);
                REP(j, segNo) {
                M0[j] = _mm_subs_epu8(M0[j], mshift);
                V[j] = _mm_subs_epu8(V[j], mshift);
                }
                shift += c0;
                opt -= c0;
                lastMax -= c0;
                curMax -= c0;
                }
                */

                if (curMax + (m - i) * mm < opt) {
                    return -1;
                }
            } else {
                lastMax = findMax8(mMax);
                if (lastMax > 255 + mi - mm) return i;
            }

        }
    }

    if (qec == 0 || lastMax >= opt) {
        updateResultLast<uint8_t, 0>();
    }

    return -1;
}

NOINLINE  void Solution4::convertTable(__m128i *T0, __m128i *T1) {
    const int segNoX = segNo / 2;
    const __m128i mlo = _mm_set1_epi16(0x00FF);
    const __m128i mminval = _mm_set1_epi16(MIN_VAL);
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

template <class T, int SIZE> NOINLINE  void Solution4::calcInvalidPos(int value) {
    INVALID_POS_NO = len - n;
    if (INVALID_POS_NO <= segNo) {
        REP(i, INVALID_POS_NO) {
            INVALID_POS[i] = len - 1 - i * SIZE;
            ((T)XP[0])[INVALID_POS[i]] = value;
            ((T)XP[1])[INVALID_POS[i]] = value;
            ((T)XP[2])[INVALID_POS[i]] = value;
            ((T)XP[3])[INVALID_POS[i]] = value;
        }
    } else {
        INVALID_POS_NO = 0;
        REP(i, len) if (POS[i] >= n) {
            INVALID_POS[INVALID_POS_NO] = i;
            ((T)XP[0])[i] = value;
            ((T)XP[1])[i] = value;
            ((T)XP[2])[i] = value;
            ((T)XP[3])[i] = value;
            INVALID_POS_NO++;
        }
    }
}

NOINLINE void Solution4::convert16Bit(string &b, int qsc, int mm, int mi) {

    segNo *= 2;
    int len64 = len;

    int ptr = 0;
    __m128i *XM0 = (__m128i*)&BUFFER[ptr]; ptr += len64;
    __m128i *XM1 = (__m128i*)&BUFFER[ptr]; ptr += len64;
    __m128i *XV = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[0] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[1] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[2] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[3] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    POS = &BUFFER[ptr]; ptr += len64;
    int16_t* STARG = &BUFFER[ptr];

    if (M0 > M1) swap(XM0, XM1);

    convertTable(V, XV);
    convertTable(M0, XM0);

    M0 = XM0;
    M1 = XM1;
    V = XV;

    __m128i mPos = _mm_setr_epi16(0, segNo, segNo * 2, segNo * 3, segNo * 4, segNo * 5, segNo * 6, segNo * 7);
    __m128i mOne = _mm_set1_epi16(1);
    REP(k, segNo) {
        _mm_store_si128((__m128i*)&POS[k * 8], mPos);
        mPos = _mm_add_epi16(mPos, mOne);
    }

    const int16_t SCHAR[4] = {'A', 'C', 'G', 'T'};
    __m128i mmul = _mm_set1_epi16(mm - mi);
    __m128i madd = _mm_set1_epi16(mi);
    REP(i, len) STARG[i] = b[POS[i]];
    REP(c, 4) {
        __m128i mChar = _mm_set1_epi16(SCHAR[c]);
        REP(k, segNo) {
            __m128i mTarg = _mm_load_si128((__m128i*)&STARG[k * 8]);
            __m128i mnew = _mm_cmpeq_epi16(mTarg, mChar);
            mnew = _mm_srli_epi16(mnew, 1);
            mnew = _mm_min_epi16(mnew, mmul);
            mnew = _mm_add_epi16(mnew, madd);
            _mm_store_si128(&XP[c][k], mnew);
        }
    }

    calcInvalidPos<int16_t*, 8>(mi);

}

NOINLINE void Solution4::preprocess16Bit(string &b, int qsc, int mm, int mi) {
    segNo = (n + 7) / 8;
    len = segNo * 8;
    int len64 = len;

    int ptr = 0;
    M0 = (__m128i*)&BUFFER[ptr]; ptr += len64;
    M1 = (__m128i*)&BUFFER[ptr]; ptr += len64;
    V = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[0] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[1] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[2] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[3] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    POS = &BUFFER[ptr]; ptr += len64;
    int16_t* STARG = &BUFFER[ptr];

      {
        __m128i mPos = _mm_setr_epi16(0, segNo, segNo * 2, segNo * 3, segNo * 4, segNo * 5, segNo * 6, segNo * 7);
        __m128i mOne = _mm_set1_epi16(1);
        REP(k, segNo) {
            _mm_store_si128((__m128i*)&POS[k * 8], mPos);
            mPos = _mm_add_epi16(mPos, mOne);
        }

        const int16_t SCHAR[4] = {'A', 'C', 'G', 'T'};
        __m128i mmul = _mm_set1_epi16(mm - mi);
        __m128i madd = _mm_set1_epi16(mi);
        REP(i, len) {
            if(n <= POS[i]) STARG[i] = 0;
            else STARG[i] = b[POS[i]];
        }
        REP(c, 4) {
            __m128i mChar = _mm_set1_epi16(SCHAR[c]);
            REP(k, segNo) {
                __m128i mTarg = _mm_load_si128((__m128i*)&STARG[k * 8]);
                __m128i mnew = _mm_cmpeq_epi16(mTarg, mChar);
                mnew = _mm_srli_epi16(mnew, 1);
                mnew = _mm_min_epi16(mnew, mmul);
                mnew = _mm_add_epi16(mnew, madd);
                _mm_store_si128(&XP[c][k], mnew);
            }
        }
      }

    calcInvalidPos<int16_t*, 8>(mi);        
}

void Solution4::preprocess8Bit(string &b, int qsc, int mm, int mi) {
    segNo = (n + 15) / 16;
    len = segNo * 16;
    int len64 = len / 2;

    int ptr = 0;
    M0 = (__m128i*)&BUFFER[ptr]; ptr += len64;
    M1 = (__m128i*)&BUFFER[ptr]; ptr += len64;
    V = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[0] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[1] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[2] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    XP[3] = (__m128i*)&BUFFER[ptr]; ptr += len64;
    POS = &BUFFER[ptr]; ptr += len64 * 2;
    int8_t* STARG = (int8_t*)&BUFFER[ptr];

      {
        __m128i mPos0 = _mm_setr_epi16(0, segNo, segNo * 2, segNo * 3, segNo * 4, segNo * 5, segNo * 6, segNo * 7);
        __m128i mPos1 = _mm_setr_epi16(segNo * 8, segNo * 9, segNo * 10, segNo * 11, segNo * 12, segNo * 13, segNo * 14, segNo * 15);
        __m128i mOne = _mm_set1_epi16(1);
        REP(k, segNo) {
            _mm_store_si128((__m128i*)&POS[k * 16 + 0], mPos0);
            _mm_store_si128((__m128i*)&POS[k * 16 + 8], mPos1);
            mPos0 = _mm_add_epi16(mPos0, mOne);
            mPos1 = _mm_add_epi16(mPos1, mOne);
        }

        const int8_t SCHAR[4] = {'A', 'C', 'G', 'T'};
        __m128i mmul = _mm_set1_epi8(mm - mi);
        REP(i, len) STARG[i] = b[POS[i]];
        REP(c, 4) {
            __m128i mChar = _mm_set1_epi8(SCHAR[c]);
            REP(k, segNo) {
                __m128i mTarg = _mm_load_si128((__m128i*)&STARG[k * 16]);
                __m128i mnew = _mm_cmpeq_epi8(mTarg, mChar);
                mnew = _mm_min_epu8(mnew, mmul);
                _mm_store_si128(&XP[c][k], mnew);
            }
        }
      }

    calcInvalidPos<int8_t*, 16>(0);

}

NOINLINE int Solution4::process(string &b, string &a, int qsc, int qec, 
                                int mm, int mi, int o, int e, int dir,
                                int *_opt, int *_te, int *_qe, int *_n_best) {
    n = b.S;
    m = a.S;

    resize(m, n);

    opt = n_best = res_min_pos = res_max_pos = 0;


    int iter = -1;

#ifdef USE_HASHING
    uint64_t curHashA = hashDNA(a, m);
    uint64_t curHashB = hashDNA(b, n);

    uint64_t hash = curHashA ^ curHashB;
    if (qsc) hash ^= 0x5555555555555555ULL;
    if (qec) hash ^= 0xCCCCCCCCCCCCCCCCULL;
    uint32_t hh = hash >> 8;
    int htpos = hash & (MAX_DIMA - 1); // MAX_DIMA must be a power of two

    if (HTDATA[htpos].hash == hh) {
        qres_t &res = HTDATA[htpos];
        opt = res.opt;
        n_best = res.n_best;
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
            n_best = corr;
            opt = mm * m;
            res_min_pos = (m << 10) + c0 + m - 1 - (1 << 10);
            res_max_pos = (m << 10) + c1 + m - 1 - (1 << 10);
            goto similar;
        }
    }
#endif

    REP(i, m) a[i] = DNA_CONV[(int)a[i]];

    opt = MIN_VAL;
    lastMax = MIN_VAL;

    if (qsc == 0) goto process16bit; 

//process8bit:        


    preprocess8Bit(b, qsc, mm, mi);


    if (qsc) {
        if (qec) {
            iter = processFastVariantA8Bit<1>(a, mm, mi, o, e);
        } else {
            iter = processFastVariantA8Bit<0>(a, mm, mi, o, e);
        }
    } else {
    }

    if (iter < 0) goto finish;
#ifdef CONVERT
    convert16Bit(b, qsc, mm, mi);
    opt += MIN_VAL;
    lastMax += MIN_VAL;
    goto go16bit;
#else
    opt = MIN_VAL;
    lastMax = MIN_VAL;
    iter = -1;
#endif

process16bit:        
    preprocess16Bit(b, qsc, mm, mi);

go16bit:


    if (qsc) {
        if (qec) {
            processFastVariantA16Bit<1>(a, mm, mi, o, e, iter);
        } else {
            processFastVariantA16Bit<0>(a, mm, mi, o, e, iter);
        }
    } else {
        if (qec) {
            processFastVariantB16BitA<1>(a, mm, mi, o, e);
        } else {
            processFastVariantB16BitB<0>(a, mm, mi, o, e);
        }
    }

finish:
    res_min_pos -= 1 << 10;
    res_max_pos -= 1 << 10;

#ifdef USE_HASHING
      {
        qres_t &res = HTDATA[htpos];
        res.hash = hh;
        res.opt = opt;
        res.n_best = n_best;
        res.res_min_pos = res_min_pos;
        res.res_max_pos = res_max_pos;
      }
#endif

similar:
    int pos = dir == 0 ? res_min_pos : res_max_pos;

    int target_end = pos & 1023;
    int query_end = pos >> 10;

    (*_opt) = opt;
    (*_te) = target_end;
    (*_qe) = query_end;
    (*_n_best) = n_best;

    return opt;
}
