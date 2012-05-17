// Coder: venco
// Submission: 13 
// URL: http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=274023&subnum=13

#define SUBMISSION 12/66+
#define SCORE /167878+
//#define RESULT_STATS
//#define CLOCK
//#define SUPER_CHECK

#include <string>
#include <vector>
#include <stdexcept>
#include <map>
#include <list>
#include <set>
#include <queue>
#include <bitset>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <complex>
#include <iostream>
#include <iomanip>
#include <sstream>
//#include <strstream>
#include <sys/time.h>
#include <ext/hash_set>
#include <ext/hash_map>
#include <ext/pool_allocator.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <stdint.h>
#include <limits>
#include "Solution8.h"

#define _mm_loadu_si128 _mm_lddqu_si128
//#define volatile

using namespace std;
/*
typedef long long i64;
typedef unsigned long long u64;
typedef unsigned u32;
typedef int i32;
typedef unsigned short u16;
typedef short i16;
typedef unsigned char u8;
typedef signed char i8;
*/
typedef int64_t i64;
typedef uint64_t u64;
typedef int32_t i32;
typedef uint32_t u32;
typedef int16_t i16;
typedef uint16_t u16;
typedef int8_t i8;
typedef uint8_t u8;

static inline ostream& operator<<(ostream& out, const __m128i& v)
{
  out << "( " << _mm_extract_epi16(v, 0)
    << ", " << _mm_extract_epi16(v, 1)
    << ", " << _mm_extract_epi16(v, 2)
    << ", " << _mm_extract_epi16(v, 3)
    << ", " << _mm_extract_epi16(v, 4)
    << ", " << _mm_extract_epi16(v, 5)
    << ", " << _mm_extract_epi16(v, 6)
    << ", " << _mm_extract_epi16(v, 7)
    << " )";
  return out;
}

////////////////////////// macros and typedefs ///////////////////////////////

#define NAME2x(a,b) a##b
#define NAME2(a,b) NAME2x(a,b)
#define ALL(C) (C).begin(), (C).end()
#define forIter(I,C) for(typeof((C).end()) I=(C).begin(); I!=(C).end(); ++I)
#define forIterE(I,C) for(typeof((C).end()) I=(C).begin(), E=(C).end(); I!=E; ++I)
#define forNF(I,F,C) for(unsigned I=(F); I<(C); ++I)
#define forN(I,C) forNF(I,0,C)
#define forEach(I,C) for(unsigned I=0; I<((C).size()); ++I)
#define NOP() do{}while(0)
#define EATSEMICOLON extern int eat_semicolon_variable
#define BRA { EATSEMICOLON
#define KET } EATSEMICOLON
#define ALIGN16 __attribute__((aligned(16)))
#ifdef HOME_RUN
# define CLOCK_FREQ 2400000000ULL
# define CLOCK_DEC 63
#else
# define CLOCK_FREQ 3600000000ULL
# define CLOCK_DEC 99
# define NDEBUG
# ifdef assert
#  undef assert
# endif
# define assert(v) NOP()
# define TR(v) NOP()
#endif

namespace SNS BRA;

#define expect(a, b) (b<0?a:__builtin_expect(a,b))

static inline int max(int a, int b) { return a > b? a: b; }

// query
const char* a;
size_t m;
// target
const char* b;
size_t n;
// scores
int ma, mi, go, ge;
__m128i m_e, m_oe, m_mi, m_ma_over_mi, m_o, m_ma, m_one;
// result
int opt;
size_t n_best, best_i, best_j;

const i16 INF = 0x7f00;

#define TMP_ALIGN __attribute__((aligned(16)))
const size_t TS = 1024+64+64;
const size_t TO = TS/8;

i16 tmp[11][TS] TMP_ALIGN;
#define bv tmp[1]
#define bn tmp[2]
#define br tmp[3]
#define Hr tmp[4]
#define Vr tmp[5]
#define R0r tmp[6]
#define R1r tmp[7]
#define bf tmp[8]
#define ar tmp[10]
#define Rtr(r) ((r)&1?R1r:R0r)

char tmp_offset[1][16] TMP_ALIGN;

char tmp_align_64[64]  __attribute__((aligned(64)));

#ifdef CLOCK
#else
struct ClockStat {
};
# define START_CLOCK() NOP()
# define END_CLOCK(c,m,n) NOP()
# define INIT_CLOCK(c) NOP()
# define PRINT_CLOCK(c) NOP()
# define LOCAL_CLOCK(c,n,m) NOP()
#endif

ClockStat clock[2][2][2];
ClockStat clock_any, clock_full;

#ifdef SUPER_CHECK
namespace SUPER BRA;
int M[TS][TS];
int H[TS][TS];
int V[TS][TS];
int R[TS][TS];
void process0(int n, int m, int qsc, int ,
              int mm, int mi, int o, int e)
{
  if (qsc) {
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
  }
  else {                                      
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
  for ( int i = 0; i <= m; ++i ) for ( int j = 0; j <= n; ++j )
    R[i][j] = max(M[i][j], max(H[i][j], V[i][j]));
}
void check_row(int r, int n, int m, int qsc, int qec,
               const int arr[TS][TS],
               const i16* row)
{
  for ( int j = arr == H; j <= min(r, n); ++j ) {
      int i = r-j;
      if ( i > m ) continue;
      if ( arr[i][j] != row[j] ) {
          const char* name = "?";
          if ( arr == M ) name = "M";
          if ( arr == H ) name = "H";
          if ( arr == V ) name = "V";
          if ( arr == R ) name = "R";
          fTR(qsc|qec|ma|mi|go|ge);
          fTR(name|r|n|m|i|j|arr[i][j]|row[j]);
          forN ( i, 10 ) fTR(i|A(H[i],10));
          forN ( i, 10 ) fTR(i|A(V[i],10));
          forN ( i, 10 ) fTR(i|A(M[i],10));
          throw runtime_error("xx");
      }
  }
}
KET;
# define CHECK_ROW(r,arr,row) SUPER::check_row(r,n,m,qsc,qec,SUPER::arr,row-1)
#else
# define CHECK_ROW(r,arr,row) NOP()
#endif

template<int qsc, int qec, int dir>
void process(size_t n, size_t m)
{
#ifdef SUPER_CHECK
  SUPER::process0(n, m, qsc, qec, ma, mi, go, ge);
#endif
  START_CLOCK();
  int o = go, e = ge, oe = o+e;

  if ( 0 ) {
      __m128i minf = _mm_set1_epi16(-INF);
      if ( qec ) {
          for ( size_t j = 0; j < n; j += 8 ) {
              _mm_store_si128((__m128i*)(bv+j), minf);
              _mm_store_si128((__m128i*)(bn+j), minf);
              _mm_store_si128((__m128i*)(br+j), minf);
          }
      }
  }

#define STEPqec0()                                                      \
    {                                                                   \
      __m128i M = _mm_loadu_si128((__m128i*)(ar-(r-2)+j));            \
      __m128i wb = *((volatile __m128i*)(bf+j));                      \
      __m128i m_ma_over_mi = *((volatile __m128i*)&SNS::m_ma_over_mi); \
      __m128i m_mi = *((volatile __m128i*)&SNS::m_mi);                \
      __m128i Rp = _mm_loadu_si128((__m128i*)(Rr+j-1));               \
      __m128i H = _mm_loadu_si128((__m128i*)(Hr+j-1));                \
      M = _mm_cmpeq_epi16(M, wb);                                     \
      M = _mm_and_si128(M, m_ma_over_mi);                             \
      M = _mm_add_epi16(M, m_mi);                                     \
      M = _mm_add_epi16(M, Rp);                                       \
      if ( qsc ) M = _mm_max_epi16(M, _mm_setzero_si128());           \
      __m128i R = *((volatile __m128i*)(Vr+j));                       \
      __m128i V = *((volatile __m128i*)(Vr+j));                       \
      __m128i m_oe = *((volatile __m128i*)&SNS::m_oe);                \
      __m128i m_e = *((volatile __m128i*)&SNS::m_e);                  \
      R = _mm_max_epi16(R, H);                                        \
      R = _mm_max_epi16(R, M);                                        \
      _mm_store_si128((__m128i*)(Rr+j), R);                           \
      M = _mm_add_epi16(M, m_oe);                                     \
      V = _mm_add_epi16(V, m_e);                                      \
      H = _mm_add_epi16(H, m_e);                                      \
      V = _mm_max_epi16(V, M);                                        \
      H = _mm_max_epi16(H, M);                                        \
      _mm_store_si128((__m128i*)(Vr+j), V);                           \
      _mm_store_si128((__m128i*)(Hr+j), H);                           \
    }
#define STEPqec0p()                                                     \
    {                                                                   \
      __m128i M = _mm_loadu_si128(arp);                               \
      __m128i wb = *((volatile __m128i*)(Hrp+TO*4));                  \
      __m128i m_ma_over_mi = *((volatile __m128i*)&SNS::m_ma_over_mi); \
      __m128i m_mi = *((volatile __m128i*)&SNS::m_mi);                \
      __m128i Rp = _mm_loadu_si128((__m128i*)((i16*)Rrp-1));          \
      __m128i H = _mm_loadu_si128((__m128i*)((i16*)Hrp-1));           \
      M = _mm_cmpeq_epi16(M, wb);                                     \
      M = _mm_and_si128(M, m_ma_over_mi);                             \
      M = _mm_add_epi16(M, m_mi);                                     \
      M = _mm_add_epi16(M, Rp);                                       \
      if ( qsc ) M = _mm_max_epi16(M, _mm_setzero_si128());           \
      __m128i R = *((volatile __m128i*)(Hrp+TO));                     \
      __m128i V = *((volatile __m128i*)(Hrp+TO));                     \
      __m128i m_oe = *((volatile __m128i*)&SNS::m_oe);                \
      __m128i m_e = *((volatile __m128i*)&SNS::m_e);                  \
      R = _mm_max_epi16(R, H);                                        \
      R = _mm_max_epi16(R, M);                                        \
      _mm_store_si128(Rrp, R);                                        \
      M = _mm_add_epi16(M, m_oe);                                     \
      V = _mm_add_epi16(V, m_e);                                      \
      H = _mm_add_epi16(H, m_e);                                      \
      V = _mm_max_epi16(V, M);                                        \
      H = _mm_max_epi16(H, M);                                        \
      _mm_store_si128(Hrp+TO, V);                                     \
      _mm_store_si128(Hrp, H);                                        \
    }
#define STEPqec0p2()                                                    \
    {                                                                   \
      __m128i M1 = _mm_loadu_si128(arp-1);                            \
      __m128i wb1 = *((volatile __m128i*)(Hrp+TO*4-1));               \
      __m128i M2 = _mm_loadu_si128(arp  );                            \
      __m128i wb2 = *((volatile __m128i*)(Hrp+TO*4  ));               \
      __m128i m_ma_over_mi = *((volatile __m128i*)&SNS::m_ma_over_mi); \
      __m128i m_mi = *((volatile __m128i*)&SNS::m_mi);                \
      M1 = _mm_cmpeq_epi16(M1, wb1);                                  \
      M2 = _mm_cmpeq_epi16(M2, wb2);                                  \
      __m128i Rp1 = _mm_loadu_si128((__m128i*)((i16*)Rrp-9));         \
      __m128i Rp2 = _mm_loadu_si128((__m128i*)((i16*)Rrp-1));         \
      M1 = _mm_and_si128(M1, m_ma_over_mi);                           \
      M2 = _mm_and_si128(M2, m_ma_over_mi);                           \
      M1 = _mm_add_epi16(M1, m_mi);                                   \
      M2 = _mm_add_epi16(M2, m_mi);                                   \
      M1 = _mm_add_epi16(M1, Rp1);                                    \
      M2 = _mm_add_epi16(M2, Rp2);                                    \
      __m128i H1 = _mm_loadu_si128((__m128i*)((i16*)Hrp-9));          \
      __m128i H2 = _mm_loadu_si128((__m128i*)((i16*)Hrp-1));          \
      if ( qsc ) {                                                    \
          __m128i zero = _mm_setzero_si128();                         \
          M1 = _mm_max_epi16(M1, zero);                               \
          M2 = _mm_max_epi16(M2, zero);                               \
      }                                                               \
      __m128i R1 = *((volatile __m128i*)(Hrp+TO-1));                  \
      __m128i R2 = *((volatile __m128i*)(Hrp+TO  ));                  \
      R1 = _mm_max_epi16(R1, H1);                                     \
      R2 = _mm_max_epi16(R2, H2);                                     \
      R1 = _mm_max_epi16(R1, M1);                                     \
      R2 = _mm_max_epi16(R2, M2);                                     \
      _mm_store_si128(Rrp-1, R1);                                     \
      _mm_store_si128(Rrp  , R2);                                     \
      __m128i m_oe = *((volatile __m128i*)&SNS::m_oe);                \
      __m128i m_e = *((volatile __m128i*)&SNS::m_e);                  \
      __m128i V1 = *((volatile __m128i*)(Hrp+TO-1));                  \
      __m128i V2 = *((volatile __m128i*)(Hrp+TO  ));                  \
      M1 = _mm_add_epi16(M1, m_oe);                                   \
      M2 = _mm_add_epi16(M2, m_oe);                                   \
      H1 = _mm_add_epi16(H1, m_e);                                    \
      H2 = _mm_add_epi16(H2, m_e);                                    \
      V1 = _mm_add_epi16(V1, m_e);                                    \
      V2 = _mm_add_epi16(V2, m_e);                                    \
      H1 = _mm_max_epi16(H1, M1);                                     \
      H2 = _mm_max_epi16(H2, M2);                                     \
      V1 = _mm_max_epi16(V1, M1);                                     \
      V2 = _mm_max_epi16(V2, M2);                                     \
      _mm_store_si128(Hrp-1, H1);                                     \
      _mm_store_si128(Hrp  , H2);                                     \
      _mm_store_si128(Hrp+TO-1, V1);                                  \
      _mm_store_si128(Hrp+TO  , V2);                                  \
    }
#define STEPqec1p()                                                     \
    {                                                                   \
      __m128i M = _mm_loadu_si128(arp);                               \
      __m128i wb = *((volatile __m128i*)(Hrp+TO*4));                  \
      __m128i m_ma_over_mi = *((volatile __m128i*)&SNS::m_ma_over_mi); \
      __m128i m_mi = *((volatile __m128i*)&SNS::m_mi);                \
      __m128i Rp = _mm_loadu_si128((__m128i*)((i16*)Rrp-1));          \
      __m128i H = _mm_loadu_si128((__m128i*)((i16*)Hrp-1));           \
      M = _mm_cmpeq_epi16(M, wb);                                     \
      M = _mm_and_si128(M, m_ma_over_mi);                             \
      M = _mm_add_epi16(M, m_mi);                                     \
      M = _mm_add_epi16(M, Rp);                                       \
      if ( qsc ) M = _mm_max_epi16(M, _mm_setzero_si128());           \
      __m128i R = *((volatile __m128i*)(Hrp+TO));                     \
      __m128i V = *((volatile __m128i*)(Hrp+TO));                     \
      __m128i m_oe = *((volatile __m128i*)&SNS::m_oe);                \
      __m128i m_e = *((volatile __m128i*)&SNS::m_e);                  \
      R = _mm_max_epi16(R, H);                                        \
      R = _mm_max_epi16(R, M);                                        \
      _mm_store_si128(Rrp, R);                                        \
      M = _mm_add_epi16(M, m_oe);                                     \
      V = _mm_add_epi16(V, m_e);                                      \
      H = _mm_add_epi16(H, m_e);                                      \
      V = _mm_max_epi16(V, M);                                        \
      H = _mm_max_epi16(H, M);                                        \
      _mm_store_si128(Hrp+TO, V);                                     \
      _mm_store_si128(Hrp, H);                                        \
      __m128i n = R;                                                  \
      __m128i vc = R;                                                 \
      __m128i v = _mm_load_si128(bvp);                                \
      __m128i rp = _mm_load_si128(bvp+TO*2);                          \
      __m128i np = _mm_load_si128(bvp+TO);                            \
      if ( dir ) {                                                    \
          n = _mm_cmpgt_epi16(n, v);                                  \
          v = _mm_max_epi16(v, vc);                                   \
          _mm_store_si128(bvp, v);                                    \
          vc = _mm_cmpeq_epi16(vc, v);                                \
          n = _mm_andnot_si128(n, np);                                \
          n = _mm_sub_epi16(n, vc);                                   \
          _mm_store_si128(bvp+TO, n);                                 \
          vc = _mm_and_si128(vc, r4);                                 \
          rp = _mm_max_epi16(rp, vc);                                 \
          _mm_store_si128(bvp+TO*2, rp);                              \
      }                                                               \
      else {                                                          \
          n = _mm_cmpgt_epi16(n, v);                                  \
          v = _mm_max_epi16(v, vc);                                   \
          _mm_store_si128(bvp, v);                                    \
          vc = _mm_cmpeq_epi16(vc, v);                                \
          rp = _mm_or_si128(rp, n);                                   \
          n = _mm_andnot_si128(n, np);                                \
          n = _mm_sub_epi16(n, vc);                                   \
          _mm_store_si128(bvp+TO, n);                                 \
          rp = _mm_max_epi16(rp, r4);                                 \
          _mm_store_si128(bvp+TO*2, rp);                              \
      }                                                               \
    }
#define STEPqec1p2()                                                    \
    {                                                                   \
      __m128i M1 = _mm_loadu_si128(arp);                              \
      __m128i M2 = _mm_loadu_si128(arp+1);                            \
      __m128i wb1 = *((volatile __m128i*)(Hrp+TO*4));                 \
      __m128i wb2 = *((volatile __m128i*)(Hrp+TO*4+1));               \
      __m128i m_ma_over_mi = *((volatile __m128i*)&SNS::m_ma_over_mi); \
      __m128i m_mi = *((volatile __m128i*)&SNS::m_mi);                \
      M1 = _mm_cmpeq_epi16(M1, wb1);                                  \
      M2 = _mm_cmpeq_epi16(M2, wb2);                                  \
      __m128i Rp1 = _mm_loadu_si128((__m128i*)((i16*)Rrp-1));         \
      __m128i Rp2 = _mm_loadu_si128((__m128i*)((i16*)Rrp+7));         \
      M1 = _mm_and_si128(M1, m_ma_over_mi);                           \
      M2 = _mm_and_si128(M2, m_ma_over_mi);                           \
      M1 = _mm_add_epi16(M1, m_mi);                                   \
      M2 = _mm_add_epi16(M2, m_mi);                                   \
      M1 = _mm_add_epi16(M1, Rp1);                                    \
      M2 = _mm_add_epi16(M2, Rp2);                                    \
      __m128i H1 = _mm_loadu_si128((__m128i*)((i16*)Hrp-1));          \
      __m128i H2 = _mm_loadu_si128((__m128i*)((i16*)Hrp+7));          \
      if ( qsc ) {                                                    \
          __m128i zero = _mm_setzero_si128();                         \
          M1 = _mm_max_epi16(M1, zero);                               \
          M2 = _mm_max_epi16(M2, zero);                               \
      }                                                               \
      __m128i R1 = *((volatile __m128i*)(Hrp+TO));                    \
      __m128i R2 = *((volatile __m128i*)(Hrp+TO+1));                  \
      R1 = _mm_max_epi16(R1, H1);                                     \
      R2 = _mm_max_epi16(R2, H2);                                     \
      R1 = _mm_max_epi16(R1, M1);                                     \
      R2 = _mm_max_epi16(R2, M2);                                     \
      _mm_store_si128(Rrp, R1);                                       \
      _mm_store_si128(Rrp+1, R2);                                     \
      __m128i m_oe = *((volatile __m128i*)&SNS::m_oe);                \
      __m128i m_e = *((volatile __m128i*)&SNS::m_e);                  \
      __m128i V1 = *((volatile __m128i*)(Hrp+TO));                    \
      __m128i V2 = *((volatile __m128i*)(Hrp+TO+1));                  \
      M1 = _mm_add_epi16(M1, m_oe);                                   \
      M2 = _mm_add_epi16(M2, m_oe);                                   \
      H1 = _mm_add_epi16(H1, m_e);                                    \
      H2 = _mm_add_epi16(H2, m_e);                                    \
      V1 = _mm_add_epi16(V1, m_e);                                    \
      V2 = _mm_add_epi16(V2, m_e);                                    \
      H1 = _mm_max_epi16(H1, M1);                                     \
      H2 = _mm_max_epi16(H2, M2);                                     \
      V1 = _mm_max_epi16(V1, M1);                                     \
      V2 = _mm_max_epi16(V2, M2);                                     \
      _mm_store_si128(Hrp, H1);                                       \
      _mm_store_si128(Hrp+1, H2);                                     \
      _mm_store_si128(Hrp+TO, V1);                                    \
      _mm_store_si128(Hrp+TO+1, V2);                                  \
      abort();                                                        \
      __m128i n = R;                                                  \
      __m128i vc = R;                                                 \
      __m128i v = _mm_load_si128(bvp);                                \
      __m128i rp = _mm_load_si128(bvp+TO*2);                          \
      __m128i np = _mm_load_si128(bvp+TO);                            \
      if ( dir ) {                                                    \
          n = _mm_cmpgt_epi16(n, v);                                  \
          v = _mm_max_epi16(v, vc);                                   \
          _mm_store_si128(bvp, v);                                    \
          vc = _mm_cmpeq_epi16(vc, v);                                \
          n = _mm_andnot_si128(n, np);                                \
          n = _mm_sub_epi16(n, vc);                                   \
          _mm_store_si128(bvp+TO, n);                                 \
          vc = _mm_and_si128(vc, r4);                                 \
          rp = _mm_max_epi16(rp, vc);                                 \
          _mm_store_si128(bvp+TO*2, rp);                              \
      }                                                               \
      else {                                                          \
          n = _mm_cmpgt_epi16(n, v);                                  \
          v = _mm_max_epi16(v, vc);                                   \
          _mm_store_si128(bvp, v);                                    \
          vc = _mm_cmpeq_epi16(vc, v);                                \
          rp = _mm_or_si128(rp, n);                                   \
          n = _mm_andnot_si128(n, np);                                \
          n = _mm_sub_epi16(n, vc);                                   \
          _mm_store_si128(bvp+TO, n);                                 \
          rp = _mm_max_epi16(rp, r4);                                 \
          _mm_store_si128(bvp+TO*2, rp);                              \
      }                                                               \
    }

#define STEPqec1()                                              \
    {                                                           \
      __m128i M = _mm_loadu_si128((__m128i*)(ar-(r-2)+j));    \
      __m128i wb = _mm_load_si128((__m128i*)(bf+j));          \
      __m128i Rp = _mm_loadu_si128((__m128i*)(Rr+j-1));       \
      __m128i H = _mm_loadu_si128((__m128i*)(Hr+j-1));        \
      M = _mm_cmpeq_epi16(M, wb);                             \
      M = _mm_and_si128(M, m_ma_over_mi);                     \
      M = _mm_add_epi16(M, m_mi);                             \
      M = _mm_add_epi16(M, Rp);                               \
      if ( qsc ) M = _mm_max_epi16(M, _mm_setzero_si128());   \
      __m128i R = _mm_load_si128((__m128i*)(Vr+j));           \
      __m128i V = _mm_load_si128((__m128i*)(Vr+j));           \
      R = _mm_max_epi16(R, H);                                \
      R = _mm_max_epi16(R, M);                                \
      _mm_store_si128((__m128i*)(Rr+j), R);                   \
      M = _mm_add_epi16(M, m_oe);                             \
      V = _mm_add_epi16(V, m_e);                              \
      H = _mm_add_epi16(H, m_e);                              \
      V = _mm_max_epi16(V, M);                                \
      H = _mm_max_epi16(H, M);                                \
      _mm_store_si128((__m128i*)(Vr+j), V);                   \
      _mm_store_si128((__m128i*)(Hr+j), H);                   \
        {                                                       \
          __m128i rv = _mm_load_si128((__m128i*)(bv+j));      \
          __m128i rr = _mm_load_si128((__m128i*)(br+j));      \
          __m128i gt = _mm_cmpgt_epi16(R, rv);                \
          rv = _mm_max_epi16(rv, R);                          \
          _mm_store_si128((__m128i*)(bv+j), rv);              \
          __m128i lt = _mm_cmplt_epi16(R, rv);                \
          __m128i inc_n = _mm_load_si128(&m_one);             \
          inc_n = _mm_add_epi16(inc_n, lt);                   \
          if ( dir ) {                                        \
              lt = _mm_or_si128(lt, r4);                      \
              rr = _mm_max_epi16(rr, lt);                     \
          }                                                   \
          else {                                              \
              rr = _mm_or_si128(rr, gt);                      \
              rr = _mm_max_epi16(rr, r4);                     \
          }                                                   \
          __m128i rn = _mm_load_si128((__m128i*)(bn+j));      \
          rn = _mm_andnot_si128(gt, rn);                      \
          _mm_store_si128((__m128i*)(br+j), rr);              \
          rn = _mm_add_epi16(rn, inc_n);                      \
          _mm_store_si128((__m128i*)(bn+j), rn);              \
        }                                                       \
    }
#define STEPr()                                         \
  if ( dir ) {                                        \
      __m128i n = _mm_load_si128((__m128i*)(Rr+j));   \
      __m128i vc = _mm_load_si128((__m128i*)(Rr+j));  \
      __m128i v = _mm_load_si128((__m128i*)(bv+j));   \
      __m128i rp = _mm_load_si128((__m128i*)(br+j));  \
      __m128i np = _mm_load_si128((__m128i*)(bn+j));  \
      n = _mm_cmpgt_epi16(n, v);                      \
      v = _mm_max_epi16(v, vc);                       \
      _mm_store_si128((__m128i*)(bv+j), v);           \
      vc = _mm_cmpeq_epi16(vc, v);                    \
      n = _mm_andnot_si128(n, np);                    \
      n = _mm_sub_epi16(n, vc);                       \
      _mm_store_si128((__m128i*)(bn+j), n);           \
      vc = _mm_and_si128(vc, r4);                     \
      rp = _mm_max_epi16(rp, vc);                     \
      _mm_store_si128((__m128i*)(br+j), rp);          \
  }                                                   \
  else {                                              \
      __m128i n = _mm_load_si128((__m128i*)(Rr+j));   \
      __m128i vc = _mm_load_si128((__m128i*)(Rr+j));  \
      __m128i v = _mm_load_si128((__m128i*)(bv+j));   \
      __m128i rp = _mm_load_si128((__m128i*)(br+j));  \
      __m128i np = _mm_load_si128((__m128i*)(bn+j));  \
      n = _mm_cmpgt_epi16(n, v);                      \
      v = _mm_max_epi16(v, vc);                       \
      _mm_store_si128((__m128i*)(bv+j), v);           \
      vc = _mm_cmpeq_epi16(vc, v);                    \
      rp = _mm_or_si128(rp, n);                       \
      n = _mm_andnot_si128(n, np);                    \
      n = _mm_sub_epi16(n, vc);                       \
      _mm_store_si128((__m128i*)(bn+j), n);           \
      rp = _mm_max_epi16(rp, r4);                     \
      _mm_store_si128((__m128i*)(br+j), rp);          \
  }
#define STEPrp()                                \
  if ( dir ) {                                \
      __m128i n = _mm_load_si128(Rrp);        \
      __m128i vc = _mm_load_si128(Rrp);       \
      __m128i v = _mm_load_si128(bvp);        \
      __m128i rp = _mm_load_si128(bvp+TO*2);  \
      __m128i np = _mm_load_si128(bvp+TO);    \
      n = _mm_cmpgt_epi16(n, v);              \
      v = _mm_max_epi16(v, vc);               \
      _mm_store_si128(bvp, v);                \
      vc = _mm_cmpeq_epi16(vc, v);            \
      n = _mm_andnot_si128(n, np);            \
      n = _mm_sub_epi16(n, vc);               \
      _mm_store_si128(bvp+TO, n);             \
      vc = _mm_and_si128(vc, r4);             \
      rp = _mm_max_epi16(rp, vc);             \
      _mm_store_si128(bvp+TO*2, rp);          \
  }                                           \
  else {                                      \
      __m128i n = _mm_load_si128(Rrp);        \
      __m128i vc = _mm_load_si128(Rrp);       \
      __m128i v = _mm_load_si128(bvp);        \
      __m128i rp = _mm_load_si128(bvp+TO*2);  \
      __m128i np = _mm_load_si128(bvp+TO);    \
      n = _mm_cmpgt_epi16(n, v);              \
      v = _mm_max_epi16(v, vc);               \
      _mm_store_si128(bvp, v);                \
      vc = _mm_cmpeq_epi16(vc, v);            \
      rp = _mm_or_si128(rp, n);               \
      n = _mm_andnot_si128(n, np);            \
      n = _mm_sub_epi16(n, vc);               \
      _mm_store_si128(bvp+TO, n);             \
      rp = _mm_max_epi16(rp, r4);             \
      _mm_store_si128(bvp+TO*2, rp);          \
  }
#define STEPrp2()                                   \
  if ( dir ) {                                    \
      __m128i n1 = _mm_load_si128(Rrp-1);         \
      __m128i n2 = _mm_load_si128(Rrp  );         \
      __m128i v1 = _mm_load_si128(bvp-1);         \
      __m128i v2 = _mm_load_si128(bvp  );         \
      __m128i vc1 = _mm_load_si128(Rrp-1);        \
      __m128i vc2 = _mm_load_si128(Rrp  );        \
      n1 = _mm_cmpgt_epi16(n1, v1);               \
      n2 = _mm_cmpgt_epi16(n2, v2);               \
      v1 = _mm_max_epi16(v1, vc1);                \
      v2 = _mm_max_epi16(v2, vc2);                \
      _mm_store_si128(bvp-1, v1);                 \
      _mm_store_si128(bvp  , v2);                 \
      vc1 = _mm_cmpeq_epi16(vc1, v1);             \
      vc2 = _mm_cmpeq_epi16(vc2, v2);             \
      __m128i np1 = _mm_load_si128(bvp+TO-1);     \
      __m128i np2 = _mm_load_si128(bvp+TO  );     \
      n1 = _mm_andnot_si128(n1, np1);             \
      n2 = _mm_andnot_si128(n2, np2);             \
      n1 = _mm_sub_epi16(n1, vc1);                \
      n2 = _mm_sub_epi16(n2, vc2);                \
      __m128i rp1 = _mm_load_si128(bvp+TO*2-1);   \
      __m128i rp2 = _mm_load_si128(bvp+TO*2  );   \
      _mm_store_si128(bvp+TO-1, n1);              \
      _mm_store_si128(bvp+TO  , n2);              \
      vc1 = _mm_and_si128(vc1, r4);               \
      vc2 = _mm_and_si128(vc2, r4);               \
      rp1 = _mm_max_epi16(rp1, vc1);              \
      rp2 = _mm_max_epi16(rp2, vc2);              \
      _mm_store_si128(bvp+TO*2-1, rp1);           \
      _mm_store_si128(bvp+TO*2  , rp2);           \
  }                                               \
  else {                                          \
      __m128i n1 = _mm_load_si128(Rrp-1);         \
      __m128i n2 = _mm_load_si128(Rrp  );         \
      __m128i v1 = _mm_load_si128(bvp-1);         \
      __m128i v2 = _mm_load_si128(bvp  );         \
      __m128i vc1 = _mm_load_si128(Rrp-1);        \
      __m128i vc2 = _mm_load_si128(Rrp  );        \
      n1 = _mm_cmpgt_epi16(n1, v1);               \
      n2 = _mm_cmpgt_epi16(n2, v2);               \
      v1 = _mm_max_epi16(v1, vc1);                \
      v2 = _mm_max_epi16(v2, vc2);                \
      _mm_store_si128(bvp-1, v1);                 \
      _mm_store_si128(bvp  , v2);                 \
      vc1 = _mm_cmpeq_epi16(vc1, v1);             \
      vc2 = _mm_cmpeq_epi16(vc2, v2);             \
      __m128i rp1 = _mm_load_si128(bvp+TO*2-1);   \
      __m128i rp2 = _mm_load_si128(bvp+TO*2  );   \
      __m128i np1 = _mm_load_si128(bvp+TO-1);     \
      __m128i np2 = _mm_load_si128(bvp+TO  );     \
      rp1 = _mm_or_si128(rp1, n1);                \
      rp2 = _mm_or_si128(rp2, n2);                \
      n1 = _mm_andnot_si128(n1, np1);             \
      n2 = _mm_andnot_si128(n2, np2);             \
      n1 = _mm_sub_epi16(n1, vc1);                \
      n2 = _mm_sub_epi16(n2, vc2);                \
      _mm_store_si128(bvp+TO-1, n1);              \
      _mm_store_si128(bvp+TO  , n2);              \
      rp1 = _mm_max_epi16(rp1, r4);               \
      rp2 = _mm_max_epi16(rp2, r4);               \
      _mm_store_si128(bvp+TO*2-1, rp1);           \
      _mm_store_si128(bvp+TO*2  , rp2);           \
  }
#define STEP() if ( qec ) { STEPqec1(); } else { STEPqec0(); }

  if ( m < n ) {
      const int use_qec1 = 0;
      const int use_stepqec0p2 = 1;
      const int use_steprp2 = 0;

      int n_oe = oe;
      R0r[-1] = 0;
      CHECK_ROW(0,R,R0r);
      if ( qsc ) {
          R1r[-1] = 0; // left
      }
      else {
          R1r[-1] = n_oe; // left
      }
      Hr[-1] = -INF; // top
      Vr[-1] = -INF; // left
      Vr[0] = -INF; // top
      R1r[0] = 0;
      CHECK_ROW(1,H,Hr-1);
      CHECK_ROW(1,V,Vr);
      CHECK_ROW(1,R,R1r);
      if ( qec ) {
          bv[0] = -INF;
          br[0] = -INF;
      }
      for ( size_t r = 2; r <= m; ++r ) {
          i16* Rr = Rtr(r);
          if ( qsc ) {
              Hr[-1] = oe;
          }
          else {
              n_oe += e;
              Hr[-1] = n_oe; // left and j==1
          }
          Hr[r-2] = -INF; // top i==0
          Vr[r-2] = oe; // i==1
          Vr[r-1] = -INF; // top i==0
          CHECK_ROW(r,H,Hr-1);
          CHECK_ROW(r,V,Vr);
          CHECK_ROW(r-2,R,Rr);
          const int js = ((r-2)&-8);
          const int je = 0;
            {
              __m128i* arp = (__m128i*)(ar-(r-2)+js);
              __m128i* Rrp = (__m128i*)(Rr+js);
              __m128i* Hrp = (__m128i*)(Hr+js);
              __m128i* Hre = (__m128i*)(Hr+je);
              if ( use_qec1 && qec ) {
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  do {
                      STEPqec1p();
                      --Hrp;
                      --arp;
                      --Rrp;
                      --bvp;
                  } while ( Hrp >= Hre );
                  bv[r-1] = -INF;
                  br[r-1] = -INF;
              }
              else {
                  do {
                      if ( use_stepqec0p2 && Hrp > Hre ) {
                          STEPqec0p2();
                          Hrp -= 2;
                          arp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPqec0p();
                          --Hrp;
                          --arp;
                          --Rrp;
                      }
                  } while ( Hrp >= Hre );
              }
            }
          if ( !qsc ) {
              Rr[-1] = n_oe; // left
          }
          Rr[r-1] = 0;
          CHECK_ROW(r,R,Rr);
          if ( qec ) {
              if ( !use_qec1 ) {
                  const int js = ((r-3)&-8);
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  __m128i* bve = (__m128i*)(bv+je);
                  __m128i* Rrp = (__m128i*)(Rr+js);
                  do {
                      if ( use_steprp2 && bvp > bve ) {
                          STEPrp2();
                          bvp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPrp();
                          --bvp;
                          --Rrp;
                      }
                  } while ( bvp >= bve );
                  bv[r-1] = -INF;
                  br[r-1] = -INF;
              }
          }
      }

      if ( qsc ) {
          Hr[-1] = oe;
      }
      else {
          Hr[-1] = n_oe+e; // j==1
      }
      for ( size_t r = m+1; r <= n; ++r ) {
          i16* Rr = Rtr(r);
          Hr[r-2] = -INF; // top i==0
          Vr[r-2] = oe; // i==1
          Vr[r-1] = -INF; // top i==0
          CHECK_ROW(r,H,Hr-1);
          CHECK_ROW(r,V,Vr);
          CHECK_ROW(r-2,R,Rr);
          const int js = ((r-2)&-8);
          const int je = (r-m-1)&-8;
            {
              __m128i* arp = (__m128i*)(ar-(r-2)+js);
              __m128i* Rrp = (__m128i*)(Rr+js);
              __m128i* Hrp = (__m128i*)(Hr+js);
              __m128i* Hre = (__m128i*)(Hr+je);
              if ( use_qec1 && qec ) {
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  do {
                      STEPqec1p();
                      --Hrp;
                      --arp;
                      --Rrp;
                      --bvp;
                  } while ( Hrp >= Hre );
                  bv[r-1] = -INF;
                  br[r-1] = -INF;
              }
              else {
                  do {
                      if ( use_stepqec0p2 && Hrp > Hre ) {
                          STEPqec0p2();
                          Hrp -= 2;
                          arp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPqec0p();
                          --Hrp;
                          --arp;
                          --Rrp;
                      }
                  } while ( Hrp >= Hre );
              }
              Rr[r-1] = 0;
            }
          CHECK_ROW(r,R,Rr);
          if ( qec ) {
              if ( !use_qec1 ) {
                  const int js = ((r-3)&-8);
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  __m128i* bve = (__m128i*)(bv+je);
                  __m128i* Rrp = (__m128i*)(Rr+js);
                  do {
                      if ( use_steprp2 && bvp > bve ) {
                          STEPrp2();
                          bvp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPrp();
                          --bvp;
                          --Rrp;
                      }
                  } while ( bvp >= bve );
                  bv[r-1] = -INF;
                  br[r-1] = -INF;
              }
          }
          else {
              size_t j = r-m, i = m, k = j-1;
              int v = Rr[k];
              if ( v > opt ) {
                  opt = v;
                  n_best = 1;
                  best_i = i;
                  best_j = j;
              }
              else if ( v == opt ) {
                  n_best += 1;
                  /*
                  if ( dir ) {
                      best_i = i;
                      best_j = j;
                  }
                  */
                  best_i = i;
                  best_j = j;
              }
          }
      }

      Vr[n-1] = oe;
      for ( size_t r = n+1; r <= n+m; ++r ) {
          i16* Rr = Rtr(r);
          CHECK_ROW(r,H,Hr-1);
          CHECK_ROW(r,V,Vr);
          CHECK_ROW(r-2,R,Rr);
          const int js = ((n-1)&-8);
          const int je = (r-m-1)&-8;
            {
              __m128i* arp = (__m128i*)(ar-(r-2)+js);
              __m128i* Rrp = (__m128i*)(Rr+js);
              __m128i* Hrp = (__m128i*)(Hr+js);
              __m128i* Hre = (__m128i*)(Hr+je);
              if ( use_qec1 && qec ) {
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  do {
                      STEPqec1p();
                      --Hrp;
                      --arp;
                      --Rrp;
                      --bvp;
                  } while ( Hrp >= Hre );
              }
              else {
                  do {
                      if ( use_stepqec0p2 && Hrp > Hre ) {
                          STEPqec0p2();
                          Hrp -= 2;
                          arp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPqec0p();
                          --Hrp;
                          --arp;
                          --Rrp;
                      }
                  } while ( Hrp >= Hre );
              }
            }
          CHECK_ROW(r,R,Rr);
          if ( qec ) {
              if ( !use_qec1 ) {
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  __m128i* bve = (__m128i*)(bv+je);
                  __m128i* Rrp = (__m128i*)(Rr+js);
                  do {
                      if ( use_steprp2 && bvp > bve ) {
                          STEPrp2();
                          bvp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPrp();
                          --bvp;
                          --Rrp;
                      }
                  } while ( bvp >= bve );
              }
          }
          else {
              size_t j = r-m, i = m, k = j-1;
              int v = Rr[k];
              if ( v > opt ) {
                  opt = v;
                  n_best = 1;
                  best_i = i;
                  best_j = j;
              }
              else if ( v == opt ) {
                  n_best += 1;
                  /*
                  if ( dir ) {
                      best_i = i;
                      best_j = j;
                  }
                  */
                  best_i = i;
                  best_j = j;
              }
          }
      }
  }
  else {
      const int use_qec1 = 0;
      const int use_stepqec0p2 = 1;
      const int use_steprp2 = 0;

      int n_oe = oe;
      R0r[-1] = 0; // left
      CHECK_ROW(0,R,R0r);
      if ( qsc ) {
          R1r[-1] = 0; // left
      }
      else {
          R1r[-1] = n_oe; // left
      }
      Hr[-1] = -INF; // top
      Vr[-1] = -INF; // left
      Vr[0] = -INF; // top
      R1r[0] = 0; // top
      CHECK_ROW(1,H,Hr-1);
      CHECK_ROW(1,V,Vr);
      CHECK_ROW(1,R,R1r);
      if ( qec ) {
          bv[0] = -INF;
          br[0] = -INF;
      }
      for ( size_t r = 2; r <= n; ++r ) {
          i16* Rr = Rtr(r);
          if ( qsc ) {
              Hr[-1] = oe;
          }
          else {
              n_oe += e;
              Hr[-1] = n_oe; // left and j==1
          }
          Hr[r-2] = -INF; // top i==0
          Vr[r-2] = oe; // i==1
          Vr[r-1] = -INF; // top i==0
          CHECK_ROW(r,H,Hr-1);
          CHECK_ROW(r,V,Vr);
          CHECK_ROW(r-2,R,Rr);
          const int js = ((r-2)&-8);
          const int je = 0;
            {
              __m128i* arp = (__m128i*)(ar-(r-2)+js);
              __m128i* Rrp = (__m128i*)(Rr+js);
              __m128i* Hrp = (__m128i*)(Hr+js);
              __m128i* Hre = (__m128i*)(Hr+je);
              if ( use_qec1 && qec ) {
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  do {
                      STEPqec1p();
                      --Hrp;
                      --arp;
                      --Rrp;
                      --bvp;
                  } while ( Hrp >= Hre );
                  bv[r-1] = -INF;
                  br[r-1] = -INF;
              }
              else {
                  do {
                      if ( use_stepqec0p2 && Hrp > Hre ) {
                          STEPqec0p2();
                          Hrp -= 2;
                          arp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPqec0p();
                          --Hrp;
                          --arp;
                          --Rrp;
                      }
                  } while ( Hrp >= Hre );
              }
            }
          if ( !qsc ) {
              Rr[-1] = n_oe; // left
          }
          Rr[r-1] = 0;
          CHECK_ROW(r,R,Rr);
          if ( qec ) {
              bv[r-1] = -INF;
              br[r-1] = -INF;
          }
          if ( qec ) {
              if ( !use_qec1 ) {
                  const int js = ((r-3)&-8);
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  __m128i* bve = (__m128i*)(bv+je);
                  __m128i* Rrp = (__m128i*)(Rr+js);
                  do {
                      if ( use_steprp2 && bvp > bve ) {
                          STEPrp2();
                          bvp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPrp();
                          --bvp;
                          --Rrp;
                      }
                  } while ( bvp >= bve );
                  bv[r-1] = -INF;
                  br[r-1] = -INF;
              }
          }
      }

      Vr[n-1] = oe; // i==1
      for ( size_t r = n+1; r <= m; ++r ) {
          i16* Rr = Rtr(r);
          if ( !qsc ) {
              n_oe += e;
              Hr[-1] = n_oe; // left and j==1
          }
          CHECK_ROW(r,H,Hr-1);
          CHECK_ROW(r,V,Vr);
          CHECK_ROW(r-2,R,Rr);
          const int js = ((n-1)&-8);
          const int je = 0;
            {
              __m128i* arp = (__m128i*)(ar-(r-2)+js);
              __m128i* Rrp = (__m128i*)(Rr+js);
              __m128i* Hrp = (__m128i*)(Hr+js);
              __m128i* Hre = (__m128i*)(Hr+je);
              if ( use_qec1 && qec ) {
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  do {
                      STEPqec1p();
                      --Hrp;
                      --arp;
                      --Rrp;
                      --bvp;
                  } while ( Hrp >= Hre );
                  bv[r-1] = -INF;
                  br[r-1] = -INF;
              }
              else {
                  do {
                      if ( use_stepqec0p2 && Hrp > Hre ) {
                          STEPqec0p2();
                          Hrp -= 2;
                          arp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPqec0p();
                          --Hrp;
                          --arp;
                          --Rrp;
                      }
                  } while ( Hrp >= Hre );
              }
            }
          if ( !qsc ) {
              Rr[-1] = n_oe; // left
          }
          CHECK_ROW(r,R,Rr);
          if ( qec ) {
              if ( !use_qec1 ) {
                  const int js = ((r-3)&-8);
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  __m128i* bve = (__m128i*)(bv+je);
                  __m128i* Rrp = (__m128i*)(Rr+js);
                  do {
                      if ( use_steprp2 && bvp > bve ) {
                          STEPrp2();
                          bvp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPrp();
                          --bvp;
                          --Rrp;
                      }
                  } while ( bvp >= bve );
                  bv[r-1] = -INF;
                  br[r-1] = -INF;
              }
          }
      }

      if ( !qsc ) {
          Hr[-1] = n_oe+e; // j==1
      }
      for ( size_t r = m+1; r <= n+m; ++r ) {
          i16* Rr = Rtr(r);
          CHECK_ROW(r,H,Hr-1);
          CHECK_ROW(r,V,Vr);
          CHECK_ROW(r-2,R,Rr);
          const int js = ((n-1)&-8);
          const int je = (r-m-1)&-8;
            {
              __m128i* arp = (__m128i*)(ar-(r-2)+js);
              __m128i* Rrp = (__m128i*)(Rr+js);
              __m128i* Hrp = (__m128i*)(Hr+js);
              __m128i* Hre = (__m128i*)(Hr+je);
              if ( use_qec1 && qec ) {
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  do {
                      STEPqec1p();
                      --Hrp;
                      --arp;
                      --Rrp;
                      --bvp;
                  } while ( Hrp >= Hre );
              }
              else {
                  do {
                      if ( use_stepqec0p2 && Hrp > Hre ) {
                          STEPqec0p2();
                          Hrp -= 2;
                          arp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPqec0p();
                          --Hrp;
                          --arp;
                          --Rrp;
                      }
                  } while ( Hrp >= Hre );
              }
            }
          CHECK_ROW(r,R,Rr);
          if ( qec ) {
              if ( !use_qec1 ) {
                  __m128i r4 = _mm_set1_epi16(dir? r: INF-r);
                  __m128i* bvp = (__m128i*)(bv+js);
                  __m128i* bve = (__m128i*)(bv+je);
                  __m128i* Rrp = (__m128i*)(Rr+js);
                  do {
                      if ( use_steprp2 && bvp > bve ) {
                          STEPrp2();
                          bvp -= 2;
                          Rrp -= 2;
                      }
                      else {
                          STEPrp();
                          --bvp;
                          --Rrp;
                      }
                  } while ( bvp >= bve );
              }
          }
          else {
              size_t j = r-m, i = m, k = j-1;
              int v = Rr[k];
              if ( v > opt ) {
                  opt = v;
                  n_best = 1;
                  best_i = i;
                  best_j = j;
              }
              else if ( v == opt ) {
                  n_best += 1;
                  /*
                  if ( dir ) {
                      best_i = i;
                      best_j = j;
                  }
                  */
                  best_i = i;
                  best_j = j;
              }
          }
      }
  }

  if ( qec ) {
      for ( size_t k = 0; k < n; ++k ) {
          size_t r = dir? br[k]: INF-br[k];
          size_t j = k+1, i = r-j;
          int v = bv[k];
          if ( v > opt ) {
              opt = v;
              n_best = bn[k];
              best_i = i;
              best_j = j;
          }
          else if ( v == opt ) {
              n_best += bn[k];
              /*
              if ( dir ) {
                  if ( i > best_i || i == best_i && j > best_j ) {
                      best_i = i;
                      best_j = j;
                  }
              }
              else {
                  if ( i < best_i || i == best_i && j < best_j ) {
                      best_i = i;
                      best_j = j;
                  }
              }
              */
              if ( dir ) {
                  if ( i > best_i || i == best_i && j > best_j ) {
                      best_i = i;
                      best_j = j;
                  }
              }
              else {
                  if ( i < best_i || i == best_i && j > best_j ) {
                      best_i = i;
                      best_j = j;
                  }
              }
          }
      }
  }
  END_CLOCK(clock[qsc][qec][dir],n,m);
}

inline char* format_int(char* p, int v, char sep)
{
  if ( sep ) *--p = sep;
  bool neg = v < 0;
  if ( neg ) v = -v;
  do {
      int d = v%10;
      v = v/10;
      *--p = '0'+d;
  } while ( v );
  if ( neg ) *--p = '-';
  return p;
}

char ret_buffer[1024];
char* const ret_end = ret_buffer+sizeof(ret_buffer)-1;
char* ret_begin;

inline void process(const string& t, const string& q, int qsc, int qec,
                    int ma, int mi, int go, int ge, int dir,
                    int *_opt, int *_te, int *_qe, int *_n_best) {
  m = q.size();
  a = q.data();

  n = t.size();
  b = t.data();

  SNS::ma = ma;
  SNS::mi = mi;
  SNS::go = go;
  SNS::ge = ge;
  m_e = _mm_set1_epi16(ge);
  m_o = _mm_set1_epi16(go);
  m_oe = _mm_set1_epi16(go+ge);
  m_mi = _mm_set1_epi16(mi);
  m_ma = _mm_set1_epi16(ma);
  m_ma_over_mi = _mm_set1_epi16(ma-mi);
  m_one = _mm_set1_epi16(1);

  forN ( k, 8 ) {
      bf[k-8] = 'b'<<8;
      bf[k+n] = 'b'<<8;
      ar[-(k-8)] = 'a'<<8;
      ar[-(k+m)] = 'a'<<8;
  }
  forN ( i, n ) bf[i] = b[i]<<8;
  forN ( j, m ) ar[-j] = a[j]<<8;

  opt = -INF;
  n_best = 0;
  best_i = best_j = 0;

  if ( qsc ) {
      if ( qec ) {
          if ( dir ) {
              process<1,1,1>(n, m);
          }
          else {
              process<1,1,0>(n, m);
          }
      }
      else {
          if ( dir ) {
              process<1,0,1>(n, m);
          }
          else {
              process<1,0,0>(n, m);
          }
      }
  }
  else {
      if ( qec ) {
          if ( dir ) {
              process<0,1,1>(n, m);
          }
          else {
              process<0,1,0>(n, m);
          }
      }
      else {
          if ( dir ) {
              process<0,0,1>(n, m);
          }
          else {
              process<0,0,0>(n, m);
          }
      }
  }
  (*_n_best) = n_best;
  (*_qe) = best_j-1;
  (*_te) = best_i-1;
  (*_opt) = opt;
  /*
  char* p = ret_end;
  p = format_int(p, SNS::n_best, 0);
  p = format_int(p, SNS::best_j-1, ' ');
  p = format_int(p, SNS::best_i-1, ' ');
  p = format_int(p, SNS::opt, ' ');
  ret_begin = p;
  */
}

double start_time;

void init()
{
#ifdef HOME_RUN
  static int skip_params;
  if ( !skip_params++ ) {
      OUT_PARAM(SUBMISSION);
      OUT_PARAM(SCORE);
  }
#endif
#ifdef CLOCK
  start_time = get_time();
  forN ( qsc, 2 ) forN ( qec, 2 ) forN ( dir, 2 )
    INIT_CLOCK(clock[qsc][qec][dir]);
  INIT_CLOCK(clock_full);
#endif
  //cout <<" tmp: "<<tmp<<endl;
}

void close()
{
  PRINT_CLOCK(clock[0][0][0]);
  PRINT_CLOCK(clock[0][0][1]);
  PRINT_CLOCK(clock[0][1][0]);
  PRINT_CLOCK(clock[0][1][1]);
  PRINT_CLOCK(clock[1][0][0]);
  PRINT_CLOCK(clock[1][0][1]);
  PRINT_CLOCK(clock[1][1][0]);
  PRINT_CLOCK(clock[1][1][1]);
  PRINT_CLOCK(clock_full);
#ifdef CLOCK
  cout << " real time: "<<get_time()-start_time<<" secs"<<endl;
#endif
}

void cleanup()
{
}

KET;

Solution8::Solution8() : printed(0) {
    SNS::init();
    max_qlen = 512;
    max_tlen = 1024;
}

Solution8::~Solution8() {
    SNS::close();
}
int Solution8::process(const string& t, const string& q, int qsc, int qec,
               int ma, int mi, int go, int ge, int dir,
               int *opt, int *te, int *qe, int *n_best) {
#if defined(CLOCK) && !defined(HOME_RUN)
    if ( !printed && SNS::clock_full.clock >= 20*CLOCK_FREQ ) {
        printed = 1;
        SNS::close();
    }
#endif
    LOCAL_CLOCK(SNS::clock_full,t.size(),q.size());
    SNS::process(t, q, qsc, qec, ma, mi, go, ge, dir, opt, te, qe, n_best);
    return (*opt);
    //return string(SNS::ret_begin, SNS::ret_end);
}

