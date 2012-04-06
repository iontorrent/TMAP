#ifndef SOLUTION7_H_
#define SOLUTION7_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

// Coder: logicmachine
// Submission: 20
// URL: http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=22887412&subnum=20
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <xmmintrin.h>
#include <emmintrin.h>
#include "solution.h"

using namespace std;

#define mm_set1_epi8(x, a) { \
        x = _mm_insert_epi16(x, a, 0); \
        x = _mm_unpacklo_epi8(x, x); \
        x = _mm_shufflelo_epi16(x, 0); \
        x = _mm_shuffle_epi32(x, 0); \
    }
#define mm_set1_epi16(x, a) { \
        x = _mm_insert_epi16(x, a, 0); \
        x = _mm_shufflelo_epi16(x, 0); \
        x = _mm_shuffle_epi32(x, 0); \
    }
#define mm_reduction_max_epu8(dst, src) { \
        __m128i tmp = _mm_max_epu8(_mm_srli_si128(src, 1), src); \
        tmp = _mm_max_epu8(_mm_srli_si128(tmp, 2), tmp); \
        tmp = _mm_max_epu8(_mm_srli_si128(tmp, 4), tmp); \
        tmp = _mm_max_epu8(_mm_srli_si128(tmp, 8), tmp); \
        dst = static_cast<uint8_t>(_mm_extract_epi16(tmp, 0) & 0xff); \
    }
#define mm_reduction_max_epi16(dst, src) { \
        __m128i tmp = _mm_max_epi16(_mm_srli_si128(src, 2), src); \
        tmp = _mm_max_epi16(_mm_srli_si128(tmp, 4), tmp); \
        tmp = _mm_max_epi16(_mm_srli_si128(tmp, 8), tmp); \
        dst = static_cast<int16_t>(_mm_extract_epi16(tmp, 0)); \
    }

class Solution7 : public Solution {

private:
    typedef unsigned char      uint8_t;
    typedef char               int8_t;
    typedef unsigned short     uint16_t;
    typedef short              int16_t;
    typedef unsigned int       uint32_t;
    typedef int                int32_t;
    typedef unsigned long long uint64_t;
    typedef long long          int64_t;

    typedef struct ProcessResult {
        int opt, query_end, target_end, n_best;
    } ProcessResult;

    typedef struct SubBlockU8S {
        uint8_t M[16];
        uint8_t V[16];
        uint8_t H[16];
    } SubBlockU8S;

    typedef struct BlockU8S {
        uint8_t target[16];
        SubBlockU8S data[2];
    } BlockU8S;

    typedef struct BlockI16S {
        int16_t target[8];
        int16_t M[8];
        int16_t V[8];
        int16_t H[8];
    } BlockI16S;

    typedef struct BlockI16 {
        int16_t target[8];
        int16_t NM[8];
        int16_t M[8];
        int16_t V[8];
        int16_t H[8];
    } BlockI16;

    template <typename T, int ALIGNMENT>
    inline T *malign(void *p){
        return reinterpret_cast<T *>((
            reinterpret_cast<unsigned long long>(p) + ALIGNMENT - 1) & ~(ALIGNMENT - 1));
    }

    static const int ALIGN_SIZE = 16;
    static const int MAX_TARGET = 1024;
    static const int MAX_QUERY = 512;
    static const int NEG_INF = -0x8000;

    static const int MAX_BLOCK_NUM_U8 = (MAX_TARGET + 31) >> 4;
    static const int MAX_SIMD_QUERY_SIZE_U8 = (MAX_QUERY + 31) & ~15;
    static const int WORKAREA_SIZE_U8 = MAX_TARGET * sizeof(BlockU8S) + MAX_TARGET;
    uint8_t m_workarea_u8[WORKAREA_SIZE_U8 + ALIGN_SIZE];

    static const int MAX_BLOCK_NUM_I16 = MAX_BLOCK_NUM_U8 * 2;
    static const int MAX_SIMD_QUERY_SIZE_I16 = (MAX_QUERY + 15) & ~7;
    static const int WORKAREA_SIZE_I16 =
        MAX_BLOCK_NUM_I16 * sizeof(BlockI16) +
        MAX_TARGET * sizeof(int16_t) +
        MAX_SIMD_QUERY_SIZE_I16 * sizeof(int16_t) * 8;
    int16_t m_workarea_i16[WORKAREA_SIZE_I16 + ALIGN_SIZE];

    void process_with_start_clip(
        const string &target, const string &query,
        int query_end_clip, int match_score, int mismatch_score,
        int gap_open, int gap_extension, int dir, ProcessResult &result);

    void process_without_start_clip(
        const string &target, const string &query,
        int query_end_clip, int match_score, int mismatch_score,
        int gap_open, int gap_extension, int dir,
        ProcessResult &result);

public:
    virtual int process(
        string& target, string& query,
        int query_start_clip, int query_end_clip,
        int match_score, int mismatch_score,
        int gap_open, int gap_extension, int dir,
        int *opt, int *te, int *qe, int *n_best);
};
#endif
