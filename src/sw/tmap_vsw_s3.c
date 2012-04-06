// Coder: folsena
// Submission: 18
// URL: http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=22926110&subnum=18

#include <stdlib.h>
#include <stdint.h>
#include <mmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h> // REQUIRES SSE3
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "tmap_vsw_definitions.h"
#include "tmap_vsw_s3.h"

#include "string.h"

const int32_t INF = 1073741824;
const int16_t NINF = -30000;

typedef struct {
    int32_t opt;
    int32_t query_end;
    int32_t target_end;
    int32_t n_best;
} result_t;

inline int32_t max(int a, int32_t b) {
    return a > b ? a : b;
}

static void
fill(int16_t *first, int16_t *last, int16_t value)
{
  while(first != last) *first++ = value;
}

int16_t*
reverse_copy(const uint8_t *first, const uint8_t *last, int16_t *result)
{
  while (first!=last) *result++ = *--last;
  return result;
}

static void
fill_n(int16_t *first, int32_t n, int16_t value)
{
  for (; n>0; --n)  *first++ = value;
}

static void
process(tmap_vsw_data_s3_t *vsw, const uint8_t *bs, const int32_t n, const uint8_t *as, int32_t m,
        int32_t qsc, int32_t qec,
        int32_t mm, int32_t mi, int32_t oe, int32_t e, int32_t dir,
        result_t *result)
{
    int32_t opt = NINF, query_end=-1, target_end=-1, n_best=0;
    int32_t bl, i, j, l, t;
    oe += e;

    int16_t *abuf = vsw->abuf;
    int16_t *B = vsw->B;
    int16_t *MV = vsw->MV;
    m128si16 *X = vsw->X;

    fill(B, B+8, 0);
    //fill_n(reverse_copy(bs.begin(), bs.end(), B+8), 8, 0);
    fill_n(reverse_copy(bs, bs+n, B+8), 8, 0);
    fill(MV+8, MV+8+n, oe);
    //copy(as.c_str(), as.c_str() + m, abuf);
    for(i=0;i<m;i++) {
        abuf[i] = as[i];
    }

    for(i=0; i<n; i++) {
        X[i+8].m = _mm_setzero_si128();
    }

    const __m128i mme = _mm_set1_epi16(e);
    const __m128i mmm = _mm_set1_epi16(mm);
    const __m128i mmmi = _mm_set1_epi16(mm-mi);
    const __m128i mmi = _mm_set1_epi16(mi);

    const __m128i mmoe = _mm_set1_epi16(oe);
    const __m128i mNINF = _mm_set1_epi16(NINF);
    const __m128i mnegones = _mm_set1_epi16(-1);

    const int32_t nb = (m+7)/8*8;
    const int16_t* b;
    const int16_t* a = abuf;
    __m128i A;

    if (qec == 0) {
        if (qsc == 0) {
            for (bl=0; bl < nb; bl+=8) {
                A = _mm_lddqu_si128((__m128i*)a);
                a += 8;
                b = B + n+7;
                  {
                    __m128i tmp = _mm_set1_epi16(oe + e * bl);
                    X[0].m = _mm_sub_epi16(tmp, mme);
                    X[1].m = tmp;
                    __m128i h = _mm_add_epi16(tmp, mme);
                    __m128i v = mNINF;
                    __m128i maskF = mnegones;

                    for (t=0; t < 7; t++) {
                        maskF = _mm_slli_si128(maskF, 2);

                        __m128i x = _mm_or_si128(_mm_slli_si128(X[t].m, 2), _mm_srli_si128(X[t+8].m, 14)); //common
                        v = _mm_insert_epi16(v, MV[t+8], 0); //common

                        __m128i mask = _mm_cmpeq_epi16(A, _mm_lddqu_si128((__m128i*)b)); //common
                        __m128i comp = _mm_or_si128(_mm_and_si128(mask, mmm), _mm_andnot_si128(mask, mmi)); //common

                        comp = _mm_or_si128(_mm_andnot_si128(maskF, comp), _mm_and_si128(maskF, mNINF));

                        const __m128i M = _mm_add_epi16(x, comp);
                        x = _mm_max_epi16(_mm_max_epi16(h, v), M); //common

                        X[t+2].m = x; //common

                        const __m128i Moe = _mm_add_epi16(M, mmoe); //common
                        h = _mm_max_epi16(_mm_add_epi16(h, mme), Moe); //common
                        v = _mm_max_epi16(_mm_add_epi16(v, mme), Moe); //common

                        MV[t+1] = _mm_extract_epi16(v, 7); //common
                        v = _mm_slli_si128(v, 2); //common
                        --b; //common
                    }

                    for (t=7; t < n+7; t++) {
                        __m128i x = _mm_or_si128(_mm_slli_si128(X[t].m, 2), _mm_srli_si128(X[t+8].m, 14)); //common
                        v = _mm_insert_epi16(v, MV[t+8], 0); //common

                        __m128i mask = _mm_cmpeq_epi16(A, _mm_lddqu_si128((__m128i*)b)); //common
                        const __m128i comp = _mm_add_epi16(x, mmi);
                        mask = _mm_and_si128(mask, mmmi);

                        const __m128i M = _mm_add_epi16(mask, comp);
                        x = _mm_max_epi16(_mm_max_epi16(h, v), M); //common

                        X[t+2].m = x; //common

                        const __m128i Moe = _mm_add_epi16(M, mmoe); //common
                        h = _mm_max_epi16(_mm_add_epi16(h, mme), Moe); //common
                        v = _mm_max_epi16(_mm_add_epi16(v, mme), Moe); //common

                        MV[t+1] = _mm_extract_epi16(v, 7); //common
                        v = _mm_slli_si128(v, 2); //common
                        --b; //common
                    }
                  }
            }
        }
        else {
            for (bl=0; bl < nb; bl+=8) {
                A = _mm_lddqu_si128((__m128i*)a);
                a += 8;
                b = B + n+7;
                  {
                    X[0].m = _mm_setzero_si128();
                    X[1].m = _mm_setzero_si128();
                    __m128i h = mmoe;
                    __m128i v = mmoe;

                    for (t=0; t < n+7; t++) {
                        __m128i x = _mm_or_si128(_mm_slli_si128(X[t].m, 2), _mm_srli_si128(X[t+8].m, 14)); //common
                        v = _mm_insert_epi16(v, MV[t+8], 0); //common

                        __m128i mask = _mm_cmpeq_epi16(A, _mm_lddqu_si128((__m128i*)b)); //common
                        const __m128i comp = _mm_add_epi16(x, mmi);
                        mask = _mm_and_si128(mask, mmmi);

                        const __m128i M = _mm_max_epi16(_mm_setzero_si128(), _mm_add_epi16(mask, comp));
                        x = _mm_max_epi16(_mm_max_epi16(h, v), M); //common

                        X[t+2].m = x; //common

                        const __m128i Moe = _mm_add_epi16(M, mmoe); //common
                        h = _mm_max_epi16(_mm_add_epi16(h, mme), Moe); //common
                        v = _mm_max_epi16(_mm_add_epi16(v, mme), Moe); //common

                        MV[t+1] = _mm_extract_epi16(v, 7); //common
                        v = _mm_slli_si128(v, 2); //common
                        --b; //common
                    }
                  }
            }
        }

          {
            const int32_t l = (m-1)%8;
            opt = NINF;
            n_best = 1;
            query_end = m-1;

            for (j=0; j < n; j++)
              opt = max(opt, X[j+l+2].s[l]);

            if (dir == 0) {
                int32_t j;
                for (j=0;; j++) {
                    if (X[j+l+2].s[l] == opt) {
                        break;
                    }
                }
                target_end = j;
                j++;
                for (; j < n; j++) {
                    if (X[j+l+2].s[l] == opt) {
                        n_best++;
                    }
                }
            }
            else {
                int32_t j;
                for (j=n-1;; j--) {
                    if (X[j+l+2].s[l] == opt) {
                        break;
                    }
                }
                target_end = j;
                j--;
                for (; j >= 0; j--) {
                    if (X[j+l+2].s[l] == opt) {
                        n_best++;
                    }
                }
            }
          }
    }
    else {
        if (qsc == 0) {
            if (n < 7) {
                for (bl=0; bl < nb; bl+=8) {
                    m128si16 tmp_opt;
                    A = _mm_lddqu_si128((__m128i*)a);
                    a += 8;
                    b = B + n+7;
                      {
                        __m128i tmp = _mm_set1_epi16(oe + e * bl);
                        X[0].m = _mm_sub_epi16(tmp, mme);
                        X[1].m = tmp;
                        tmp_opt.m = mNINF;
                        __m128i h = _mm_add_epi16(tmp, mme);
                        __m128i v = mNINF;
                        __m128i maskF = mnegones;


                        for (t=0; t < n; t++) {
                            maskF = _mm_slli_si128(maskF, 2);

                            __m128i x = _mm_or_si128(_mm_slli_si128(X[t].m, 2), _mm_srli_si128(X[t+8].m, 14)); //common
                            v = _mm_insert_epi16(v, MV[t+8], 0); //common

                            __m128i mask = _mm_cmpeq_epi16(A, _mm_lddqu_si128((__m128i*)b)); //common
                            __m128i comp = _mm_or_si128(_mm_and_si128(mask, mmm), _mm_andnot_si128(mask, mmi)); //common

                            comp = _mm_or_si128(_mm_andnot_si128(maskF, comp), _mm_and_si128(maskF, mNINF));

                            const __m128i M = _mm_add_epi16(x, comp);
                            x = _mm_max_epi16(_mm_max_epi16(h, v), M); //common

                            tmp_opt.m = _mm_max_epi16(tmp_opt.m, x); //common
                            X[t+2].m = x; //common

                            const __m128i Moe = _mm_add_epi16(M, mmoe); //common
                            h = _mm_max_epi16(_mm_add_epi16(h, mme), Moe); //common
                            v = _mm_max_epi16(_mm_add_epi16(v, mme), Moe); //common

                            MV[t+1] = _mm_extract_epi16(v, 7); //common
                            v = _mm_slli_si128(v, 2); //common
                            --b; //common
                        }

                        __m128i maskL = mnegones;

                        for (t=n; t < n+7; t++) {
                            maskF = _mm_slli_si128(maskF, 2);
                            maskL = _mm_slli_si128(maskL, 2);

                            __m128i x = _mm_or_si128(_mm_slli_si128(X[t].m, 2), _mm_srli_si128(X[t+8].m, 14)); //common
                            v = _mm_insert_epi16(v, MV[t+8], 0); //common

                            __m128i mask = _mm_cmpeq_epi16(A, _mm_lddqu_si128((__m128i*)b)); //common
                            __m128i comp = _mm_or_si128(_mm_and_si128(mask, mmm), _mm_andnot_si128(mask, mmi)); //common

                            comp = _mm_or_si128(_mm_andnot_si128(maskF, comp), _mm_and_si128(maskF, mNINF));

                            const __m128i M = _mm_add_epi16(x, comp);
                            x = _mm_max_epi16(_mm_max_epi16(h, v), M); //common

                            x = _mm_or_si128(_mm_and_si128(maskL, x), _mm_andnot_si128(maskL, mNINF));
                            tmp_opt.m = _mm_max_epi16(tmp_opt.m, x); //common
                            X[t+2].m = x; //common

                            const __m128i Moe = _mm_add_epi16(M, mmoe); //common
                            h = _mm_max_epi16(_mm_add_epi16(h, mme), Moe); //common
                            v = _mm_max_epi16(_mm_add_epi16(v, mme), Moe); //common

                            MV[t+1] = _mm_extract_epi16(v, 7); //common
                            v = _mm_slli_si128(v, 2); //common
                            --b; //common
                        }
                      }

                      {
                        int32_t lmax = bl != nb-8 ? 8 : (m+7)%8+1;
                        int32_t tmp_opt_max = NINF;
                        for (l=0; l < lmax; l++)
                          tmp_opt_max = max(tmp_opt_max, tmp_opt.s[l]);

                        if (dir == 0) {
                            if (tmp_opt_max > opt) {
                                opt = tmp_opt_max;
                                n_best = 1;

                                int32_t l;
                                for (l=0;; l++) {
                                    if (tmp_opt.s[l] == opt) {
                                        int32_t j;
                                        for (j=0; X[j+l+2].s[l] != opt; j++) {
                                            ;
                                        }
                                        target_end = j;
                                        j++;
                                        for (; j < n; j++) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                        break;
                                    }
                                }
                                query_end = l+bl;
                                l++;
                                for (; l < lmax; l++) {
                                    if (tmp_opt.s[l] == opt) {
                                        for (j=0; j < n; j++) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                    }
                                }
                            }
                            else if (tmp_opt_max == opt) {
                                for (l=0; l < lmax; l++) {
                                    if (tmp_opt.s[l] == opt) {
                                        for (j=0; j < n; j++) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else {
                            if (tmp_opt_max > opt) {
                                opt = tmp_opt_max;
                                n_best = 0;
                            }

                            if (tmp_opt_max == opt) {
                                int32_t l;
                                for (l=lmax-1;; l--) {
                                    if (tmp_opt.s[l] == opt) {
                                        query_end = l+bl;
                                        int32_t j;
                                        for (j=n-1; X[j+l+2].s[l] != opt; j--) {
                                            ;
                                        }
                                        target_end = j;
                                        j--;
                                        n_best++;
                                        for (; j >=0; j--) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                        break;
                                    }
                                }
                                l--;
                                for (; l >= 0; l--) {
                                    if (tmp_opt.s[l] == opt) {
                                        for (j=n-1; j >=0; j--) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                      }
                }
            }
            else {
                for (bl=0; bl < nb; bl+=8) {
                    m128si16 tmp_opt;
                    A = _mm_lddqu_si128((__m128i*)a);
                    a += 8;
                    b = B + n+7;
                      {
                        __m128i tmp = _mm_set1_epi16(oe + e * bl);
                        X[0].m = _mm_sub_epi16(tmp, mme);
                        X[1].m = tmp;
                        tmp_opt.m = mNINF;
                        __m128i h = _mm_add_epi16(tmp, mme);
                        __m128i v = mNINF;
                        __m128i maskF = mnegones;

                        for (t=0; t < 7; t++) {
                            maskF = _mm_slli_si128(maskF, 2);

                            __m128i x = _mm_or_si128(_mm_slli_si128(X[t].m, 2), _mm_srli_si128(X[t+8].m, 14)); //common
                            v = _mm_insert_epi16(v, MV[t+8], 0); //common

                            __m128i mask = _mm_cmpeq_epi16(A, _mm_lddqu_si128((__m128i*)b)); //common
                            __m128i comp = _mm_or_si128(_mm_and_si128(mask, mmm), _mm_andnot_si128(mask, mmi)); //common

                            comp = _mm_or_si128(_mm_andnot_si128(maskF, comp), _mm_and_si128(maskF, mNINF));

                            const __m128i M = _mm_add_epi16(x, comp);
                            x = _mm_max_epi16(_mm_max_epi16(h, v), M); //common

                            tmp_opt.m = _mm_max_epi16(tmp_opt.m, x); //common
                            X[t+2].m = x; //common

                            const __m128i Moe = _mm_add_epi16(M, mmoe); //common
                            h = _mm_max_epi16(_mm_add_epi16(h, mme), Moe); //common
                            v = _mm_max_epi16(_mm_add_epi16(v, mme), Moe); //common

                            MV[t+1] = _mm_extract_epi16(v, 7); //common
                            v = _mm_slli_si128(v, 2); //common
                            --b; //common
                        }

                        for (t=7; t < n; t++) {
                            __m128i x = _mm_or_si128(_mm_slli_si128(X[t].m, 2), _mm_srli_si128(X[t+8].m, 14)); //common
                            v = _mm_insert_epi16(v, MV[t+8], 0); //common

                            __m128i mask = _mm_cmpeq_epi16(A, _mm_lddqu_si128((__m128i*)b)); //common
                            const __m128i comp = _mm_add_epi16(x, mmi);
                            mask = _mm_and_si128(mask, mmmi);

                            const __m128i M = _mm_add_epi16(mask, comp);
                            x = _mm_max_epi16(_mm_max_epi16(h, v), M); //common

                            tmp_opt.m = _mm_max_epi16(tmp_opt.m, x); //common
                            X[t+2].m = x; //common

                            const __m128i Moe = _mm_add_epi16(M, mmoe); //common
                            h = _mm_max_epi16(_mm_add_epi16(h, mme), Moe); //common
                            v = _mm_max_epi16(_mm_add_epi16(v, mme), Moe); //common

                            MV[t+1] = _mm_extract_epi16(v, 7); //common
                            v = _mm_slli_si128(v, 2); //common
                            --b; //common
                        }

                        __m128i maskL = mnegones;

                        for (t=n; t < n+7; t++) {
                            maskL = _mm_slli_si128(maskL, 2);

                            __m128i x = _mm_or_si128(_mm_slli_si128(X[t].m, 2), _mm_srli_si128(X[t+8].m, 14)); //common
                            v = _mm_insert_epi16(v, MV[t+8], 0); //common

                            __m128i mask = _mm_cmpeq_epi16(A, _mm_lddqu_si128((__m128i*)b)); //common
                            const __m128i comp = _mm_add_epi16(x, mmi);
                            mask = _mm_and_si128(mask, mmmi);

                            const __m128i M = _mm_add_epi16(mask, comp);
                            x = _mm_max_epi16(_mm_max_epi16(h, v), M); //common

                            x = _mm_or_si128(_mm_and_si128(maskL, x), _mm_andnot_si128(maskL, mNINF));
                            tmp_opt.m = _mm_max_epi16(tmp_opt.m, x); //common
                            X[t+2].m = x; //common

                            const __m128i Moe = _mm_add_epi16(M, mmoe); //common
                            h = _mm_max_epi16(_mm_add_epi16(h, mme), Moe); //common
                            v = _mm_max_epi16(_mm_add_epi16(v, mme), Moe); //common

                            MV[t+1] = _mm_extract_epi16(v, 7); //common
                            v = _mm_slli_si128(v, 2); //common
                            --b; //common
                        }
                      }

                      {
                        int32_t lmax = bl != nb-8 ? 8 : (m+7)%8+1;
                        int32_t tmp_opt_max = NINF;
                        for (l=0; l < lmax; l++)
                          tmp_opt_max = max(tmp_opt_max, tmp_opt.s[l]);

                        if (dir == 0) {
                            if (tmp_opt_max > opt) {
                                opt = tmp_opt_max;
                                n_best = 1;

                                int32_t l;
                                for (l=0;; l++) {
                                    if (tmp_opt.s[l] == opt) {
                                        int32_t j;
                                        for (j=0; X[j+l+2].s[l] != opt; j++) {
                                            ;
                                        }
                                        target_end = j;
                                        j++;
                                        for (; j < n; j++) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                        break;
                                    }
                                }
                                query_end = l+bl;
                                l++;
                                for (; l < lmax; l++) {
                                    if (tmp_opt.s[l] == opt) {
                                        for (j=0; j < n; j++) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                    }
                                }
                            }
                            else if (tmp_opt_max == opt) {
                                for (l=0; l < lmax; l++) {
                                    if (tmp_opt.s[l] == opt) {
                                        for (j=0; j < n; j++) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else {
                            if (tmp_opt_max > opt) {
                                opt = tmp_opt_max;
                                n_best = 0;
                            }

                            if (tmp_opt_max == opt) {
                                int32_t l;
                                for (l=lmax-1;; l--) {
                                    if (tmp_opt.s[l] == opt) {
                                        query_end = l+bl;
                                        int32_t j;
                                        for (j=n-1; X[j+l+2].s[l] != opt; j--) {
                                            ;
                                        }
                                        target_end = j;
                                        j--;
                                        n_best++;
                                        for (; j >=0; j--) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                        break;
                                    }
                                }
                                l--;
                                for (; l >= 0; l--) {
                                    if (tmp_opt.s[l] == opt) {
                                        for (j=n-1; j >=0; j--) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                      }
                }
            }
        }
        else {
            for (bl=0; bl < nb; bl+=8) {
                m128si16 tmp_opt;
                A = _mm_lddqu_si128((__m128i*)a);
                a += 8;
                b = B + n+7;
                  {
                    X[0].m = _mm_setzero_si128();
                    X[1].m = _mm_setzero_si128();
                    tmp_opt.m = _mm_setzero_si128();
                    __m128i h = mmoe;
                    __m128i v = mmoe;

                    for (t=0; t < n; t++) {
                        __m128i x = _mm_or_si128(_mm_slli_si128(X[t].m, 2), _mm_srli_si128(X[t+8].m, 14)); //common
                        v = _mm_insert_epi16(v, MV[t+8], 0); //common

                        __m128i mask = _mm_cmpeq_epi16(A, _mm_lddqu_si128((__m128i*)b)); //common
                        const __m128i comp = _mm_add_epi16(x, mmi);
                        mask = _mm_and_si128(mask, mmmi);

                        const __m128i M = _mm_max_epi16(_mm_setzero_si128(), _mm_add_epi16(mask, comp));
                        x = _mm_max_epi16(_mm_max_epi16(h, v), M); //common

                        X[t+2].m = x; //common
                        tmp_opt.m = _mm_max_epi16(tmp_opt.m, x); //common

                        const __m128i Moe = _mm_add_epi16(M, mmoe); //common
                        h = _mm_max_epi16(_mm_add_epi16(h, mme), Moe); //common
                        v = _mm_max_epi16(_mm_add_epi16(v, mme), Moe); //common

                        MV[t+1] = _mm_extract_epi16(v, 7); //common
                        v = _mm_slli_si128(v, 2); //common
                        --b; //common
                    }

                    __m128i maskL = mnegones;

                    for (t=n; t < n+7; t++) {
                        maskL = _mm_slli_si128(maskL, 2);

                        __m128i x = _mm_or_si128(_mm_slli_si128(X[t].m, 2), _mm_srli_si128(X[t+8].m, 14)); //common
                        v = _mm_insert_epi16(v, MV[t+8], 0); //common

                        __m128i mask = _mm_cmpeq_epi16(A, _mm_lddqu_si128((__m128i*)b)); //common
                        const __m128i comp = _mm_add_epi16(x, mmi);
                        mask = _mm_and_si128(mask, mmmi);

                        const __m128i M = _mm_max_epi16(_mm_setzero_si128(), _mm_add_epi16(mask, comp));
                        x = _mm_max_epi16(_mm_max_epi16(h, v), M); //common

                        X[t+2].m = x; //common
                        x = _mm_and_si128(x, maskL);
                        tmp_opt.m = _mm_max_epi16(tmp_opt.m, x); //common

                        const __m128i Moe = _mm_add_epi16(M, mmoe); //common
                        h = _mm_max_epi16(_mm_add_epi16(h, mme), Moe); //common
                        v = _mm_max_epi16(_mm_add_epi16(v, mme), Moe); //common

                        MV[t+1] = _mm_extract_epi16(v, 7); //common
                        v = _mm_slli_si128(v, 2); //common
                        --b; //common
                    }
                  }

                  {
                    int32_t lmax = bl != nb-8 ? 8 : (m+7)%8+1;
                    int32_t tmp_opt_max = NINF;
                    for (l=0; l < lmax; l++)
                      tmp_opt_max = max(tmp_opt_max, tmp_opt.s[l]);

                    if (dir == 0) {
                        if (tmp_opt_max > opt) {
                            opt = tmp_opt_max;
                            n_best = 1;

                            int32_t l;
                            for (l=0;; l++) {
                                if (tmp_opt.s[l] == opt) {
                                    int32_t j;
                                    for (j=0; X[j+l+2].s[l] != opt; j++) {
                                        ;
                                    }
                                    target_end = j;
                                    j++;
                                    for (; j < n; j++) {
                                        if (X[j+l+2].s[l] == opt) {
                                            n_best++;
                                        }
                                    }
                                    break;
                                }
                            }
                            query_end = l+bl;
                            l++;
                            for (; l < lmax; l++) {
                                if (tmp_opt.s[l] == opt) {
                                    for (j=0; j < n; j++) {
                                        if (X[j+l+2].s[l] == opt) {
                                            n_best++;
                                        }
                                    }
                                }
                            }
                        }
                        else if (tmp_opt_max == opt) {
                            for (l=0; l < lmax; l++) {
                                if (tmp_opt.s[l] == opt) {
                                    for (j=0; j < n; j++) {
                                        if (X[j+l+2].s[l] == opt) {
                                            n_best++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else {
                        if (tmp_opt_max > opt) {
                            opt = tmp_opt_max;
                            n_best = 0;
                        }

                        if (tmp_opt_max == opt) {
                            int32_t l;
                            for (l=lmax-1;; l--) {
                                if (tmp_opt.s[l] == opt) {
                                    query_end = l+bl;
                                    int32_t j;
                                    for (j=n-1; X[j+l+2].s[l] != opt; j--) {
                                        ;
                                    }
                                    target_end = j;
                                    j--;
                                    n_best++;
                                    for (; j >=0; j--) {
                                        if (X[j+l+2].s[l] == opt) {
                                            n_best++;
                                        }
                                    }
                                    break;
                                }
                            }
                            l--;
                            for (; l >= 0; l--) {
                                if (tmp_opt.s[l] == opt) {
                                    for (j=n-1; j >=0; j--) {
                                        if (X[j+l+2].s[l] == opt) {
                                            n_best++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                  }
            }
        }
    }

    result->opt = opt;
    result->query_end = query_end;
    result->target_end = target_end;
    result->n_best = n_best;
}
static tmap_vsw_data_s3_t*
tmap_vsw_data_init_s3_helper(const uint8_t *query, int32_t qlen, int32_t tlen, int32_t query_start_clip, int32_t query_end_clip, tmap_vsw_opt_t *opt)
{
  tmap_vsw_data_s3_t *vsw = NULL;
  vsw = tmap_calloc(1, sizeof(tmap_vsw_data_s3_t), "vsw");
  vsw->mem_qlen = qlen;
  vsw->mem_tlen = tlen;
  tmap_roundup32(vsw->mem_qlen);
  tmap_roundup32(vsw->mem_tlen);
  vsw->max_qlen = INT32_MAX;
  vsw->max_tlen = INT32_MAX;
  vsw->query_start_clip = query_start_clip;
  vsw->query_end_clip = query_end_clip;
  vsw->abuf = tmap_calloc(vsw->mem_qlen, sizeof(int16_t), "vsw->abuf");
  vsw->B = tmap_calloc(vsw->mem_tlen+16, sizeof(int16_t), "vsw->B");
  vsw->MV = tmap_calloc(vsw->mem_tlen+16, sizeof(int16_t), "vsw->MV");
  vsw->X = tmap_calloc(vsw->mem_tlen+16, sizeof(m128si16), "vsw->X");
  return vsw;
}

tmap_vsw_data_s3_t*
tmap_vsw_data_init_s3(const uint8_t *query, int32_t qlen, int32_t query_start_clip, int32_t query_end_clip, tmap_vsw_opt_t *opt)
{
  return tmap_vsw_data_init_s3_helper(query, qlen, qlen, query_start_clip, query_end_clip, opt);
}

tmap_vsw_data_s3_t*
tmap_vsw_data_update_s3(tmap_vsw_data_s3_t *vsw, const uint8_t *query, int32_t qlen, const uint8_t *target, int32_t tlen)
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
      tmap_vsw_data_destroy_s3(vsw); // destroy
      vsw = tmap_vsw_data_init_s3_helper(query, mem_qlen, mem_tlen, qsc, qec, opt);
  }
  return vsw;
}

void
tmap_vsw_data_destroy_s3(tmap_vsw_data_s3_t *vsw)
{
  if(NULL == vsw) return;
  free(vsw->abuf);
  free(vsw->B);
  free(vsw->MV);
  free(vsw->X);
  free(vsw);
}

int32_t
tmap_vsw_process_s3(tmap_vsw_data_s3_t *vsw,
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
