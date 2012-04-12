// Coder: folsena
// Submission: 18
// URL: http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=22926110&subnum=18
#include <cstring>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <mmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include "Solution6.h"

using namespace std;

const int INF = 1073741824;
const short NINF = -30000;

inline int max(int a, int b) {
    return a > b ? a : b;
}

Solution6::Solution6() {
  max_qlen = 512;
  max_tlen = 1024;
}

int Solution6::process(const string& bs, const string& as, int qsc, int qec,
                          int mm, int mi, int oe, int e, int dir,
                          int *_opt, int *_te, int *_qe, int *_n_best) {
    const int n = bs.size(), m = as.size();
    int opt = NINF, query_end=-1, target_end=-1, n_best=0;
    oe += e;

    short abuf[512];
    short B[1024+16];
    short MV[1024+16];
    m128si16 X[1024+16];

    fill(B, B+8, 0);
    fill_n(reverse_copy(bs.begin(), bs.end(), B+8), 8, 0);
    /*
    for(int i=0;i<n+8;i++) {
        fprintf(stderr, "i=%d B[n]=%d\n", i, B[i]);
    } 
    */
    fill(MV+8, MV+8+n, oe);
    copy(as.c_str(), as.c_str() + m, abuf);

    for (int i=0; i<n; i++)
      X[i+8].m = _mm_setzero_si128();

    const __m128i mme = _mm_set1_epi16(e);
    const __m128i mmm = _mm_set1_epi16(mm);
    const __m128i mmmi = _mm_set1_epi16(mm-mi);
    const __m128i mmi = _mm_set1_epi16(mi);

    const __m128i mmoe = _mm_set1_epi16(oe);
    const __m128i mNINF = _mm_set1_epi16(NINF);
    const __m128i mnegones = _mm_set1_epi16(-1);

    const int nb = (m+7)/8*8;
    const short* b;
    const short* a = abuf;
    __m128i A;

    if (qec == 0) {
        if (qsc == 0) {
            for (int bl=0; bl < nb; bl+=8) {
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

                    for (int t=0; t < 7; t++) {
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

                    for (int t=7; t < n+7; t++) {
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
            for (int bl=0; bl < nb; bl+=8) {
                A = _mm_lddqu_si128((__m128i*)a);
                a += 8;
                b = B + n+7;
                  {
                    X[0].m = _mm_setzero_si128();
                    X[1].m = _mm_setzero_si128();
                    __m128i h = mmoe;
                    __m128i v = mmoe;

                    for (int t=0; t < n+7; t++) {
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
            const int l = (m-1)%8;
            opt = NINF;
            n_best = 1;
            query_end = m-1;

            for (int j=0; j < n; j++)
              opt = max(opt, X[j+l+2].s[l]);

            if (dir == 0) {
                int j;
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
                int j;
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
                for (int bl=0; bl < nb; bl+=8) {
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

                        for (int t=0; t < n; t++) {
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

                        for (int t=n; t < n+7; t++) {
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
                        int lmax = bl != nb-8 ? 8 : (m+7)%8+1;
                        int tmp_opt_max = NINF;
                        for (int l=0; l < lmax; l++)
                          tmp_opt_max = max(tmp_opt_max, tmp_opt.s[l]);

                        if (dir == 0) {
                            if (tmp_opt_max > opt) {
                                opt = tmp_opt_max;
                                n_best = 1;

                                int l;
                                for (l=0;; l++) {
                                    if (tmp_opt.s[l] == opt) {
                                        int j;
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
                                        for (int j=0; j < n; j++) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                    }
                                }
                            }
                            else if (tmp_opt_max == opt) {
                                for (int l=0; l < lmax; l++) {
                                    if (tmp_opt.s[l] == opt) {
                                        for (int j=0; j < n; j++) {
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
                                int l;
                                for (l=lmax-1;; l--) {
                                    if (tmp_opt.s[l] == opt) {
                                        query_end = l+bl;
                                        int j;
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
                                        for (int j=n-1; j >=0; j--) {
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
                for (int bl=0; bl < nb; bl+=8) {
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

                        for (int t=0; t < 7; t++) {
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

                        for (int t=7; t < n; t++) {
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

                        for (int t=n; t < n+7; t++) {
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
                        int lmax = bl != nb-8 ? 8 : (m+7)%8+1;
                        int tmp_opt_max = NINF;
                        for (int l=0; l < lmax; l++)
                          tmp_opt_max = max(tmp_opt_max, tmp_opt.s[l]);

                        if (dir == 0) {
                            if (tmp_opt_max > opt) {
                                opt = tmp_opt_max;
                                n_best = 1;

                                int l;
                                for (l=0;; l++) {
                                    if (tmp_opt.s[l] == opt) {
                                        int j;
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
                                        for (int j=0; j < n; j++) {
                                            if (X[j+l+2].s[l] == opt) {
                                                n_best++;
                                            }
                                        }
                                    }
                                }
                            }
                            else if (tmp_opt_max == opt) {
                                for (int l=0; l < lmax; l++) {
                                    if (tmp_opt.s[l] == opt) {
                                        for (int j=0; j < n; j++) {
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
                                int l;
                                for (l=lmax-1;; l--) {
                                    if (tmp_opt.s[l] == opt) {
                                        query_end = l+bl;
                                        int j;
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
                                        for (int j=n-1; j >=0; j--) {
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
            for (int bl=0; bl < nb; bl+=8) {
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

                    for (int t=0; t < n; t++) {
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

                    for (int t=n; t < n+7; t++) {
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
                    int lmax = bl != nb-8 ? 8 : (m+7)%8+1;
                    int tmp_opt_max = NINF;
                    for (int l=0; l < lmax; l++)
                      tmp_opt_max = max(tmp_opt_max, tmp_opt.s[l]);

                    if (dir == 0) {
                        if (tmp_opt_max > opt) {
                            opt = tmp_opt_max;
                            n_best = 1;

                            int l;
                            for (l=0;; l++) {
                                if (tmp_opt.s[l] == opt) {
                                    int j;
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
                                    for (int j=0; j < n; j++) {
                                        if (X[j+l+2].s[l] == opt) {
                                            n_best++;
                                        }
                                    }
                                }
                            }
                        }
                        else if (tmp_opt_max == opt) {
                            for (int l=0; l < lmax; l++) {
                                if (tmp_opt.s[l] == opt) {
                                    for (int j=0; j < n; j++) {
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
                            int l;
                            for (l=lmax-1;; l--) {
                                if (tmp_opt.s[l] == opt) {
                                    query_end = l+bl;
                                    int j;
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
                                    for (int j=n-1; j >=0; j--) {
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

    //char buf[32];
    //sprintf(buf, "%d %d %d %d", opt, query_end, target_end, n_best);
    //return string(buf);
    (*_opt) = opt;
    (*_te) = target_end;
    (*_qe) = query_end;
    (*_n_best) = n_best;
    return opt;
}
