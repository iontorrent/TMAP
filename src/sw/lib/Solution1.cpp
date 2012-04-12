#include <cstring>
#include <sstream>
#include <stdint.h>
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_definitions.h"
#include "vsw16.h"
#include "vsw.h"
#include "vsw16.h"
#include "Solution1.h"

// TMAP Vectorized Smith Waterman

using namespace std;

#define SCORE_THR 0

// Input: ASCII character
// Output: 2-bit DNA value
static uint8_t nt_char_to_int[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

Solution1::Solution1() {
    // dummy ctor
    vsw_opt = NULL;
    vsw_query = NULL;
    target = query = NULL;
    target_len = query_len = 0;
    q_max = t_max = 0;
    max_qlen = 1 << 31;
    max_tlen = 1 << 31;
}

Solution1::~Solution1() {
    free(query);
    free(target);
    vsw_opt_destroy(vsw_opt);
    vsw_opt = NULL;
    if(NULL != vsw_query) vsw_query_destroy(vsw_query);
    vsw_query = NULL;
}

int Solution1::process(const string& b, const string& a, int qsc, int qec,
                                     int mm, int mi, int o, int e, int dir,
                                     int *_opt, int *_te, int *_qe, int *_n_best) {
    int i, score;
    int n = b.size(), m = a.size();
    int16_t query_end, target_end;
    int32_t overflow, n_best;
    // set opt
    if(NULL == vsw_opt) vsw_opt = vsw_opt_init(mm, -mi, -o, -e, SCORE_THR); // NB: mm/mi/o/e are fixed
    // set query
    if(0 == query_cmp(a, m, query, query_len)) {
        query_len = m;
        if(q_max < query_len) {
            q_max = query_len;
            tmap_roundup32(q_max);
            query = (uint8_t*)tmap_realloc(query, sizeof(uint8_t) * q_max, "query");
        }
        for(i=0;i<m;i++) {
            query[i] = nt_char_to_int[(int)a[i]];
        }
        if(NULL != vsw_query) vsw_query_destroy(vsw_query);
        vsw_query = vsw_query_init(query, query_len, query_len, qsc, qec, vsw_opt); 
    }

    // reset query
    vsw_query->query16 = vsw16_query_init(vsw_query->query16, query, query_len, qsc, qec, vsw_opt);

    // set target
    target_len = n;
    if(t_max < target_len) {
        t_max = target_len;
        tmap_roundup32(t_max);
        target = (uint8_t*)tmap_realloc(target, sizeof(uint8_t) * t_max, "target");
    }
    for(i=0;i<(int)target_len;i++) {
        target[i] = nt_char_to_int[(int)b[i]];
    }

    // run SW
    overflow = 0;
    score = vsw16_sse2_forward(vsw_query->query16, target, target_len,
                               qsc, qec,
                               vsw_opt, &query_end, &target_end,
                               dir, &overflow, &n_best, SCORE_THR);

    // return results
    (*_opt) = score;
    (*_te) = target_end;
    (*_qe) = query_end;
    (*_n_best) = n_best;
    //fprintf(stderr, "%s SCORE=%d %d %d %d %d\n", __func__, score, target_end, query_end, n_best, overflow);
    return score;
}

/*
   vsw_opt_t *vsw_opt = NULL;
   vsw_query_t* vsw_query = NULL;
   uint8_t query[Q_MAX];
   uint32_t query_len = 0;
   uint8_t target[T_MAX];
   uint32_t target_len = 0;
   */

int32_t Solution1::query_cmp(string q, int n, uint8_t *query, int32_t query_len) {
    int i;
    if(NULL == query) return 0;
    if(query_len != n) return 0;
    for(i=0;i<n;i++) {
        if(q[i] != "ACGTN"[query[i]]) return 0;
    }
    return 1;
}
