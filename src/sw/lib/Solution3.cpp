#include <cstring>
#include <sstream>
#include "sw-vector.h"
#include "solution3.h"

// SHRiMP2's Vectorized Smith Waterman

using namespace std;

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

Solution3::Solution3() {
    n = 0;
    // dummy ctor
}

int Solution3::process(string& b, string& a, int qsc, int qec,
                          int mm, int mi, int o, int e, int dir,
                          int *_opt, int *_te, int *_qe, int *_n_best) {
    int i, bl, al;
    int opt, n_best, query_end, target_end;

    bl = b.length();
    al = a.length();
    if(0 == n) {
        sw_vector_setup(T_MAX, Q_MAX, o, e, o, e, mm, mi, 0, true);
    }
    // target
    for(i=0;i<al;i++) {
        target[i] = nt_char_to_int[(int)a[i]];
    }
    // query
    for(i=0;i<bl;i++) {
        query[i] = nt_char_to_int[(int)b[i]]; 
    }
    
    opt = n_best = query_end = target_end = 0;

    opt = sw_vector(target, bl, query, al);

    (*_opt) = opt;
    (*_te) = target_end;
    (*_qe) = query_end;
    (*_n_best) = n_best;

    return opt;
}
