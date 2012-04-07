#ifndef SOLUTION6_H_
#define SOLUTION6_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

#include <cstring>
#include <sstream>
#include "Solution.h"

using namespace std;

class Solution6 : public Solution {
    typedef union { short s[8]; __m128i m; } m128si16;
public:
    Solution6();

    virtual int process(string& bs, string& as, int qsc, int qec,
                        int mm, int mi, int oe, int e, int dir,
                        int *opt, int *te, int *qe, int *n_best);
};

#endif
