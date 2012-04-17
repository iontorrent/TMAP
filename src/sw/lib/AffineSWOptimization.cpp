#include <cstring>
#include <sstream>
#include <iostream>
#include "Solution1.h"
#include "Solution2.h"
#include "Solution3.h"
#include "Solution4.h"
#include "Solution5.h"
#include "Solution6.h"
#include "Solution7.h"
#include "Solution8.h"
#include "Solution9.h"
#include "Solution10.h"
#include "AffineSWOptimization.h"

using namespace std;

AffineSWOptimization::AffineSWOptimization(int type) {
    myType = type;
    switch(type) {
      case 1:
        s = new Solution1();
        break;
      case 2:
        s = new Solution2();
        break;
      case 3:
        s = new Solution3();
        break;
      case 4:
        s = new Solution4();
        break;
      case 5:
        s = new Solution5();
        break;
      case 6:
        s = new Solution6();
        break;
      case 7:
        s = new Solution7();
        break;
      case 8:
        s = new Solution8();
        break;
      case 9:
        s = new Solution9();
        break;
      case 10:
        s = new Solution10();
        break;
      default:
        fprintf(stderr, "Error: unknown type: %d\n", type);
        exit(1);
    }
#ifdef AFFINESWOPTIMIZATION_USE_HASH
    hash = new AffineSWOptimizationHash();
#else
    hash = NULL;
#endif
    a.reserve(512);
    b.reserve(1024);
}

int AffineSWOptimization::process(const uint8_t *target, int32_t tlen,
                                  const uint8_t *query, int32_t qlen,
                                  int qsc, int qec,
                                  int mm, int mi, int o, int e, int dir,
                                  int *opt, int *te, int *qe, int *n_best) {
    int i;
    // resize
    b.resize(tlen);
    a.resize(qlen);
    // copy
    for(i=0;i<tlen;i++) b[i] = "ACGTN"[target[i]];
    for(i=0;i<qlen;i++) a[i] = "ACGTN"[query[i]];
#ifdef AFFINESWOPTIMIZATION_USE_HASH
    // try the hash
    if(!hash->process(b, a, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best)) {
        s->process(b, a, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best);
        // add to the hash
        hash->add(b, a, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best);
    }
    return (*opt);
#else
    return s->process(b, a, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best);
#endif
}

AffineSWOptimization::~AffineSWOptimization()
{
#ifdef AFFINESWOPTIMIZATION_USE_HASH
  delete hash;
#endif
  delete s;
}

