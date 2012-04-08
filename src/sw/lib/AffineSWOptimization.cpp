#include <cstring>
#include <sstream>
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
    hash = new AffineSWOptimizationHash();
}

int AffineSWOptimization::process(string b, string a, int qsc, int qec,
                                     int mm, int mi, int o, int e, int dir,
                                     int *opt, int *te, int *qe, int *n_best) {
    // try the hash
    if(!hash->process(b, a, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best)) {
        s->process(b, a, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best);
        // add to the hash
        hash->add(b, a, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best);
    }
    return (*opt);
}

AffineSWOptimization::~AffineSWOptimization()
{
  //s->~Solution();
  delete hash;
  delete s;
}

