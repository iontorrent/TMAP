#ifndef AFFINESWOPTIMIZATION_H_
#define AFFINESWOPTIMIZATION_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

#include <cstring>
#include <sstream>
#include "Solution.h"

using namespace std;

class AffineSWOptimization {
public:
  AffineSWOptimization(int type);

  int32_t process(string b, string a, int qsc, int qec,
                 int mm, int mi, int o, int e, int dir,
                 int *opt, int *te, int *ts, int *n_best);

  ~AffineSWOptimization();
private:
  int myType;
  Solution *s;
};
#endif
