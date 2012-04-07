#ifndef SOLUTION_H_
#define SOLUTION_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

#include <cstring>
#include <sstream>

using namespace std;

class Solution {
public:
  Solution();

  virtual int process(string& b, string& a, int qsc, int qec,
                 int mm, int mi, int o, int e, int dir,
                 int *opt, int *te, int *qe, int *n_best) = 0;

  virtual ~Solution();
};

#endif
