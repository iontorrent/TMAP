#ifndef SOLUTION8_H_
#define SOLUTION8_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

#include <cstring>
#include <sstream>
#include "solution.h"

using namespace std;

class Solution8 : public Solution {
public:
  Solution8();
  ~Solution8();

  virtual int process(string& b, string& a, int qsc, int qec,
                 int mm, int mi, int o, int e, int dir,
                 int *opt, int *te, int *qe, int *n_best);
private:
  bool printed;
};

#endif
