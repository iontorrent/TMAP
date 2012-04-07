#ifndef SOLUTION2_H_
#define SOLUTION2_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

#include <cstring>
#include <sstream>
#include "solution.h"

using namespace std;

const int INF = 1073741824;
const int MAX_DIM = 1025;

#define max(a, b) ((a)>(b)?a:b)

class Solution2 : public Solution {
public:
  Solution2();
  ~Solution2();

  virtual int process(string& b, string& a, int qsc, int qec,
                 int mm, int mi, int o, int e, int dir,
                 int *opt, int *te, int *qe, int *n_best);
private:
  int **M;
  int **H;
  int **V;
  int mem;
  void resize(int n);
};

#endif
