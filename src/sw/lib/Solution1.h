#ifndef SOLUTION1_H_
#define SOLUTION1_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

#include <cstring>
#include <sstream>
#include "solution.h"
#include "vsw16.h"
#include "vsw.h"

using namespace std;

#define Q_MAX 512
#define T_MAX 1024

class Solution1 : public Solution {
public:
  Solution1();

  ~Solution1();

  virtual int process(string& b, string& a, int qsc, int qec,
                 int mm, int mi, int o, int e, int dir,
                 int *opt, int *te, int *qe, int *n_best);

  int32_t query_cmp(string q, int n, uint8_t *query, int32_t query_len);
private:
  vsw_opt_t *vsw_opt;
  vsw_query_t* vsw_query;
  uint8_t query[Q_MAX];
  int32_t query_len;
  uint8_t target[T_MAX];
  int32_t target_len;
};
#endif
