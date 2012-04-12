#ifndef AFFINESWOPTIMIZATIONHASH_H_
#define AFFINESWOPTIMIZATIONHASH_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

#include <cstring>
#include <sstream>

using namespace std;

typedef struct {
    unsigned long long int hash;
    int tlen;
    int dir;
    string b;
    int opt;
    int te;
    int qe;
    int n_best;
} hash_t;

class AffineSWOptimizationHash {
public:
  AffineSWOptimizationHash();

  bool process(string b, string a, int qsc, int qec,
                  int mm, int mi, int o, int e, int dir,
                  int *opt, int *te, int *qe, int *n_best);

  void add(string b, string a, int _qsc, int _qec,
           int mm, int mi, int o, int e, int dir,
           int *opt, int *te, int *qe, int *n_best);

  ~AffineSWOptimizationHash();
private:
  int qsc;
  int qec;
  int size;
  string query;
  hash_t *hash;
};
#endif
