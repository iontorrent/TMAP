#ifndef SOLUTION9_H_
#define SOLUTION9_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

#include <cstring>
#include <sstream>
#include "Solution.h"

using namespace std;

class Solution9Reference {
    const static int MAX_DIM = 1025;
    const static int INF = 1073741824;

public:
    int M[MAX_DIM][MAX_DIM];
    int H[MAX_DIM][MAX_DIM];
    int V[MAX_DIM][MAX_DIM];

    int process(string b, string a, int qsc, int qec, int mm, int mi, int o, int e, int dir,
                         int *opt, int *te, int *qe, int *n_best);
};

class Solution9 : public Solution {

public:
  Solution9Reference Reference;
  int MatchScore, MismatchScore, GapOpen, GapExtension;
  int Opt, QueryEnd, TargetEnd, NBest;
  int QuerySize, TargetSize, TargetSizeInSIMDBlocks;

  Solution9();
  ~Solution9();

  virtual int process(string& b, string& a, int qsc, int qec,
                         int mm, int mi, int o, int e, int dir,
                         int *opt, int *te, int *qe, int *n_best);
private:
  template<bool DIRECTION>
    void inline CheckOpt8(const unsigned char m, const int i, const int j);
  template<bool DIRECTION>
    void inline CheckOpt16(const short int m, const int i, const int j);
  template<bool DIRECTION>
    void ExtractOpts8(const __m128i newM, const int ltMask, const int i, const int jStart);
  template<bool DIRECTION>
    void ExtractOpts16(const __m128i newM, const int ltMask, const int i, const int jStart);
  template<int VARIANT, bool DIRECTION, bool DETECTREDUNDANT> void DoVariant8(int threshold);
  template<int VARIANT, bool DIRECTION, bool DETECTREDUNDANT> void DoVariant16();
  void DoSetup8(const string& target, const string& query);
  void DoSetup16(const string& target, const string& query);
};

#endif
