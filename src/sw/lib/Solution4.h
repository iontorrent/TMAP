#ifndef SOLUTION4_H_
#define SOLUTION4_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

#include <cstring>
#include <sstream>
#include "Solution.h"

using namespace std;

typedef struct {
    uint32_t hash;
    short opt;
    short n_best;
    int res_min_pos;
    int res_max_pos;
} qres_t;

class Solution4 : public Solution {
public:
  Solution4();
  ~Solution4();

  virtual int process(string &b, string &a, int qsc, int qec, 
                      int mm, int mi, int o, int e, int dir,
                      int *opt, int *te, int *qe, int *n_best);
private:
  int MAX_DIMA;
  int MAX_DIMB;
  int8_t DNA_CONV[128];
  //int16_t BUFFER[MAX_DIMB * 9] __attribute__((aligned(64)));
  int16_t *BUFFER __attribute__((aligned(64)));
  int n, m;
  int segNo;
  int len;
  __m128i* XP[4];
  __m128i* M0;
  __m128i* M1;
  __m128i* V;
  int16_t* POS;
  int16_t INVALID_POS[16];
  int INVALID_POS_NO;
  int lastMax;
  int opt, n_best, res_min_pos, res_max_pos;
  qres_t *HTDATA;
  int64_t count0, count1;

  // Utility functions
  int16_t findMax16(const __m128i &mMax);
  int16_t findMax16Simple(const __m128i &mMax);
  uint8_t findMax8(const __m128i &mMax);
  uint64_t hashDNA(const string &s, const int len);
  int resize(int a, int b);
  template <class T, bool BYTE> void updateResult(int i, int curMax);
  template <class T, int DEFAULT_VALUE> void updateResultLast();
  template <int qec> void processFastVariantB16BitA(string &a, int mm, int mi, int o, int e);
  template <int qec> void processFastVariantB16BitB(string &a, int mm, int mi, int o, int e);
  template <int qec> void processFastVariantA16Bit(string &a, int mm, int mi, int o, int e, int iter);
  template <int qec> int processFastVariantA8Bit(const string &a, const int mm, const int mi, const int o, const int e);
  void convertTable(__m128i *T0, __m128i *T1);
  template <class T, int SIZE> void calcInvalidPos(int value);
  void convert16Bit(string &b, int qsc, int mm, int mi);
  void preprocess16Bit(string &b, int qsc, int mm, int mi);
  void preprocess8Bit(string &b, int qsc, int mm, int mi);
};

#endif
