#ifndef SOLUTION5_H_
#define SOLUTION5_H_

#ifndef __cplusplus
#error a C++ compiler is required
#endif

#include <cstring>
#include <sstream>
#include "solution.h"

using namespace std;
  
typedef struct { 
  int idx,f,h;
} dpstate_t;

typedef struct {
    int f, h;
} dpstate_large_t;

class Solution5 : public Solution {
public:
  Solution5();
  ~Solution5();

  virtual int process(string& b, string& a, int qsc, int qec,
                      int mm, int mi, int o, int e, int dir,
                      int *opt, int *te, int *qe, int *n_best);
private:

  int opt, query_end,target_end,n_best;
  int MAX_LENGTH;

  int m,n;
  char *target, *query;
  int *a, *b, *ca, *cb;
  int task_id;
  int match_score, mismatch_score, gap_open, gap_extension;
  int direction;
  uint32_t batch8[1024];
  uint32_t mismatch_score_batch16;

  int sizeq[2];
  dpstate_t *dp_buffer;

  int *q[3];
  dpstate_large_t **dp_large;
  int counter_task2;
  bool **vf, **vh;

  short *hbuffer16 __attribute__ ((aligned (16)));
  uint8_t *hbuffer8 __attribute__ ((aligned (16)));

  int int_buffer[4] __attribute__ ((aligned (16)));
  short short_buffer[8] __attribute__ ((aligned (16)));
  uint8_t uchar_buffer[16] __attribute__ ((aligned (16)));

  int vdeg[256], *vgraph[256];
  int vfirst[256], *vnext;

  int testcase_id;
  int measure_id;

  int est_opt;
  int first_j, last_j;

  bool valid0[256];
  int encode_idx[128];

  int *shuffle;
  uint8_t *cost8[4] __attribute__ ((aligned (16)));
  short *cost16[4] __attribute__ ((aligned (16)));

  short *transfer_h, *transfer_e;

  int *best_end_points;

  // Utility functions
  void resize(int m);
  void init();
  void brute_force();
  void solve_task1();
  void solve_task2();
  void solve_task3();
  void solve();
  void shift_right16(__m128i &a, short first);
  void shift_right8(__m128i &a, uint8_t first);


};

#endif
