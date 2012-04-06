#include <cstring>
#include <sstream>
#include "solution2.h"

// Ivan's solution

using namespace std;

void Solution2::resize(int n)
{
  int i;
  if(n <= mem) return;
  M = (int**)realloc(M, sizeof(int*) * n);
  H = (int**)realloc(H, sizeof(int*) * n);
  V = (int**)realloc(V, sizeof(int*) * n);
  for(i=0;i<mem;i++) {
      M[i] = (int*)realloc(M[i], sizeof(int) * n);
      H[i] = (int*)realloc(H[i], sizeof(int) * n);
      V[i] = (int*)realloc(V[i], sizeof(int) * n);
  }
  for(;i<n;i++) {
      M[i] = (int*)malloc(sizeof(int) * n);
      H[i] = (int*)malloc(sizeof(int) * n);
      V[i] = (int*)malloc(sizeof(int) * n);
  }
  mem = n;
}

Solution2::Solution2()
{
  int i;
  mem = MAX_DIM;
  M = (int**)malloc(sizeof(int*) * MAX_DIM);
  H = (int**)malloc(sizeof(int*) * MAX_DIM);
  V = (int**)malloc(sizeof(int*) * MAX_DIM);
  for(i=0;i<MAX_DIM;i++) {
      M[i] = (int*)malloc(sizeof(int) * MAX_DIM);
      H[i] = (int*)malloc(sizeof(int) * MAX_DIM);
      V[i] = (int*)malloc(sizeof(int) * MAX_DIM);
  }
}

int Solution2::process(string& b, string& a, int qsc, int qec,
                                     int mm, int mi, int o, int e, int dir,
                                     int *_opt, int *_te, int *_qe, int *_n_best) {
    int n = b.size(), m = a.size();

    int id = -1;
    if (qsc == 1 && qec == 1) id = 1;
    if (qsc == 0 && qec == 1) id = 2;
    if (qsc == 1 && qec == 0) id = 3;
    if (qsc == 0 && qec == 0) id = 4;

    resize(n);
    resize(m);

    if (id == 1 || id == 3) {
        for (int i=0; i <= m; i++) M[i][0] = 0;
        for (int j=0; j <= n; j++) M[0][j] = 0;
        for (int i=0; i <= m; i++) H[i][0] = V[i][0] = -INF;
        for (int j=0; j <= n; j++) H[0][j] = V[0][j] = -INF;

        for (int i=1; i <= m; i++)
          for (int j=1; j <= n; j++) {
              V[i][j] = max(M[i-1][j] + o + e, V[i-1][j] + e);
              H[i][j] = max(M[i][j-1] + o + e, H[i][j-1] + e);
              int mx = max(max(M[i-1][j-1], V[i-1][j-1]), H[i-1][j-1]);
              M[i][j] = max(0, mx + (a[i-1] == b[j-1] ? mm : mi));
          }
    } else {                                      
        for (int j=0; j <= n; j++) M[0][j] = 0;
        for (int i=1; i <= m; i++) M[i][0] = -INF;
        for (int j=0; j <= n; j++) H[0][j] = V[0][j] = -INF;
        for (int i=1; i <= m; i++) {
            V[i][0] = -INF;
            H[i][0] = o + e * i;
        }

        for (int i=1; i <= m; i++)
          for (int j=1; j <= n; j++) {
              V[i][j] = max(M[i-1][j] + o + e, V[i-1][j] + e);
              H[i][j] = max(M[i][j-1] + o + e, H[i][j-1] + e);
              int mx = max(max(M[i-1][j-1], V[i-1][j-1]), H[i-1][j-1]);
              M[i][j] = mx + (a[i-1] == b[j-1] ? mm : mi);
          }
    }

    int minI = 1, maxI = m, minJ = 1, maxJ = n;
    if (id == 3 || id == 4) minI = m;

    int opt = -INF, query_end = -1, target_end = -1, n_best = 0;

    for (int i=minI; i <= maxI; i++)
      for (int j=minJ; j <= maxJ; j++) {
          opt = max(opt, M[i][j]);
          opt = max(opt, V[i][j]);
          opt = max(opt, H[i][j]);
      }

    n_best = 0;
    for (int i=minI; i <= maxI; i++)
      for (int j=minJ; j <= maxJ; j++)
        if (M[i][j] == opt || V[i][j] == opt || H[i][j] == opt)
          n_best++;

    if (dir == 0) {
        for (int i=minI; i <= maxI && query_end == -1; i++)
          for (int j=minJ; j <= maxJ && query_end == -1; j++)
            if (M[i][j] == opt || V[i][j] == opt || H[i][j] == opt) {
                query_end = i-1;
                target_end = j-1;
            }
    } else {
        for (int i=maxI; i >= minI && query_end == -1; i--)
          for (int j=maxJ; j >= minJ && query_end == -1; j--)
            if (M[i][j] == opt || V[i][j] == opt || H[i][j] == opt) {
                query_end = i-1;
                target_end = j-1;
            }
    }

    (*_opt) = opt;
    (*_te) = target_end;
    (*_qe) = query_end;
    (*_n_best) = n_best;

    return opt;
}
