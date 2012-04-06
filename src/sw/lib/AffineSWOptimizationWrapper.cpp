#include "AffineSWOptimization.h"
#include "AffineSWOptimizationWrapper.h"
#include <iostream>

using namespace std;

tmap_vsw_wrapper_t*
tmap_vsw_wrapper_init(int32_t type)
{
  return new AffineSWOptimization(type);
}

int32_t 
tmap_vsw_wrapper_process(tmap_vsw_wrapper_t *v,
                         const uint8_t *target, int32_t tlen,
                         const uint8_t *query, int32_t qlen,
                         int mm, int mi, int o, int e, int dir,
                         int qsc, int qec,
                         int *opt, int *te, int *qe, int *n_best)
{
  int i;
  string a, b;

  for(i=0;i<tlen;i++) b += target[i];
  for(i=0;i<qlen;i++) a += query[i];
  
  /*
  fprintf(stderr, "%s tlen=%d qlen=%d\n", __func__, tlen, qlen);
  fprintf(stderr, "%s mm=%d mi=%d o=%d e=%d\n", __func__, mm, mi, o, e);
  fprintf(stderr, "%s dir=%d qsc=%d qec=%d\n", __func__, dir, qsc, qec);
  cerr << "b: " << b << endl;
  cerr << "a: " << a << endl;
  */

  return v->process(b, a, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best);
}

void
tmap_vsw_wrapper_destroy(tmap_vsw_wrapper_t *v)
{
  delete v;
}

