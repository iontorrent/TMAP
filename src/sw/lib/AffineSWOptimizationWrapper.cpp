#define __STDC_LIMIT_MACROS // Seriously, I want these, now
#include <stdint.h>
#include <limits>
#include <iostream>
#include "AffineSWOptimization.h"
#include "AffineSWOptimizationWrapper.h"

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

  for(i=0;i<tlen;i++) b += "ACGTN"[target[i]];
  for(i=0;i<qlen;i++) a += "ACGTN"[query[i]];
  
  /*
  // Top coder style
  for(i=0;i<tlen;i++) fputc("ACGTN"[target[i]], stderr);
  fputc('\t', stderr);
  for(i=0;i<qlen;i++) fputc("ACGTN"[query[i]], stderr);
  fputc('\t', stderr);
  fprintf(stderr, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
          qsc, qec,
          mm, mi, o, e, dir,
          -1, -1, -1, -1);
          */

  return v->process(b, a, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best);
}

void
tmap_vsw_wrapper_destroy(tmap_vsw_wrapper_t *v)
{
  //v->~AffineSWOptimization();
  delete v;
}

