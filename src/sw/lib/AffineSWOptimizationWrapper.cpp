#define __STDC_LIMIT_MACROS // Seriously, I want these, now
#include <stdint.h>
//#include <stdlib.h>
//#include <stdio.h>
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
  /*
  fprintf(stderr, "query=%d qlen=%d target=%d tlen=%d\n",
          (NULL == query) ? 0 : 1,
          qlen,
          (NULL == target) ? 0 : 1,
          tlen);
  if(NULL == query || 0 == qlen || NULL == target || 0 == tlen) {
      exit(1);
  }
  */

  // Top coder style
  /*
  int i;
  for(i=0;i<tlen;i++) fputc("ACGTN"[target[i]], stderr);
  fputc('\t', stderr);
  for(i=0;i<qlen;i++) fputc("ACGTN"[query[i]], stderr);
  fputc('\t', stderr);
  fprintf(stderr, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
          qsc, qec,
          mm, mi, o, e, dir,
          -1, -1, -1, -1);
  */

  return v->process(target, tlen, query, qlen, qsc, qec, mm, mi, o, e, dir, opt, te, qe, n_best);
}

void
tmap_vsw_wrapper_destroy(tmap_vsw_wrapper_t *v)
{
  //v->~AffineSWOptimization();
  delete v;
}

int
tmap_vsw_wrapper_get_max_qlen(tmap_vsw_wrapper_t *v)
{
  return v->getMaxQlen();
}

int
tmap_vsw_wrapper_get_max_tlen(tmap_vsw_wrapper_t *v)
{
  return v->getMaxTlen();
}
