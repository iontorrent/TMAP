#include <stdlib.h>
#include "../util/tmap_alloc.h"
#include "tmap_vsw_definitions.h"

tmap_vsw_opt_t*
tmap_vsw_opt_init(int32_t score_match, int32_t pen_mm, int32_t pen_gapo, int32_t pen_gape, int32_t score_thr)
{
  tmap_vsw_opt_t *opt;

  opt = tmap_malloc(sizeof(tmap_vsw_opt_t), "opt");
  opt->score_match = score_match;
  opt->pen_mm = pen_mm;
  opt->pen_gapo = pen_gapo;
  opt->pen_gape = pen_gape;
  opt->score_thres = score_thr;

  return opt;
}

void
tmap_vsw_opt_destroy(tmap_vsw_opt_t *opt)
{
  free(opt);
}
