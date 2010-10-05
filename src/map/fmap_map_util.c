#include <stdlib.h>
#include "../seq/fmap_seq.h"
#include "fmap_map_util.h"

inline double
fmap_map_util_get_score(fmap_seq_t *seq, int32_t score, int32_t mode)
{
  switch(mode) {
    case FMAP_MAP_UTIL_ALN_MODE_SCORE_LEN_NORM:
      return score / (double)fmap_seq_get_bases(seq)->l; 
      break;
    case FMAP_MAP_UTIL_ALN_MODE_LEN:
      return (double)fmap_seq_get_bases(seq)->l;
      break;
    case FMAP_MAP_UTIL_ALN_MODE_SCORE:
    case FMAP_MAP_UTIL_ALN_MODE_RAND: 
    case FMAP_MAP_UTIL_ALN_MODE_ALL: 
    default:
      return (double)score; 
      break;
  }
}
