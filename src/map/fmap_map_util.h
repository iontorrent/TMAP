#ifndef FMAP_MAP_UTIL_H_
#define FMAP_MAP_UTIL_H_

/*!
  The various modes to modify the alignment score
  */
enum {
    FMAP_MAP_UTIL_ALN_MODE_RAND           = 0,  /*!< Output a random alignment > */
    FMAP_MAP_UTIL_ALN_MODE_SCORE_LEN_NORM = 1,  /*!< Output the best scoring alignment(s) normalized by alignment length */
    FMAP_MAP_UTIL_ALN_MODE_SCORE          = 2,  /*!< Output the best scoring alignment(s) */
    FMAP_MAP_UTIL_ALN_MODE_LEN            = 3,  /*!< Output the longest alignment(s) */
    FMAP_MAP_UTIL_ALN_MODE_ALL            = 4   /*!< Output all alignments */
};

/*!
  @param  seq    the sequence associated with this score
  @param  score  the unmodified alignment score
  @param  mode   the alignment scoring mode
  */
inline double 
fmap_map_util_get_score(fmap_seq_t *seq, int32_t score, int32_t mode);

#endif
