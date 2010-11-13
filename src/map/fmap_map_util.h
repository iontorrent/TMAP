#ifndef FMAP_MAP_UTIL_H_
#define FMAP_MAP_UTIL_H_

/*!
  The various modes to modify the alignment score
  */
enum {
    FMAP_MAP_UTIL_ALN_MODE_UNIQ_BEST      = 0,  /*!< output an alignment only if its score is better than all other alignments */
    FMAP_MAP_UTIL_ALN_MODE_RAND_BEST      = 1,  /*!< output a random best scoring alignment */
    FMAP_MAP_UTIL_ALN_MODE_ALL_BEST       = 2,  /*!< output all the alignments with the best score */
    FMAP_MAP_UTIL_ALN_MODE_ALL            = 3   /*!< Output all alignments */
};

// TODO: have a general struct that stores an alignment

#endif
