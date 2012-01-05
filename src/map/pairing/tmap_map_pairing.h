/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_PAIRING_H
#define TMAP_MAP_PAIRING_H

/*!
  the pairing strand orientation (strandedness)
  */
enum {
    TMAP_MAP_PAIRING_SAME_STRAND=0, /*!< the pairs should be on the same strand */
    TMAP_MAP_PAIRING_OPPOSITE_STRAND=1 /*!< the pairs should be on the opposite strand */
};

/*!
  the pairing position orientation (positioning)
  */
enum {
    TMAP_MAP_PAIRING_POSITIONING_AB=0, /*!< the first end (A) should be before the second end (B) */
    TMAP_MAP_PAIRING_POSITIONING_BA, /*!< the second end (B) should be before the first end (A) */
    TMAP_MAP_PAIRING_POSITIONING_NONE /*!< no positioning is required */
};

/*!
  given sets of seeds for two ends of pair, scores all possible pairs of seeds and
  filters accordingly.
  @param  one      the seeds for the first end (A)
  @param  two      the seeds for the second end (B)
  @param  one_seq  the sequence for the first end (A)
  @param  two_seq  the sequence for the second end (B)
  @param  rand     the random number generator
  @param  opt      the program parameters
 */
void
tmap_map_pairing_pick_pairs(tmap_map_sams_t *one, tmap_map_sams_t *two, tmap_seq_t *one_seq, tmap_seq_t *two_seq,
                            tmap_rand_t *rand, tmap_map_opt_t *opt);

#endif
