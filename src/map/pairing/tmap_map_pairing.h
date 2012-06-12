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
  @param  one  the mapping for the first end
  @param  two  the mapping for the second end
  @param  strandedness the strandedness orientation
  @return  1 if the two mappings are consistent with the strandedness parameter, 0 otherwise 
 */
int32_t
tmap_map_pairing_get_strand_diff(tmap_map_sam_t *one, tmap_map_sam_t *two, int32_t strandedness);

/*!
  Gets the positional difference between two alignments
  @param  one  the mapping for the first end
  @param  two  the mapping for the second end
  @param  one_len  the sequence length for the first end
  @param  two_len  the sequence length for the second end
  @param  strandedness the strandedness orientation
  @param  positioning the positioning orientation
  @returns the position difference assuming the two reads are from the same contig and match strandedness
 */
int32_t
tmap_map_pairing_get_position_diff(tmap_map_sam_t *one, tmap_map_sam_t *two, int32_t one_len, int32_t two_len,
                                   int32_t strandedness, int32_t positioning);

/*!
  performs read rescue
  @param  refseq   the reference sequence
  @param  one      the seeds for the first end (A)
  @param  two      the seeds for the second end (B)
  @param  one_seq  the sequence for the first end (foward/reverse compliment) (A)
  @param  two_seq  the sequence for the second end (forward/reverse compliment) (B)
  @param  rand     the random number generator
  @param  opt      the program parameters
  @return          0 if no reads are rescued, 1 if only end 1 was rescued, 2 if only end 2 was 
  rescued, and 3 if both ends were rescued
  */
int32_t
tmap_map_pairing_read_rescue(tmap_refseq_t *refseq, 
                             tmap_map_sams_t *one, tmap_map_sams_t *two, 
                             tmap_seq_t *one_seq[2], tmap_seq_t *two_seq[2], 
                             tmap_rand_t *rand, tmap_map_opt_t *opt);

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
