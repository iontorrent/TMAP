#ifndef FMAP_FSW_H_
#define FMAP_FSW_H_

// We have 6-bits total, so 3-bits for above, and 3-bits for below
#define FMAP_FSW_MAX_OFFSET 7

#define FMAP_FSW_SET_SCORE_INF(s) (s).match_score = (s).ins_score = (s).del_score = FMAP_SW_MINOR_INF
#define FMAP_FSW_SET_FROM(s, from) (s).match_from = (s).ins_from = (s).del_from = from 
#define FMAP_FSW_SET_OFFSET(s, offset) (s).match_offset = (s).ins_offset = (s).del_offset = offset
#define FMAP_FSW_INIT_CELL(s, offset, from) (FMAP_FSW_SET_FROM(s, from), FMAP_FSW_SET_OFFSET(s, offset))

/*!
 *   From which cell; helps recovert the best scoring path.
 *     */
enum {
    FMAP_FSW_FROM_M = 0, /*!< from a mismatch cell*/
    FMAP_FSW_FROM_I = 1, /*!< from an insertion cell */
    FMAP_FSW_FROM_D = 2, /*!< from a deletion cell */
    FMAP_FSW_FROM_S = 3  /*!< from a start cell */
};

/*!
  @details  the valid offset in this implementation is only +-7
  */
typedef struct {
    uint8_t match_from:2; /*!< from cell for match */
    uint8_t match_offset:6; /*!< the difference between the original base call and the aligned base call, shifted by the offset*/
    uint8_t ins_from:2; /*!< from cell for insertion */
    uint8_t ins_offset:6; /*!< the difference between the original base call and the aligned base call, shifted by the offset*/
    uint8_t del_from:2; /*!< from cell for deletion */
    uint8_t del_offset:6; /*!< the difference between the original base call and the aligned base call, shifted by the offset*/
} fmap_fsw_dpcell_t;

/*!
  Stores the score for the current cell
  */
typedef struct
{
  int64_t match_score; /*!< match score */
  int64_t ins_score; /*!< insertion score */
  int64_t del_score; /*! <deletion score */
} fmap_fsw_dpscore_t;

/*!
  Parameters for the Smith-Waterman alignment.
  @details  The gap penalties should be positive (they will be subtracted)
  while the substitution score should be positive for matches, and negative for
  mismatches (they will be added).
  */
typedef struct
{
  int32_t gap_open; /*!< gap open penalty (positive) */
  int32_t gap_ext; /*!< gap extension penalty (positive) */
  int32_t gap_end; /*!< gap end penalty (positive */
  int32_t *matrix; /*!< substitution matrix (see details)*/
  int32_t fscore; /*!< the flow score (positive) */
  uint8_t offset; 
  int32_t row; /*!< the alphabet size */
  int32_t band_width; /*!< for Smith-Waterman banding */
} fmap_fsw_param_t;

/*!
  The best-scoring alignment path.
  */
typedef struct
{
  int32_t i; /*!< the flowgram1 index (0-based) */
  int32_t j; /*!< the seq2 index (0-based) */
  uint8_t ctype; /*!< the edit operator applied */
} fmap_fsw_path_t;

/*
   Fills in a row for one flow-space Smith-Waterman alignment
   @param  flow_base     the base associated with this flow
   @param  base_call    the number of bases called
   @param  flow_signal   the flow signal of this flow (100*signal)
   @param  seq2         the 2-bit DNA reference sequence 
   @param  len2         the length of the second sequence
   @param  ap           the alignment parameters
   @param  sub_dpcell   pre-allocated DP cells of minimum dimensions [base_call+2*(ap->offset+1),len2+1]
   @param  sub_score    pre-allocated DP scores of minimum dimensions [base_call+2*(ap->offset+1),len2+1]
   @param  dpcell_last  the last cell row in the DP matrix
   @param  score_last   the last score row in the DP matrix
   @param  dpcell_curr  the current cell row in the DP matrix
   @param  score_curr   the current score row in the DP matrix
   @details             this assumes that the ap parameter scores have been multiplied by 100
   */
inline void
fmap_fsw_sub_core(uint8_t flow_base, uint8_t base_call, uint16_t flow_signal,
                  uint8_t *seq2, int32_t len2,
                  const fmap_fsw_param_t *ap,
                  fmap_fsw_dpcell_t **sub_dpcell,
                  fmap_fsw_dpscore_t **sub_score,
                  fmap_fsw_dpcell_t *dpcell_last,
                  fmap_fsw_dpscore_t *score_last,
                  fmap_fsw_dpcell_t *dpcell_curr,
                  fmap_fsw_dpscore_t *score_curr);

/*
   Performs global flow-space Smith-Waterman alignment
   @param  flow        for each flow, the 2-bit DNA base flowed
   @param  bc1        for each flow, the number of bases called [0-255]
   @param  flowgram1   for each flow, the flowgram signal (100*signal) 
   @param  flowlen1    the number of flows
   @param  seq2       the 2-bit DNA reference sequence 
   @param  len2       the length of the second sequence
   @param  ap         the alignment parameters
   @param  path       the returned alignment path with maximum length of (1 + [len2 * (flowlen1 + 1) * (ap->offset + 1)])
   @param  path_len   the returned path_len
   @return            the returned alignment score
   @details           this assumes that the ap parameter scores have been multiplied by 100
   */
int64_t
fmap_fsw_global_core(uint8_t *flow, uint8_t *bc1, uint16_t *flowgram1, int32_t flowlen1,
                     uint8_t *seq2, int32_t len2,
                     const fmap_fsw_param_t *ap,
                     fmap_fsw_path_t *path, int32_t *path_len);

/*
   Performs flow-space Smith-Waterman alignment extension
   @param  flow          for each flow, the 2-bit DNA base flowed
   @param  bc1          for each flow, the number of bases called [0-255]
   @param  flowgram1     for each flow, the flowgram signal (100*signal) 
   @param  flowlen1      the number of flows
   @param  seq2         the 2-bit DNA reference sequence 
   @param  len2         the length of the second sequence
   @param  min_flowlen1  the minimum number of flows to align from the start
   @param  min_len2     the minimum number of bases from seq2 to align from the start
   @param  prev_score   the alignment score up to this point
   @param  ap           the alignment parameters
   @param  path         the returned alignment path with maximum length of (1 + [len2 * (flowlen1 + 1) * (ap->offset + 1)])
   @param  path_len     the returned path_len
   @return              the returned alignment score
   @details             this assumes that the ap parameter scores have been multiplied by 100.  
   */
int64_t
fmap_fsw_extend_core(uint8_t *flow, uint8_t *bc1, uint16_t *flowgram1, int32_t flowlen1,
                     uint8_t *seq2, int32_t len2,
                     int32_t min_flowlen1, int32_t min_len2, int64_t prev_score,
                     const fmap_fsw_param_t *ap,
                     fmap_fsw_path_t *path, int32_t *path_len);

/*!
  Creates a cigar array from an alignment path
  @param  path      the Smith-Waterman alignment path
  @param  path_len  the Smith-Waterman alignment path length
  @param  n_cigar   pointer to the returned number of cigar operations
  @return           the cigar array, NULL if the path is NULL or the path length is zero 
  */
uint32_t *
emap_fsw_path2cigar(const fmap_fsw_path_t *path, int32_t path_len, int32_t *n_cigar);

/*! 
  main-like function for 'fmap fsw'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */

int
fmap_fsw_main(int argc, char *argv[]);

#endif
