/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_FSW_H
#define TMAP_FSW_H

#include "../seq/tmap_sff.h"

// We have 6-bits total, so 3-bits for above, and 3-bits for below
#define TMAP_FSW_MAX_OFFSET 7

#define TMAP_FSW_SET_SCORE_INF(s) (s).match_score = (s).ins_score = (s).del_score = TMAP_SW_MINOR_INF
#define TMAP_FSW_SET_FROM(s, from) (s).match_from = (s).ins_from = (s).del_from = from 
#define TMAP_FSW_SET_BC(s, bc) (s).match_bc = (s).ins_bc = (s).del_bc = bc
#define TMAP_FSW_INIT_CELL(s) (TMAP_FSW_SET_FROM(s,TMAP_FSW_FROM_S), TMAP_FSW_SET_BC(s, 0))

#define TMAP_FSW_MAX_PATH_LENGTH(ref_len, flow_len, offset) ((1 + (ref_len * (flow_len + 1) * (offset + 1))))

#define __tmap_fsw_gen_ap1(par, score_match, pen_mm, pen_gapo, pen_gape, fscore) do { \
    int32_t i; \
    for(i=0;i<25;i++) { \
        (par).matrix[i] = -100 * pen_mm; \
    } \
    for(i=0;i<4;i++) { \
        (par).matrix[i*5+i] = 100 * score_match; \
    } \
    (par).gap_open = 100 * pen_gapo; \
    (par).gap_ext = 100 *pen_gape; \
    (par).gap_end = 100 * pen_gape; \
    (par).fscore = 100 * fscore; \
    (par).row = 5; \
} while(0)

#define __tmap_fsw_gen_ap(par, opt) (__tmap_fsw_gen_ap1(par, (opt)->score_match, (opt)->pen_mm, (opt)->pen_gapo, (opt)->pen_gape), (opt)->fscore)

/*!
 *   From which cell; helps recovert the best scoring path.
 *     */
enum {
    TMAP_FSW_FROM_M = 0, /*!< from a match/mismatch cell*/
    TMAP_FSW_FROM_I = 1, /*!< from an insertion cell */
    TMAP_FSW_FROM_D = 2, /*!< from a deletion cell */
    TMAP_FSW_FROM_S = 3, /*!< from a start cell */
    TMAP_FSW_FROM_HP_PLUS = 4, /*!< from a match/mismatch cell, but a hp overcall error */
    TMAP_FSW_FROM_HP_MINUS = 5 /*!< from a match/mismatch cell, but a hp undercall error */
};

enum {
    TMAP_FSW_NO_JUSTIFY = 0, /*!< do not perform indel justification */
    TMAP_FSW_JUSTIFY_LEFT_REF = 1, /*!< justify 5' on the reference strand */
    TMAP_FSW_JUSTIFY_LEFT_READ = 2  /*!< justify 5' on the read strand */
};

/*!
  The path for the current cell
  */
typedef struct {
    uint8_t match_bc; /*!< the base call for a match */
    uint8_t match_from; /*!< the from cell in the lower 2 bits, and the column offset in the upper 6 bits */
    uint8_t ins_bc; /*!< the base call for a insertion */
    uint8_t ins_from; /*!< the from cell in the lower 2 bits, and the column offset in the upper 6 bits */
    uint8_t del_bc; /*!< the base call for a deletion */
    uint8_t del_from; /*!< the from cell in the lower 2 bits, and the column offset in the upper 6 bits */
} tmap_fsw_dpcell_t;

/*!
  Stores the score for the current cell
  */
typedef struct {
    int64_t match_score; /*!< match score */
    int64_t ins_score; /*!< insertion score */
    int64_t del_score; /*! <deletion score */
} tmap_fsw_dpscore_t;

/*!
  Parameters for the Smith-Waterman alignment.
  @details  The gap penalties should be positive (they will be subtracted)
  while the substitution score should be positive for matches, and negative for
  mismatches (they will be added).
  */
typedef struct {
    int32_t gap_open; /*!< gap open penalty (positive) */
    int32_t gap_ext; /*!< gap extension penalty (positive) */
    int32_t gap_end; /*!< gap end penalty (positive */
    int32_t *matrix; /*!< substitution matrix (see details)*/
    int32_t fscore; /*!< the flow score (positive) */
    uint8_t offset; 
    int32_t row; /*!< the alphabet size */
    int32_t band_width; /*!< for Smith-Waterman banding */
} tmap_fsw_param_t;

/*!
  The best-scoring alignment path.
  */
typedef struct {
    // TODO: make 1-based to match tmap_sw.{c,h}
    int32_t i; /*!< the flowgram index (0-based) */
    int32_t j; /*!< the seq index (0-based) */
    uint8_t ctype; /*!< the edit operator applied */
} tmap_fsw_path_t;

typedef struct {
    uint8_t *flow_order; /*!< for each of the four flows, the 2-bit DNA base flowed */
    int32_t flow_order_len; /*!< the flow order length */
    uint8_t *base_calls; /*!< for each flow, the number of bases called [0-255] */
    uint16_t *flowgram; /*!< for each flow, the flowgram signal (100*signal)  */
    int32_t num_flows; /*!< the number of flows */
    int32_t key_index; /*!< the 0-based index of the key base if the key is part of the first or last flow (should be -1, 0, or num_flow-1). */
    int32_t key_bases; /*!< the number of bases part of the key_index flow that are explained by the key sequence */
} tmap_fsw_flowseq_t;

/*!
  Stores parameters for flow-space Smith-Waterman
  @param  flow         for each of the four flows, the 2-bit DNA base flowed 
  @param  base_calls  for each flow, the number of bases called [0-255] 
  @param  flowgram     for each flow, the flowgram signal (100*signal)  
  @param  num_flows    the number of flows 
  @param  key_index   the 0-based index of the key base if the key is part of the first or last flow (should be -1, 0, or num_flow-1). 
  @param  key_bases   the number of bases part of the key_index flow that are explained by the key sequence 
  @return             pointer to the initialized memory
  */
tmap_fsw_flowseq_t *
tmap_fsw_flowseq_init(uint8_t *flow_order, int32_t flow_order_len, uint8_t *base_calls, uint16_t *flowgram,
                      int32_t num_flows, int32_t key_index, int32_t key_bases);

/*!
  @param  fp      the file stream in which to print
  @param  flowseq  the flow sequence to print
  */
void
tmap_fsw_flowseq_print(tmap_file_t *fp, tmap_fsw_flowseq_t *flowseq);

/*!
  @param  flowseq  pointer to the flow sequence to destroy
  @details        this does a shallow-destroy
  */
void
tmap_fsw_flowseq_destroy_shallow(tmap_fsw_flowseq_t *flowseq);

/*!
  @param  flowseq  pointer to the flow sequence to destroy
  */
void
tmap_fsw_flowseq_destroy(tmap_fsw_flowseq_t *flowseq);

/*!
  Fills in a row for one flow-space Smith-Waterman alignment
  @param  seq                the 2-bit DNA reference sequence 
  @param  len                the length of the second sequence
  @param  flow_base           the base associated with this flow
  @param  base_call          the number of bases called (do not inclue the key base(s))
  @param  flow_signal         the flow signal of this flow (100*signal)
  @param  ap                 the alignment parameters
  @param  sub_dpcell         pre-allocated DP cells of minimum dimensions [base_call+2*(ap->offset+1),len+1]
  @param  sub_score          pre-allocated DP scores of minimum dimensions [base_call+2*(ap->offset+1),len+1]
  @param  dpcell_last        the last cell row in the DP matrix
  @param  score_last         the last score row in the DP matrix
  @param  dpcell_curr        the current cell row in the DP matrix
  @param  score_curr         the current score row in the DP matrix
  @param  path               the sub-alignment path, NULL if not required
  @param  path_len           the returned sub-alignment path length, 0 if path is NULL
  @param  best_ctype         the sub-cell from which to backtrace if path is not NULL
  @param  key_bases          the number of key bases that are part of this flow
  @param  flowseq_start_clip  1 to allow clipping at the start of the flow sequence, 0 otherwise
  @param  flowseq_end_clip    1 to allow clipping at the end of the flow sequence, 0 otherwise
  @details                   this assumes that the ap parameter scores have been multiplied by 100; only include non-key flows
  */
inline void
tmap_fsw_sub_core(uint8_t *seq, int32_t len,
                  uint8_t flow_base, uint8_t base_call, uint16_t flow_signal,
                  const tmap_fsw_param_t *ap,
                  tmap_fsw_dpcell_t **sub_dpcell,
                  tmap_fsw_dpscore_t **sub_score,
                  tmap_fsw_dpcell_t *dpcell_last,
                  tmap_fsw_dpscore_t *score_last,
                  tmap_fsw_dpcell_t *dpcell_curr,
                  tmap_fsw_dpscore_t *score_curr,
                  tmap_fsw_path_t *path, int32_t *path_len, int32_t best_ctype,
                  uint8_t key_bases,
                  int32_t flowseq_start_clip, int32_t flowseq_end_clip);

/*!
  Performs flow-space Smith-Waterman alignment extension
  @param  seq         the 2-bit DNA reference sequence 
  @param  len         the length of the second sequence
  @param  flowseq      the flow sequence structure to align
  @param  ap          the alignment parameters
  @param  path        the returned alignment path with maximum length of (1 + [len * (num_flows + 1) * (ap->offset + 1)])
  @param  path_len    the returned path_len
  @param  prev_score  the alignment score up to this point
  @return             the returned alignment score
  @details            this assumes that the ap parameter scores have been multiplied by 100; only include non-key flows
  */
int64_t
tmap_fsw_extend_core(uint8_t *seq, int32_t len,
                     tmap_fsw_flowseq_t *flowseq,
                     const tmap_fsw_param_t *ap,
                     tmap_fsw_path_t *path, int32_t *path_len, 
                     int32_t prev_score);

/*!
  Performs flow-space Smith-Waterman alignment extension and aligns the entire flowgram
  @param  seq         the 2-bit DNA reference sequence 
  @param  len         the length of the second sequence
  @param  flowseq      the flow sequence structure to align
  @param  ap          the alignment parameters
  @param  path        the returned alignment path with maximum length of (1 + [len * (num_flows + 1) * (ap->offset + 1)])
  @param  path_len    the returned path_len
  @param  prev_score  the alignment score up to this point
  @return             the returned alignment score
  @details            this assumes that the ap parameter scores have been multiplied by 100; only include non-key flows
  */
int64_t
tmap_fsw_extend_fitting_core(uint8_t *seq, int32_t len,
                             tmap_fsw_flowseq_t *flowseq,
                             const tmap_fsw_param_t *ap,
                             tmap_fsw_path_t *path, int32_t *path_len, 
                             int32_t prev_score);

/*!
  Performs flow-space Smith-Waterman alignment 
  @param  seq                the 2-bit DNA reference sequence 
  @param  len                the length of the second sequence
  @param  flowseq             the flow sequence structure to align
  @param  ap                 the alignment parameters
  @param  path               the returned alignment path with maximum length of (1 + [len * (num_flows + 1) * (ap->offset + 1)])
  @param  path_len           the returned path_len
  @param  flowseq_start_clip  1 to allow clipping at the start of the flow sequence, 0 otherwise
  @param  flowseq_end_clip    1 to allow clipping at the end of the flow sequence, 0 otherwise
  @return                    the returned alignment score
  @details                   this assumes that the ap parameter scores have been multiplied by 100; only include non-key flows
  */
int64_t
tmap_fsw_clipping_core(uint8_t *seq, int32_t len,
                       tmap_fsw_flowseq_t *flowseq,
                       const tmap_fsw_param_t *ap,
                       int32_t flowseq_start_clip, int32_t flowseq_end_clip,
                       tmap_fsw_path_t *path, int32_t *path_len);

/*!
  Creates a cigar array from an alignment path
  @param  path      the Smith-Waterman alignment path
  @param  path_len  the Smith-Waterman alignment path length
  @param  n_cigar   pointer to the returned number of cigar operations
  @param  rm_hp     1 if we are to remove hp edits (merge into indels), 0 otherwise
  @return           the cigar array, NULL if the path is NULL or the path length is zero 
  */
uint32_t *
tmap_fsw_path2cigar(const tmap_fsw_path_t *path, int32_t path_len, int32_t *n_cigar, int32_t rm_hp);

/*!
  Gets the pretty-print alignment
  @param  path           the alignment path
  @param  path_len       the alignment path length
  @param  flow_order      the flow order
  @param  flow_order_len  the flow order length
  @param  target         the 2-bit DNA reference sequence 
  @param  strand         0 for the forward strand, 1 for the reverse
  @param  ref            pointer to the returned reference string
  @param  read           pointer to the returned read string
  @param  aln            pointer to the returned alignment string
  @param  j_type         the indel justification method 
  */
void
tmap_fsw_get_aln(tmap_fsw_path_t *path, int32_t path_len,
                 uint8_t *flow_order, int32_t flow_order_len, uint8_t *target, uint8_t strand,
                 char **ref, char **read, char **aln, int32_t j_type);

/*!
  Pretty-prints an alignment
  @param  fp             the file pointer to print
  @param  score          the flow space alignment score
  @param  path           the alignment path
  @param  path_len       the alignment path length
  @param  flow_order      the flow order
  @param  flow_order_len  the flow order length
  @param  target         the 2-bit DNA reference sequence 
  @param  strand         0 for the forward strand, 1 for the reverse
  @param  j_type         the indel justification method 
  @param  sep            the field separator
  */
void 
tmap_fsw_print_aln(tmap_file_t *fp, int64_t score, tmap_fsw_path_t *path, int32_t path_len,
                   uint8_t *flow_order, int32_t flow_order_len, uint8_t *target, uint8_t strand, int32_t j_type, char sep);

/*!
  Create a structure for flow-space Smith Waterman from an sequence structure
  @param  fq             the sequence structure
  @param  flow_order      the flow order (in integer format); not used for an SFF
  @param  flow_order_len  the flow order length; not used for an SFF
  @return                pointer to the flowseq structure
  */
tmap_fsw_flowseq_t *
tmap_fsw_seq_to_flowseq(tmap_seq_t *seq, uint8_t *flow_order, int32_t flow_order_len);

/*!
  Frees the memory associated with this structure
  @param  flowseq  pointer the flow sequence to free
  */
void
tmap_fsw_flowseq_destroy(tmap_fsw_flowseq_t *flowseq);

#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
/*! 
  main-like function for 'tmap fsw'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */

int
tmap_fsw_main(int argc, char *argv[]);
#endif

#endif
