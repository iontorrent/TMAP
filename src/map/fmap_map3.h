#ifndef FMAP_MAP3_H_
#define FMAP_MAP3_H_

#include <config.h>
#include <sys/types.h>

/*! 
  SSAHA2-like Mapping Algorithm
  */


/*!
  Structure for a final hit
  */ 
typedef struct {
    uint16_t strand:1; /*!< the strand */
    uint32_t seqid;  /*!< the sequence index (0-based) */
    uint32_t pos; /*!< the position (0-based) */
    int32_t score; /*!< the alignment score */
    int32_t score_subo; /*!< the alignment score of the sub-optimal hit */
    uint8_t mapq; /*!< the mapping quality */
    uint16_t n_seeds:15; /*!< the number seeds in this hit */
    int32_t n_cigar; /*!< the number of cigar operators */
    uint32_t *cigar; /*!< the cigar operator array */
} fmap_map3_hit_t;

#ifdef HAVE_LIBPTHREAD
/*! 
  data to be passed to a thread
  */
typedef struct {
    fmap_seq_t **seq_buffer;  /*!< the buffer of sequences */
    fmap_map_sams_t **sams;  /*!< the alignments to output */
    int32_t seq_buffer_length;  /*!< the buffer length */
    fmap_refseq_t *refseq; /*< pointer to the packed referene sequence (forward) */
    fmap_bwt_t *bwt;  /*!< pointer to the BWT indices (reverse) */
    fmap_sa_t *sa;  /*< pointer to the SA (reverse) */
    int32_t tid;  /*!< the zero-based thread id */
    fmap_map_opt_t *opt;  /*!< the options to this program */
} fmap_map3_thread_data_t;
#endif

/*!
  Returns the inferred seed length given the reference length
  @param  ref_len  the reference length
  @return          the estimated seed length
  */
int32_t
fmap_map3_get_seed_length(uint64_t ref_len);

/*! 
  main-like function for 'fmap map3'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
fmap_map3_main(int argc, char *argv[]);
#endif
