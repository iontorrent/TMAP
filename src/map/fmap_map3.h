#ifndef FMAP_MAP3_H_
#define FMAP_MAP3_H_

#include <config.h>

/*! 
  SSAHA2-like Mapping Algorithm
  */

/*! 
  structure to store the command line options for 'fmap map3'
  */
typedef struct {
    char **argv;  /*!< the command line argv structure */
    int argc;  /*!< the number of command line arguments passed */
    char *fn_fasta;  /*!< the fasta reference file name (-f) */
    char *fn_reads;  /*!< the reads file name (-r) */
    int32_t reads_format;  /*!< the reads file format (-F)  */

    int32_t seed_length; /*!< the kmer seed length (-l) */
    int32_t max_seed_hits; /*!< the maximum number of hits returned by a seed (-S) */
    int32_t max_seed_band; /*!< the band to group seeds (-b)*/
    int32_t sw_offset; /*!< the extra bases to add before and after the target during Smith-Waterman (-w) */
    // TODO: vary the offsets of the kmers

    int32_t score_match;  /*!< the match score (-A) */
    int32_t pen_mm;  /*!< the mismatch penalty (-M) */
    int32_t pen_gapo;  /*!< the indel open penalty (-O) */
    int32_t pen_gape;  /*!< the indel extension penalty (-E) */
    int32_t score_thr;  /*!< the score threshold (match-score-scaled) (-T) */
    int32_t reads_queue_size;  /*!< the reads queue size (-q) */
    int32_t num_threads;  /*!< the number of threads (-n) */
    int32_t aln_output_mode;  /*!< specifies how to choose alignments (-a)  */
    int32_t input_compr;  /*!< the input compression type (-j and -z) */
    int32_t output_compr;  /*!< the output compression type (-J and -Z) */
    key_t shm_key;  /*!< the shared memory key (-s) */
} fmap_map3_opt_t;

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

/*!
 Stucture for holding alignment hits
 */
typedef struct {
    int32_t n; /*!< the number of hits */
    fmap_map3_hit_t *hits; /*!< array of hits */
} fmap_map3_aln_t;

#ifdef HAVE_LIBPTHREAD
/*! 
  data to be passed to a thread
  */
typedef struct {
    // TODO: return alignments
    fmap_seq_t **seq_buffer;  /*!< the buffer of sequences */
    fmap_map3_aln_t **alns;  /*!< the alignments to output */
    int32_t seq_buffer_length;  /*!< the buffer length */
    fmap_refseq_t *refseq; /*< pointer to the packed referene sequence (forward) */
    fmap_bwt_t *bwt;  /*!< pointer to the BWT indices (reverse) */
    fmap_sa_t *sa;  /*< pointer to the SA (reverse) */
    int32_t tid;  /*!< the zero-based thread id */
    fmap_map3_opt_t *opt;  /*!< the options to this program */
} fmap_map3_thread_data_t;
#endif

/*! 
  main-like function for 'fmap map3'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
fmap_map3_main(int argc, char *argv[]);
#endif
