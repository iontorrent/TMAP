#ifndef FMAP_MAP1_H_
#define FMAP_MAP1_H_

#include <config.h>

/*! 
  BWA-like (short-read) Mapping Algorithm
  */

/*! 
  @details  determines how to output multiple alignments
  */
enum {
    FMAP_MAP1_ALN_OUTPUT_MODE_BEST = 0, /*!< Output an alignment only if it is uniquely the best */
    FMAP_MAP1_ALN_OUTPUT_MODE_BEST_RAND = 1, /*!< Output a random best scoring alignment */
    FMAP_MAP1_ALN_OUTPUT_MODE_BEST_ALL = 2, /*!< Output all the alignments with the best score */
    FMAP_MAP1_ALN_OUTPUT_MODE_ALL = 3, /*!< Output all alignments */
};

/*! 
  structure to store the command line options for 'fmap map1'
  */
typedef struct {
    char **argv;  /*!< the command line argv structure */
    int argc;  /*!< the number of command line arguments passed */
    char *fn_fasta;  /*!< the fasta reference file name (-f) */
    char *fn_reads;  /*!< the reads file name (-r) */
    int32_t reads_format;  /*!< the reads file format (-F)  */
    int32_t seed_length;  /*!< the k-mer length to seed CALs (-l) */
    int32_t seed_max_mm;  /*!< maximum number of mismatches in hte seed (-k) */
    int32_t max_mm;  /*!< maximum number of mismatches (-m) */
    double max_mm_frac;  /*!< maximum (read length) fraction of mismatches (-m) */
    int32_t max_gapo;  /*!< maximum number of indel opens (-o) */
    double max_gapo_frac;  /*!< maximum (read length) fraction of indel opens (-o) */
    int32_t max_gape;  /*!< maximum number of indel extensions (-e) */
    double max_gape_frac;  /*!< maximum fraction of indel extensions (-e) */
    int32_t pen_mm;  /*!< the mismatch penalty (-M) */
    int32_t pen_gapo;  /*!< the indel open penalty (-O) */
    int32_t pen_gape;  /*!< the indel extension penalty (-E) */
    int32_t max_cals_del;  /*!< the maximum number of CALs to extend a deletion (-d) */
    int32_t indel_ends_bound;  /*!< indels are not allowed within INT number of bps from the end of the read (-i) */
    int32_t max_best_cals;  /*!< stop searching when INT optimal CALs have been found (-b) */
    int32_t reads_queue_size;  /*!< the reads queue size (-q) */
    int32_t max_entries;  /*!< maximum number of alignment nodes (-Q) */
    int32_t num_threads;  /*!< the number of threads (-n) */
    int32_t aln_output_mode;  /*!< specifies how to choose alignments (-a)  */
    int32_t input_compr;  /*!< the input compression type (-j and -z) */
    int32_t output_compr;  /*!< the output compression type (-J and -Z) */
    key_t shm_key;  /*!< the shared memory key (-s) */
} fmap_map1_opt_t;
/*! 
  @details In the CIGAR array, each element is a 32-bit integer. The
  lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
  length of a CIGAR.
  */
typedef struct {
    uint32_t score;  /*!< the current alignment score */
    uint16_t n_mm;  /*!< the current number of mismatches  */
    uint16_t n_gapo;  /*!< the current number of gap opens */
    uint16_t n_gape;  /*!< the current number of gap extensions */
    uint8_t mapq;  /*!< the mapping quality */
    uint8_t strand;  /*!< the strand of the alignment */
    uint32_t k;  /*!< the lower range of the SA interval */
    uint32_t l;  /*!< the upper range of the SA interval */
    uint32_t n_cigar;  /*!< the length of the cigar array */
    uint32_t *cigar;  /*!< the cigar array */
} fmap_map1_aln_t;
#ifdef HAVE_LIBPTHREAD
/*! 
  data to be passed to a thread
  */
typedef struct {
    fmap_seq_t **seq_buffer;  /*!< the buffer of sequences */
    int32_t seq_buffer_length;  /*!< the buffer length */
    fmap_map1_aln_t ***alns;  /*!< alignments for each sequence */
    fmap_bwt_t *bwt[2];  /*!< pointer to the BWT indices (forward/reverse) */
    int32_t tid;  /*!< the zero-based thread id */
    fmap_map1_opt_t *opt;  /*!< the options to this program */
} fmap_map1_thread_data_t;
#endif

/*! 
  main-like function for 'fmap map1'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
fmap_map1_main(int argc, char *argv[]);
#endif
