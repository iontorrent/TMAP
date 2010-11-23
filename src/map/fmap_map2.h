#ifndef FMAP_MAP2_H_
#define FMAP_MAP2_H_

#include <stdlib.h>
#include <sys/ipc.h>

/*! 
  The BWA-like (long-read) Mapping Algorithm
  */

/*! 
  structure to store the command line options for 'fmap map2'
  */
typedef struct {
    char **argv;  /*!< the command line argv structure */
    int argc;  /*!< the number of command line arguments passed */
    char *fn_fasta;  /*!< the fasta reference file name (-f) */
    char *fn_reads;  /*!< the reads file name (-r) */
    int32_t reads_format;  /*!< the reads file format (-F) */
    int32_t score_match;  /*!< the match score (-A) */
    int32_t pen_mm;  /*!< the mismatch penalty (-M) */
    int32_t pen_gapo;  /*!< the indel open penalty (-O) */
    int32_t pen_gape;  /*!< the indel extension penalty (-E) */
    int32_t fscore;  /*!< the flow score penalty (-X) */
    double yita;  /*!< the error recurrence coefficient (-y)  */
    //double mask_level;  /*!< the mask level (-m) */
    double length_coef;  /*!< the coefficient of length-threshold adjustment (-c) */
    int32_t bw;  /*!< the band width (-w)  */
    int32_t score_thr;  /*!< the score threshold (match-score-scaled) (-T) */
    int32_t max_seed_intv;  /*!< the maximum seed interval (-S) */
    int32_t z_best;  /*!< the number of top scoring hits to keep (-b) */
    int32_t seeds_rev;  /*!< the maximum number of seeds for which reverse alignment is triggered (-N) */
    int32_t aln_global; /*!< align the full read (-g) */
    int32_t reads_queue_size;  /*!< the reads queue size (-q) */
    int32_t num_threads;  /*!< the number of threads (-n) */
    int32_t aln_output_mode;  /*!< specifies how to choose alignments (-a)  */
    int32_t input_compr;  /*!< the input compression type (-j and -z) */
    int32_t output_compr;  /*!< the output compression type (-J and -Z) */
    key_t shm_key;  /*!< the shared memory key (-s) */
} fmap_map2_opt_t;

/*! 
  data to be passed to a thread
  */
typedef struct {
    fmap_seq_t **seq_buffer;  /*!< the buffer of sequences */
    int32_t seq_buffer_length;  /*!< the buffer length */
    struct __fmap_map2_sam_t **sams;  /*!< the sam alignments for each sequence */
    fmap_refseq_t *refseq;  /*!< pointer to the reference sequence (forward) */
    fmap_bwt_t *bwt[2];  /*!< pointer to the BWT indices (forward/reverse) */
    fmap_sa_t *sa[2];  /*!< pointer to the SA (forward/reverse) */
    int32_t tid;  /*!< the zero-based thread id */
    fmap_map2_opt_t *opt;  /*!< the options to this program */
} fmap_map2_thread_data_t;

/*!
  Prints the usage of map2
  @param  opt  the current options
  @return      always 1
  */
int
fmap_map2_usage(fmap_map2_opt_t *opt);

/*!
  Gets the initialized options
  @return  pointer to the initialized options
  */
fmap_map2_opt_t *
fmap_map2_opt_init();

/*!
  Destroys the memory associated with these options
  @param  opt  pointer to the options
  */
void
fmap_map2_opt_destroy(fmap_map2_opt_t *opt);

/*!
  Parses the command line options and stores them in the options structure
  @param  argc  the number of arguments
  @param  argv  the argument list
  @param  opt   pointer to the options
  @return       1 if successful, 0 otherwise
  */
int32_t
fmap_map2_opt_parse(int argc, char *argv[], fmap_map2_opt_t *opt);

/*!
  Checks that all options are within range
  @param  opt   pointer to the options
  */
void
fmap_map2_opt_check(fmap_map2_opt_t *opt);

/*! 
  main-like function for 'fmap map2'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
fmap_map2_main(int argc, char *argv[]);

#endif
