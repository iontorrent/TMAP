#ifndef FMAP_MAP_ALL_H_
#define FMAP_MAP_ALL_H_

#include "fmap_map1.h"
#include "fmap_map2.h"
#include "fmap_map3.h"

enum {
    FMAP_MAP_ALL_ALGO_NONE = 0x0,  /*!< dummy algorithm */
    FMAP_MAP_ALL_ALGO_MAP1 = 0x1,  /*!< the map1 algorithm */
    FMAP_MAP_ALL_ALGO_MAP2 = 0x2,  /*!< the map2 algorithm */
    FMAP_MAP_ALL_ALGO_MAP3 = 0x4,  /*!< the map3 algorithm */
};

typedef struct {
    int32_t algo_id; /*< the algorithm id used to obtain this hit */

    uint16_t strand:1; /*!< the strand */
    uint32_t seqid;  /*!< the sequence index (0-based) */
    uint32_t pos; /*!< the position (0-based) */
    int32_t score; /*!< the alignment score */
    uint8_t mapq; /*!< the mapping quality */
    int32_t n_cigar; /*!< the number of cigar operators */
    uint32_t *cigar; /*!< the cigar operator array */

    // map1
    uint16_t n_mm;  /*!< the current number of mismatches  */
    uint16_t n_gapo;  /*!< the current number of gap opens */
    uint16_t n_gape;  /*!< the current number of gap extensions */
    // map2
    uint8_t XF;  /*!< support for the forward/reverse alignment (1-forward 2-reverse 3-both) */
    int32_t XI;  /*!< the suffix interval size */
    // map3
    // - none
    // map2 & map3
    uint16_t n_seeds; /*!< the number seeds in this hit */
    int32_t score_subo; /*!< the alignment score of the sub-optimal hit */
} fmap_map_all_hit_t;

typedef struct {
    fmap_map_all_hit_t *hits;
    int32_t n;
} fmap_map_all_aln_t;

/*!
  Structure to store the command line options for 'fmap mapall'
  */
typedef struct {
    // global
    char **argv;  /*!< the command line argv structure */
    int argc;  /*!< the number of command line arguments passed */
    uint32_t algos;  /*!< the algorithms that should be run, bit-packed */

    // common options
    char *fn_fasta;  /*!< the fasta reference file name (-f) */
    char *fn_reads;  /*!< the reads file name (-r) */
    int32_t reads_format;  /*!< the reads file format (-F)  */
    int32_t score_match;  /*!< the match score (-A) */
    int32_t pen_mm;  /*!< the mismatch penalty (-M) */
    int32_t pen_gapo;  /*!< the indel open penalty (-O) */
    int32_t pen_gape;  /*!< the indel extension penalty (-E) */
    int32_t sw_offset;  /*!< the band width (-w)  */
    int32_t aln_global; /*!< align the full read (-g) */
    int32_t reads_queue_size;  /*!< the reads queue size (-q) */
    int32_t num_threads;  /*!< the number of threads (-n) */
    int32_t aln_output_mode;  /*!< specifies how to choose alignments (-a)  */
    int32_t input_compr;  /*!< the input compression type (-j and -z) */
    int32_t output_compr;  /*!< the output compression type (-J and -Z) */
    key_t shm_key;  /*!< the shared memory key (-s) */

    // map all specific options
    int32_t dup_window; /*!< remove duplicate alignments from different algorithms within this bp window (-W) */
    int32_t aln_output_mode_ind; /*!< apply the output filter for each algorithm separately */

    // mapping algorithm specific options
    fmap_map1_opt_t *opt_map1; /*!< map 1 options */
    fmap_map2_opt_t *opt_map2; /*!< map 2 options */
    fmap_map3_opt_t *opt_map3; /*!< map 3 options */
} fmap_map_all_opt_t;

/*! 
 data to be passed to a thread                         */
typedef struct {                                            
    fmap_seq_t **seq_buffer;  /*!< the buffer of sequences */    
    int32_t seq_buffer_length;  /*!< the buffer length */
    fmap_map_all_aln_t **alns;  /*!< the alignments for each sequence */
    fmap_refseq_t *refseq;  /*!< pointer to the reference sequence (forward) */
    fmap_bwt_t *bwt[2];  /*!< pointer to the BWT indices (forward/reverse) */
    fmap_sa_t *sa[2];  /*!< pointer to the SA (forward/reverse) */    
    int32_t tid;  /*!< the zero-based thread id */
    fmap_map_all_opt_t *opt;  /*!< the options to this program */    
} fmap_map_all_thread_data_t;



    /*!
      Prints the usage of map_all
      @param  opt  the current options
      @return      always 1
      */
    int
    fmap_map_all_usage(fmap_map_all_opt_t *opt);

    /*!
      Gets the initialized options
      @return  pointer to the initialized options
      */
    fmap_map_all_opt_t *
    fmap_map_all_opt_init();

    /*!
      Destroys the memory associated with these options
      @param  opt  pointer to the options
      */
    void
    fmap_map_all_opt_destroy(fmap_map_all_opt_t *opt);

    /*!
      Parses the command line options and stores them in the options structure
      @param  argc  the number of arguments
      @param  argv  the argument list
      @param  opt   pointer to the options
      @return       1 if successful, 0 otherwise
      */
    int32_t
    fmap_map_all_opt_parse(int argc, char *argv[], fmap_map_all_opt_t *opt);

    /*!
      Checks that all options are within range
      @param  opt   pointer to the options
      */
    void
    fmap_map_all_opt_check(fmap_map_all_opt_t *opt);

    /*! 
      main-like function for 'fmap map_all'
      @param  argc  the number of arguments
      @param  argv  the argument list
      @return       0 if executed successful
      */
    int 
    fmap_map_all_main(int argc, char *argv[]);


#endif 
