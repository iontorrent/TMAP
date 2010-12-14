#ifndef FMAP_MAP_ALL_H_
#define FMAP_MAP_ALL_H_

#include <sys/types.h>
#include "fmap_map1.h"
#include "fmap_map2.h"
#include "fmap_map3.h"

/*!
  Structure to store the command line options for 'fmap mapall'
  */
typedef struct {
    // global
    char **argv;  /*!< the command line argv structure */
    int argc;  /*!< the number of command line arguments passed */
    uint32_t algos[2];  /*!< the algorithms that should be run in stage 1 and stage 2, bit-packed */

    // common options
    char *fn_fasta;  /*!< the fasta reference file name (-f) */
    char *fn_reads;  /*!< the reads file name (-r) */
    int32_t reads_format;  /*!< the reads file format (-F)  */
    int32_t score_match;  /*!< the match score (-A) */
    int32_t pen_mm;  /*!< the mismatch penalty (-M) */
    int32_t pen_gapo;  /*!< the indel open penalty (-O) */
    int32_t pen_gape;  /*!< the indel extension penalty (-E) */
    int32_t fscore;  /*!< the flow score penalty (-X) */
    int32_t bw;  /*!< the band width (-w)  */
    int32_t aln_global; /*!< align the full read (-g) */
    int32_t reads_queue_size;  /*!< the reads queue size (-q) */
    int32_t num_threads;  /*!< the number of threads (-n) */
    int32_t aln_output_mode;  /*!< specifies how to choose alignments (-a)  */
    char *sam_rg;  /*!< specifies the RG line in the SAM header (-R) */
    int32_t input_compr;  /*!< the input compression type (-j and -z) */
    int32_t output_compr;  /*!< the output compression type (-J and -Z) */
    key_t shm_key;  /*!< the shared memory key (-s) */

    // map all specific options
    int32_t dup_window; /*!< remove duplicate alignments from different algorithms within this bp window (-W) */
    int32_t aln_output_mode_ind; /*!< apply the output filter for each algorithm separately (-I) */

    // stage 1/2 mapping algorithm specific options
    fmap_map1_opt_t *opt_map1[2]; /*!< map 1 options */
    fmap_map2_opt_t *opt_map2[2]; /*!< map 2 options */
    fmap_map3_opt_t *opt_map3[2]; /*!< map 3 options */
} fmap_map_all_opt_t;

/*! 
  data to be passed to a thread                         */
typedef struct {                                            
    fmap_seq_t **seq_buffer;  /*!< the buffer of sequences */    
    int32_t seq_buffer_length;  /*!< the buffer length */
    fmap_map_sams_t **sams;  /*!< the alignments for each sequence */
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
