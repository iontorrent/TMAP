#ifndef FMAP_MAP1_H_
#define FMAP_MAP1_H_

#include <config.h>
#include <sys/types.h>

/*! 
  BWA-like (short-read) Mapping Algorithm
  */

#define FMAP_MAP1_FSW_BW 50

/*! 
  @details  determines how to output multiple alignments
  */
enum {
    FMAP_MAP1_ALN_OUTPUT_MODE_BEST = 0, /*!< Output an alignment only if it is uniquely the best */
    FMAP_MAP1_ALN_OUTPUT_MODE_BEST_RAND = 1, /*!< Output a random best scoring alignment */
    FMAP_MAP1_ALN_OUTPUT_MODE_BEST_ALL = 2, /*!< Output all the alignments with the best score */
    FMAP_MAP1_ALN_OUTPUT_MODE_ALL = 3, /*!< Output all alignments */
};

#ifdef HAVE_LIBPTHREAD
/*! 
  data to be passed to a thread
  */
typedef struct {
    fmap_seq_t **seq_buffer;  /*!< the buffer of sequences */
    int32_t seq_buffer_length;  /*!< the buffer length */
    fmap_map_sams_t **sams;  /*!< alignments for each sequence */
    fmap_refseq_t *refseq; /*!< pointer to the reference sequence */
    fmap_bwt_t *bwt[2];  /*!< pointer to the BWT indices (forward/reverse) */
    fmap_sa_t *sa; /*!< pointer to the SA (reverse) */
    int32_t tid;  /*!< the zero-based thread id */
    fmap_map_opt_t *opt;  /*!< the options to this program */
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
