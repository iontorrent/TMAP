#ifndef TMAP_MAP2_H_
#define TMAP_MAP2_H_

#include <stdlib.h>
#include <sys/types.h>

/*! 
  The BWA-like (long-read) Mapping Algorithm
  */

#ifdef HAVE_LIBPTHREAD
/*! 
  data to be passed to a thread
  */
typedef struct {
    tmap_seq_t **seq_buffer;  /*!< the buffer of sequences */
    int32_t seq_buffer_length;  /*!< the buffer length */
    tmap_map_sams_t **sams;  /*!< the sam alignments for each sequence */
    tmap_refseq_t *refseq;  /*!< pointer to the reference sequence (forward) */
    tmap_bwt_t *bwt[2];  /*!< pointer to the BWT indices (forward/reverse) */
    tmap_sa_t *sa[2];  /*!< pointer to the SA (forward/reverse) */
    int32_t tid;  /*!< the zero-based thread id */
    tmap_map_opt_t *opt;  /*!< the options to this program */
} tmap_map2_thread_data_t;
#endif

/*! 
  main-like function for 'tmap map2'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
tmap_map2_main(int argc, char *argv[]);

#endif
