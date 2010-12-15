#ifndef FMAP_MAP_ALL_H_
#define FMAP_MAP_ALL_H_

#include <sys/types.h>

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
    fmap_map_opt_t *opt;  /*!< the options to this program */    
} fmap_map_all_thread_data_t;

/*!
  Parses the command line options and stores them in the options structure
  @param  argc  the number of arguments
  @param  argv  the argument list
  @param  opt   pointer to the options
  @return       1 if successful, 0 otherwise
  */
int32_t
fmap_map_all_opt_parse(int argc, char *argv[], fmap_map_opt_t *opt);

/*! 
  main-like function for 'fmap map_all'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
fmap_map_all_main(int argc, char *argv[]);

#endif 
