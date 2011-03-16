/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_ALL_H
#define TMAP_MAP_ALL_H

#include <sys/types.h>

#ifdef HAVE_LIBPTHREAD
#define TMAP_MAP_ALL_THREAD_BLOCK_SIZE 512
#endif

/*! 
  data to be passed to a thread                         */
typedef struct {                                            
    tmap_seq_t **seq_buffer;  /*!< the buffer of sequences */    
    int32_t seq_buffer_length;  /*!< the buffer length */
    tmap_map_sams_t **sams;  /*!< the alignments for each sequence */
    tmap_refseq_t *refseq;  /*!< pointer to the reference sequence (forward) */
    tmap_bwt_t *bwt[2];  /*!< pointer to the BWT indices (forward/reverse) */
    tmap_sa_t *sa[2];  /*!< pointer to the SA (forward/reverse) */    
    int32_t thread_block_size; /*!< the number of reads per thread to process */
    int32_t tid;  /*!< the zero-based thread id */
    tmap_map_opt_t *opt;  /*!< the options to this program */    
} tmap_map_all_thread_data_t;

/*!
 The core worker routine of mapall
 @param  seq_buffer         the buffer of sequences
 @param  sams               the sams to return
 @param  seq_buffer_length  the number of sequences in the buffer
 @param  refseq             the reference sequence
 @param  bwt                the BWT indices (forward/reverse)
 @param  sa                 the SA (forward/reverse)
 @param  thread_block_size  the number of reads to batch per thread
 @param  tid                the thread ids
 @param  opt                the program parameters 
 */
void
tmap_map_all_core_worker(tmap_seq_t **seq_buffer, tmap_map_sams_t **sams, int32_t seq_buffer_length,
                         tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2],
                         int32_t thread_block_size, int32_t tid, tmap_map_opt_t *opt);

/*!
 A wrapper around the core function of mapall
 @param  arg  the worker arguments in the type: tmap_map_all_thread_data_t
 @return      the worker arguments
 */
void *
tmap_map_all_core_thread_worker(void *arg);

/*!
  Parses the command line options and stores them in the options structure
  @param  argc  the number of arguments
  @param  argv  the argument list
  @param  opt   pointer to the options
  @return       1 if successful, 0 otherwise
  */
int32_t
tmap_map_all_opt_parse(int argc, char *argv[], tmap_map_opt_t *opt);

/*! 
  main-like function for 'tmap map_all'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
tmap_map_all_main(int argc, char *argv[]);

#endif 
