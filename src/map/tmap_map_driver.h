/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_DRIVER_H
#define TMAP_MAP_DRIVER_H

#include <sys/types.h>

#ifdef HAVE_LIBPTHREAD
#define TMAP_MAP_DRIVER_THREAD_BLOCK_SIZE 512
#endif

// TODO
// init and cleanup script for each sequence

/*!
  This function will be invoked after reading in all the reference data
  to initialize any program options and print messages.
  @param  refseq  the reference sequence
  @param  opt     the program options
  @return         0 upon success, non-zero otherwise
 */
typedef int32_t (*tmap_driver_func_init)(tmap_refseq_t *refseq, tmap_map_opt_t *opt);

/*!
  This function will be invoked before a thread begins process its sequences.
  @param  data  the thread persistent data
  @param  opt   the program options
  @return       0 upon success, non-zero otherwise
 */
typedef int32_t (*tmap_driver_func_thread_init)(void **data, tmap_map_opt_t *opt);

/*!
  This function will be invoked to map a sequence.
  @param  data    the thread persistent data
  @param  seq     the sequence to map
  @param  refseq  the reference sequence
  @param  bwt     the bwt structure
  @param  sa      the sa structure
  @param  opt     the program options
  @return         the mappings upon success, NULL otherwise
 */
typedef tmap_map_sams_t* (*tmap_driver_func_thread_map)(void **data, tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2], tmap_map_opt_t *opt);

/*!
  This function will be invoked to give a mapping quality to a set of mappings
  @param  sams     the mappings
  @param  seq_len  the sequence length
  @param  opt      the program options
  @return          0 upon success, non-zero otherwise
  */
typedef int32_t (*tmap_driver_func_mapq)(tmap_map_sams_t *sams, int32_t seq_len, tmap_map_opt_t *opt);

/*!
  This function will be invoked after a thread has process all of its sequences.
  @param  data  the thread persistent data
  @param  opt   the program options
  @return       0 upon success, non-zero otherwise
  */
typedef int32_t (*tmap_driver_func_thread_cleanup)(void **data, tmap_map_opt_t *opt);

/*! 
  data to be passed to a thread                         
  */
typedef struct {                                            
    tmap_seq_t **seq_buffer;  /*!< the buffer of sequences */    
    int32_t seq_buffer_length;  /*!< the buffer length */
    tmap_map_sams_t **sams;  /*!< the alignments for each sequence */
    tmap_refseq_t *refseq;  /*!< pointer to the reference sequence (forward) */
    tmap_bwt_t *bwt[2];  /*!< pointer to the BWT indices (forward/reverse) */
    tmap_sa_t *sa[2];  /*!< pointer to the SA (forward/reverse) */    
    tmap_driver_func_thread_init func_thread_init; /* this function will be run once per thread to initialize persistent data across that thread */
    tmap_driver_func_thread_map func_thread_map; /* this function will be run once per thread per input sequence to map the sequence */
    tmap_driver_func_mapq func_mapq; /* this function will be run to calculate the mapping quality */
    tmap_driver_func_thread_cleanup func_thread_cleanup; /* this function will be run once per thread to cleanup/destroy any persistent data across that thread */
    int32_t thread_block_size; /*!< the number of reads per thread to process */
    int32_t tid;  /*!< the zero-based thread id */
    tmap_map_opt_t *opt;  /*!< the options to this program */    
} tmap_map_driver_thread_data_t;

/*!
  The core worker routine of mapall
  @param  seq_buffer           the buffer of sequences
  @param  sams                 the sams to return
  @param  seq_buffer_length    the number of sequences in the buffer
  @param  refseq               the reference sequence
  @param  bwt                  the BWT indices (forward/reverse)
  @param  sa                   the SA (forward/reverse)
  @param  func_thread_init     the thread initialization function
  @param  func_thread_map      the thread map function
  @param  func_mapq            the mapping quality function
  @param  func_thread_cleanup  the thread cleanup function
  @param  thread_block_size    the number of reads to batch per thread
  @param  tid                  the thread ids
  @param  opt                  the program parameters 
 */
void
tmap_map_driver_core_worker(tmap_seq_t **seq_buffer, tmap_map_sams_t **sams, int32_t seq_buffer_length,
                         tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2],
                         tmap_driver_func_thread_init func_thread_init, 
                         tmap_driver_func_thread_map func_thread_map, 
                         tmap_driver_func_mapq func_mapq,
                         tmap_driver_func_thread_cleanup func_thread_cleanup,
                         int32_t thread_block_size, int32_t tid, tmap_map_opt_t *opt);

/*!
 A wrapper around the core function of mapall
 @param  arg  the worker arguments in the type: tmap_map_driver_thread_data_t
 @return      the worker arguments
 */
void *
tmap_map_driver_core_thread_worker(void *arg);

/*!
  the core routine for mapping data with or without threads
  @param  func_init            the mapping algorithm initialization function
  @param  func_thread_init     the thread initialization function
  @param  func_thread_map      the thread map function
  @param  func_mapq            the mapq function
  @param  func_thread_cleanup  the thread cleanup function
  @param  opt                  the program parameters 
  */
void
tmap_map_driver_core(tmap_driver_func_init fund_init,
                  tmap_driver_func_thread_init func_thread_init, 
                  tmap_driver_func_thread_map func_thread_map, 
                  tmap_driver_func_mapq func_mapq,
                  tmap_driver_func_thread_cleanup func_thread_cleanup,
                  tmap_map_opt_t *opt);

#endif 
