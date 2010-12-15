/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SEQ_IO_H
#define TMAP_SEQ_IO_H

#include "../seq/tmap_fq.h"
#include "../seq/tmap_sff.h"
#include "../seq/tmap_seq.h"
#include "tmap_sff_io.h"
#include "tmap_fq_io.h"

/*! 
  An Abstract DNA Sequence Reading Library
  */

/*! 
*/
typedef struct {
  int8_t type;  /*!< the type of io associated with this structure */
  union {
  tmap_fq_io_t *fqio;  /*!< the pointer to the fastq io structure */
  tmap_sff_io_t *sffio;  /*!< the pointer to the sff io structure */
  } io;
} tmap_seq_io_t;

/*! 
  initializes input/output structure
  @param  fp    a pointer to a file structure from which to read
  @param  type  the type of io associated with this structure
  @return       pointer to the initialized memory for reading in sequences
  */
inline tmap_seq_io_t *
tmap_seq_io_init(tmap_file_t *fp, int8_t type);

/*! 
  destroys input/output structure
  @param  io  a pointer to the sequence structure
  */
inline void
tmap_seq_io_destroy(tmap_seq_io_t *io);

/*! 
  reads in a reading structure
  @param  io     a pointer to a previously initialized sequence structure
  @param  seq    the sequence structure in which to store the data
  @return        the length of the sequence read, -1 indicates an a EOF, -2 indicates a truncated quality string
  */
inline int
tmap_seq_io_read(tmap_seq_io_t *io, tmap_seq_t *seq);

/*! 
  reads sequences into a buffer
  @param  io             a pointer to a previously initialized sequence structure
  @param  seq_buffer     the sequence structure in which to store the data
  @param  buffer_length  the number of sequences to read
  @return                the number of sequences read
  */
int
tmap_seq_io_read_buffer(tmap_seq_io_t *io, tmap_seq_t **seq_buffer, int32_t buffer_length);

/*! 
  main-like function for 'tmap sff2fq'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
tmap_seq_io_sff2fq_main(int argc, char *argv[]);

#endif
