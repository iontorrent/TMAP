/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SEQ_IO_H
#define TMAP_SEQ_IO_H

#include <config.h>
#include "../seq/tmap_fq.h"
#include "../seq/tmap_sff.h"
#include "../seq/tmap_seq.h"
#include "../samtools/sam_header.h"
#include "tmap_fq_io.h"
#include "tmap_sff_io.h"
#include "tmap_sam_io.h"

/*! 
  An Abstract DNA Sequence Reading Library
  */

/*! 
*/
typedef struct {
  int8_t type;  /*!< the type of io associated with this structure */
  tmap_file_t *fp; /*!< the file pointer */
  union {
      tmap_fq_io_t *fqio;  /*!< the pointer to the fastq io structure */
      tmap_sff_io_t *sffio;  /*!< the pointer to the sff io structure */
      tmap_sam_io_t *samio;  /*!< the pointer to the SAM/BAM io structure */
  } io;
} tmap_seq_io_t;

/*! 
  initializes input/output structure
  @param  fn           the file name of the input/output
  @param  seq_type     the type of io associated with this structure
  @param  out_type     the output type (0 for reading, 1 for writing)
  @param  compression  the compression type
  @return              pointer to the initialized memory for reading/writing sequences
  */
inline tmap_seq_io_t *
tmap_seq_io_init(const char *fn, int8_t seq_type, int32_t out_type, int32_t compression);

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

/*!
  @param  io  a pointer to a previously initialized sequence structure
  @return     SAM header structure 
 */
sam_header_t*
tmap_seq_io_get_sam_header(tmap_seq_io_t *io);

#endif
