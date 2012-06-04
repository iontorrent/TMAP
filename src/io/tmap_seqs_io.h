/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SEQS_IO_H
#define TMAP_SEQS_IO_H

#include <config.h>
#include "../seq/tmap_fq.h"
#include "../seq/tmap_sff.h"
#include "../seq/tmap_seq.h"
#include "../seq/tmap_seqs.h"
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
  tmap_seq_io_t **seqios; // TODO
  int32_t n; // TODO
} tmap_seqs_io_t;

/*! 
  initializes input/output structure
  @param  fns           the file name(s) of the input
  @param  fn_num       the number of file names
  @param  seq_type     the type of io associated with this structure
  @param  compression  the compression type
  @return              pointer to the initialized memory for reading/writing sequences
  */
inline tmap_seqs_io_t *
tmap_seqs_io_init(const char **fns, int32_t fn_num, int8_t seq_type, int32_t compression);

/*! 
  destroys input/output structure
  @param  io  a pointer to the sequence structure
  */
inline void
tmap_seqs_io_destroy(tmap_seqs_io_t *io);

/*! 
  reads in a reading structure
  @param  io     a pointer to a previously initialized sequence structure
  @param  seqs   the sequence structure in which to store the data
  @return        non-negative if successfull, -1 indicates an a EOF, -2 indicates a truncated quality string
  */
inline int
tmap_seqs_io_read(tmap_seqs_io_t *io, tmap_seqs_t *seqs);

/*! 
  reads sequences into a buffer
  @param  io             a pointer to a previously initialized sequence structure
  @param  seqs_buffer    the sequence structure in which to store the data
  @param  buffer_length  the number of sequences to read
  @return                the number of sequences read
  */
int
tmap_seqs_io_read_buffer(tmap_seqs_io_t *io, tmap_seqs_t **seqs_buffer, int32_t buffer_length);


#endif
