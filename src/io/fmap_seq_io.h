#ifndef FMAP_SEQ_IO_H_
#define FMAP_SEQ_IO_H_

#include "../seq/fmap_fq.h"
#include "../seq/fmap_sff.h"
#include "../seq/fmap_seq.h"
#include "fmap_sff_io.h"
#include "fmap_fq_io.h"

/*! 
  An Abstract DNA Sequence Reading Library
  */

/*! 
*/
typedef struct {
  int8_t type;  /*!< the type of io associated with this structure */
  union {
  fmap_fq_io_t *fqio;  /*!< the pointer to the fastq io structure */
  fmap_sff_io_t *sffio;  /*!< the pointer to the sff io structure */
  } io;
} fmap_seq_io_t;
/*! 
  initializes input/output structure
  @param  fp    a pointer to a file structure from which to read
  @param  type  the type of io associated with this structure
  @return       pointer to the initialized memory for reading in sequences
  */
inline fmap_seq_io_t *
fmap_seq_io_init(fmap_file_t *fp, int8_t type);

/*! 
  destroys input/output structure
  @param  io  a pointer to the sequence structure
  */
inline void
fmap_seq_io_destroy(fmap_seq_io_t *io);

/*! 
  reads in a reading structure
  @param  io     a pointer to a previously initialized sequence structure
  @param  seq    the sequence structure in which to store the data
  @return        the length of the sequence read, -1 indicates an a EOF, -2 indicates a truncated quality string
  */
inline int
fmap_seq_io_read(fmap_seq_io_t *io, fmap_seq_t *seq);

/*! 
  reads sequences into a buffer
  @param  io             a pointer to a previously initialized sequence structure
  @param  seq_buffer     the sequence structure in which to store the data
  @param  buffer_length  the number of sequences to read
  @return                the number of sequences read
  */
int
fmap_seq_io_read_buffer(fmap_seq_io_t *io, fmap_seq_t **seq_buffer, int32_t buffer_length);

/*! 
  main-like function for 'fmap sff2fq'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int
fmap_seq_io_sff2fq_main(int argc, char *argv[]);

#endif
