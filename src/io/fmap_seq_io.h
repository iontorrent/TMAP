#ifndef FMAP_SEQ_IO_H_
#define FMAP_SEQ_IO_H_

#include "../seq/fmap_fq.h"
#include "../seq/fmap_sff.h"
#include "../seq/fmap_seq.h"
#include "fmap_sff_io.h"
#include "fmap_fq_io.h"

/*! @typedef 
  @field  type  the type of io associated with this structure
  @field  io    pointer to the particular io data structure
*/
typedef struct {
    int8_t type;
    union {
        fmap_fq_io_t *fqio;
        fmap_sff_io_t *sffio;
    } io;
} fmap_seq_io_t;

/*! @function
  @abstract     initializes input/output structure
  @param  fp    a pointer to a file structure from which to read
  @param  type  the type of io associated with this structure
  @return       pointer to the initialized memory for reading in sequences
  */
inline fmap_seq_io_t *
fmap_seq_io_init(fmap_file_t *fp, int8_t type);

/*! @function
  @abstract   destroys input/output structure
  @param  io  a pointer to the sequence structure
  */
inline void
fmap_seq_io_destroy(fmap_seq_io_t *io);

/*! @function
  @abstract      reads in a reading structure
  @param  io  a pointer to a previously initialized sequence structure
  @param  seq    the sequence structure in which to store the data
  @return        the length of the sequence read, -1 indicates an a EOF, -2 indicates a truncated quality string
  */
inline int
fmap_seq_io_read(fmap_seq_io_t *io, fmap_seq_t *seq);

/*! @function
  @abstract              reads sequences into a buffer
  @param  io             a pointer to a previously initialized sequence structure
  @param  seq_buffer     the sequence structure in which to store the data
  @param  buffer_length  the number of sequences to read
  @return                the number of sequences read
  */
int
fmap_seq_io_read_buffer(fmap_seq_io_t *io, fmap_seq_t **seq_buffer, int32_t buffer_length);

#endif
