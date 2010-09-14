#include <stdlib.h>
#include <stdio.h>

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "fmap_seq_io.h"
#include "fmap_sff_io.h"
#include "fmap_seq_io.h"

inline fmap_seq_io_t *
fmap_seq_io_init(fmap_file_t *fp, int8_t type)
{
  fmap_seq_io_t *io = NULL;

  io = fmap_calloc(1, sizeof(fmap_seq_io_t), "io");
  io->type = type;

  switch(io->type) {
    case FMAP_SEQ_TYPE_FQ:
      io->io.fqio = fmap_fq_io_init(fp);
      break;
    case FMAP_SEQ_TYPE_SFF:
      io->io.sffio = fmap_sff_io_init(fp);
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }

  return io;
}

inline void
fmap_seq_io_destroy(fmap_seq_io_t *io)
{
  switch(io->type) {
    case FMAP_SEQ_TYPE_FQ:
      fmap_fq_io_destroy(io->io.fqio);
      break;
    case FMAP_SEQ_TYPE_SFF:
      fmap_sff_io_destroy(io->io.sffio);
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }

  free(io);
}

inline int
fmap_seq_io_read(fmap_seq_io_t *io, fmap_seq_t *seq)
{
  if(io->type != seq->type) {
      fmap_error("type mismatch", Exit, OutOfRange);
  }
  switch(io->type) {
    case FMAP_SEQ_TYPE_FQ:
      seq->type = io->type;
      return fmap_fq_io_read(io->io.fqio, seq->data.fq);
      break;
    case FMAP_SEQ_TYPE_SFF:
      seq->type = io->type;
      return fmap_sff_io_read(io->io.sffio, seq->data.sff);
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }

  return EOF;
}

int
fmap_seq_io_read_buffer(fmap_seq_io_t *io, fmap_seq_t **seq_buffer, int32_t buffer_length)
{
  int32_t n = 0;

  if(buffer_length <= 0) return 0;

  while(n < buffer_length) {
      if(NULL == seq_buffer[n]) {
          seq_buffer[n] = fmap_seq_init(io->type);
      }
      else if(io->type != seq_buffer[n]->type) { // check if we need to change the type
          fmap_seq_destroy(seq_buffer[n]);
          seq_buffer[n] = fmap_seq_init(io->type);
      }
      if(fmap_seq_io_read(io, seq_buffer[n]) < 0) {
          break;
      }
      n++;
  }

  return n;
}
