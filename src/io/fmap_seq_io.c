#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_progress.h"
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

static int
fmap_seq_io_print(fmap_file_t *fp, fmap_seq_t *seq)
{
  switch(seq->type) {
    case FMAP_SEQ_TYPE_FQ:
      return fmap_file_fprintf(fp, "@%s\n%s\n+%s\n%s\n",
                        seq->data.fq->name->s,
                        seq->data.fq->seq->s,
                        (0 < seq->data.fq->comment->l) ? seq->data.fq->comment->s : "",
                        seq->data.fq->qual->s);
      break;
    case FMAP_SEQ_TYPE_SFF:
      fmap_error("SFF writing is unsupported", Exit, OutOfRange);
      break;
    default:
      fmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return 0;
}

int
fmap_seq_io_sff2fq_main(int argc, char *argv[])
{
  int c, help = 0;
  fmap_file_t *fmap_file_in = NULL;
  fmap_seq_io_t *io_in = NULL, *io_out = NULL;
  fmap_seq_t *seq_in = NULL, *seq_out = NULL;

  while((c = getopt(argc, argv, "vh")) >= 0) {
      switch(c) {
        case 'v': fmap_progress_set_verbosity(1); break;
        case 'h': help = 1; break;
        default: return 1;
      }
  }
  if(1 != argc - optind || 1 == help) {
      fmap_file_fprintf(fmap_file_stderr, "Usage: %s %s [-v -h] <in.sff>\n", PACKAGE, argv[0]);
      return 1;
  }

  // input
  fmap_file_in = fmap_file_fopen(argv[optind], "rb", FMAP_FILE_NO_COMPRESSION);
  io_in = fmap_seq_io_init(fmap_file_in, FMAP_SEQ_TYPE_SFF);
  seq_in = fmap_seq_init(FMAP_SEQ_TYPE_SFF);

  // output
  fmap_file_stdout = fmap_file_fdopen(fileno(stdout), "wb", FMAP_FILE_NO_COMPRESSION);
  io_out = fmap_seq_io_init(fmap_file_stdout, FMAP_SEQ_TYPE_FQ);

  while(0 < fmap_seq_io_read(io_in, seq_in)) {
      seq_out = fmap_seq_sff2fq(seq_in);
      fmap_seq_io_print(fmap_file_stdout, seq_out);
      fmap_seq_destroy(seq_out);
  }
  fmap_seq_destroy(seq_in);

  // input
  fmap_seq_io_destroy(io_in);
  fmap_file_fclose(fmap_file_in);
  
  // output
  fmap_seq_io_destroy(io_out);
  fmap_file_fclose(fmap_file_stdout);

  return 0;
}
