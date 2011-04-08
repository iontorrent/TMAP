/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "tmap_seq_io.h"
#include "tmap_sff_io.h"
#include "tmap_seq_io.h"

static inline tmap_file_t *
tmap_seq_io_open(const char *fn, int32_t out_type, int32_t compression)
{
  if(0 == out_type) {
      if(NULL == fn || 0 == strcmp("-", fn)) {
          return tmap_file_fdopen(fileno(stdin), "rb", compression);
      }
      else {
          return tmap_file_fopen(fn, "rb", compression);
      }
  }
  else {
      if(NULL == fn || 0 == strcmp("-", fn)) {
          return tmap_file_fdopen(fileno(stdout), "wb", compression);
      }
      else {
          return tmap_file_fopen(fn, "wb", compression);
      }
  }
}

inline tmap_seq_io_t *
tmap_seq_io_init(const char *fn, int8_t seq_type, int32_t out_type, int32_t compression)
{
  tmap_seq_io_t *io = NULL;

  io = tmap_calloc(1, sizeof(tmap_seq_io_t), "io");
  io->type = seq_type;

#ifdef HAVE_SAMTOOLS
  if(TMAP_SEQ_TYPE_SAM == io->type || TMAP_SEQ_TYPE_BAM == io->type) {
      if(1 == out_type) {
          tmap_error("Writing not supported with a SAM/BAM file", Exit, OutOfRange);
      }
      if(TMAP_FILE_NO_COMPRESSION != compression) {
          tmap_error("Compression not supported with a SAM/BAM file", Exit, OutOfRange);
      }
  }
#endif

  switch(io->type) {
    case TMAP_SEQ_TYPE_FQ:
      io->fp = tmap_seq_io_open(fn, out_type, compression);
      io->io.fqio = tmap_fq_io_init(io->fp);
      break;
    case TMAP_SEQ_TYPE_SFF:
      io->fp = tmap_seq_io_open(fn, out_type, compression);
      io->io.sffio = tmap_sff_io_init(io->fp);
      break;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
      io->fp = NULL;
      io->io.samio = tmap_sam_io_init(fn);
      break;
    case TMAP_SEQ_TYPE_BAM:
      io->fp = NULL;
      io->io.samio = tmap_bam_io_init(fn);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }

  return io;
}

inline void
tmap_seq_io_destroy(tmap_seq_io_t *io)
{
  switch(io->type) {
    case TMAP_SEQ_TYPE_FQ:
      tmap_file_fclose(io->fp);
      tmap_fq_io_destroy(io->io.fqio);
      break;
    case TMAP_SEQ_TYPE_SFF:
      tmap_file_fclose(io->fp);
      tmap_sff_io_destroy(io->io.sffio);
      break;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_sam_io_destroy(io->io.samio);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }

  free(io);
}

inline int
tmap_seq_io_read(tmap_seq_io_t *io, tmap_seq_t *seq)
{
  if(io->type != seq->type) {
      tmap_error("type mismatch", Exit, OutOfRange);
  }
  switch(io->type) {
    case TMAP_SEQ_TYPE_FQ:
      seq->type = io->type;
      return tmap_fq_io_read(io->io.fqio, seq->data.fq);
      break;
    case TMAP_SEQ_TYPE_SFF:
      seq->type = io->type;
      return tmap_sff_io_read(io->io.sffio, seq->data.sff);
      break;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      seq->type = io->type;
      return tmap_sam_io_read(io->io.samio, seq->data.sam);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }

  return EOF;
}

int
tmap_seq_io_read_buffer(tmap_seq_io_t *io, tmap_seq_t **seq_buffer, int32_t buffer_length)
{
  int32_t n = 0;

  if(buffer_length <= 0) return 0;

  while(n < buffer_length) {
      if(NULL == seq_buffer[n]) {
          seq_buffer[n] = tmap_seq_init(io->type);
      }
      else if(io->type != seq_buffer[n]->type) { // check if we need to change the type
          tmap_seq_destroy(seq_buffer[n]);
          seq_buffer[n] = tmap_seq_init(io->type);
      }
      if(tmap_seq_io_read(io, seq_buffer[n]) < 0) {
          break;
      }
      n++;
  }

  return n;
}

static inline int
tmap_seq_io_print2(tmap_file_t *fp, tmap_seq_t *seq)
{
  switch(seq->type) {
    case TMAP_SEQ_TYPE_FQ:
      return tmap_file_fprintf(fp, "@%s\n%s\n+%s\n%s\n",
                        seq->data.fq->name->s,
                        seq->data.fq->seq->s,
                        (0 < seq->data.fq->comment->l) ? seq->data.fq->comment->s : "",
                        seq->data.fq->qual->s);
      break;
    case TMAP_SEQ_TYPE_SFF:
      tmap_error("SFF writing is unsupported", Exit, OutOfRange);
      break;
#ifdef HAVE_SAMTOOLS
    case TMAP_SEQ_TYPE_SAM:
    case TMAP_SEQ_TYPE_BAM:
      tmap_error("SAM/BAM writing is unsupported", Exit, OutOfRange);
      break;
#endif
    default:
      tmap_error("type is unrecognized", Exit, OutOfRange);
      break;
  }
  return 0;
}

static int
tmap_seq_io_print(tmap_seq_io_t *io, tmap_seq_t *seq)
{
  return tmap_seq_io_print2(io->fp, seq);
}

int
tmap_seq_io_sff2fq_main(int argc, char *argv[])
{
  int c, help = 0;
  tmap_seq_io_t *io_in = NULL, *io_out = NULL;
  tmap_seq_t *seq_in = NULL, *seq_out = NULL;

  while((c = getopt(argc, argv, "vh")) >= 0) {
      switch(c) {
        case 'v': tmap_progress_set_verbosity(1); break;
        case 'h': help = 1; break;
        default: return 1;
      }
  }
  if(1 != argc - optind || 1 == help) {
      tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [-v -h] <in.sff>\n", PACKAGE, argv[0]);
      return 1;
  }

  // input
  io_in = tmap_seq_io_init(argv[optind], TMAP_SEQ_TYPE_SFF, 0, TMAP_FILE_NO_COMPRESSION);
  seq_in = tmap_seq_init(TMAP_SEQ_TYPE_SFF);

  // output
  io_out = tmap_seq_io_init("-", TMAP_SEQ_TYPE_FQ, 1, TMAP_FILE_NO_COMPRESSION);

  while(0 < tmap_seq_io_read(io_in, seq_in)) {
      seq_out = tmap_seq_sff2fq(seq_in);
      tmap_seq_io_print(io_out, seq_out);
      tmap_seq_destroy(seq_out);
  }
  tmap_seq_destroy(seq_in);

  // input
  tmap_seq_io_destroy(io_in);
  
  // output
  tmap_seq_io_destroy(io_out);

  return 0;
}
