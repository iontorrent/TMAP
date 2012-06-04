/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sam_print.h"
#include "../samtools/sam_header.h"
#include "tmap_seq_io.h"
#include "tmap_sff_io.h"
#include "tmap_seq_io.h"
#include "tmap_seqs_io.h"

inline tmap_seqs_io_t*
tmap_seqs_io_init(char **fns, int32_t fn_num, int8_t seq_type, int32_t compression)
{
  tmap_seqs_io_t *io= NULL;
  int32_t i;

  io = tmap_calloc(1, sizeof(tmap_seqs_io_t), "io");
  io->type = seq_type;
      
  if(1 < io->n && (TMAP_SEQ_TYPE_SAM == io->type || TMAP_SEQ_TYPE_BAM == io->type)) {
      tmap_error("Multi-SAM/BAM not supported", Exit, OutOfRange);
  }

  if(NULL == fns) { // stdin
      io->n = 1;
      io->seqios = tmap_calloc(1, sizeof(tmap_seq_io_t*), "io->seqios");
      io->seqios[0] = tmap_seq_io_init("-", seq_type, 0, compression); // NB: always reading
  }
  else { // from file(s)
      io->n = fn_num;
      io->seqios = tmap_calloc(fn_num, sizeof(tmap_seq_io_t*), "io->seqios");
      for(i=0;i<io->n;i++) {
          io->seqios[i] = tmap_seq_io_init(fns[i], seq_type, 0, compression); // NB: always reading
      }
  }

  return io;
}

inline void
tmap_seqs_io_destroy(tmap_seqs_io_t *io)
{
  int32_t i;
  for(i=0;i<io->n;i++) {
      tmap_seq_io_destroy(io->seqios[i]);
  }
  free(io->seqios);
  free(io);
}

inline int
tmap_seqs_io_read(tmap_seqs_io_t *io, tmap_seqs_t *seqs)
{
  int32_t i;

  /*
   * Case 1 - SAM/BAM
   *    - NB: there must only be one input file
   *    - Read a record, if paired, then read the next
   * Case 2 - SFF/FQ
   *    - NB: there can be zero or more input files
   *    - Read one from each file, store in one record
   */

  if(io->type != seqs->type) {
      tmap_error("type mismatch", Exit, OutOfRange);
  }

  // reset seqs
  seqs->n = 0;
  if(TMAP_SEQ_TYPE_SAM == io->type || TMAP_SEQ_TYPE_BAM == io->type) {
      // NB: to supported paired reads, we check the paired flag
      for(i=0;i<2;i++) {
          tmap_seq_t *seq = tmap_seqs_get(seqs, i);
          if(tmap_seq_io_read(io->seqios[0], seq) < 0) return EOF; // TODO: better error checking
          tmap_seqs_add(seqs, seq); 
          // break if not paired
          if(0 == (seq->data.sam->b->core.flag & BAM_FPAIRED)) break;
      }
  }
  else {
      // read in one per file
      for(i=0;i<io->n;i++) {
          tmap_seq_t *seq = tmap_seqs_get(seqs, i);
          if(tmap_seq_io_read(io->seqios[i], seq) < 0) return EOF; // TODO: better error checking
          tmap_seqs_add(seqs, seq); 
      }
  }

  return 0;
}

int
tmap_seqs_io_read_buffer(tmap_seqs_io_t *io, tmap_seqs_t **seqs_buffer, int32_t buffer_length)
{
  int32_t n = 0;

  if(buffer_length <= 0) return 0;

  while(n < buffer_length) {
      if(NULL == seqs_buffer[n]) {
          seqs_buffer[n] = tmap_seqs_init(io->type);
      }
      else if(io->type != seqs_buffer[n]->type) { // check if we need to change the type
          tmap_seqs_destroy(seqs_buffer[n]);
          seqs_buffer[n] = tmap_seqs_init(io->type);
      }
      if(tmap_seqs_io_read(io, seqs_buffer[n]) < 0) {
          break;
      }
      n++;
  }

  return n;
}
