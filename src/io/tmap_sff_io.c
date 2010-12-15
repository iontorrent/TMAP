/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>

#include "../util/tmap_alloc.h"
#include "../seq/tmap_sff.h"
#include "tmap_file.h"
#include "tmap_sff_io.h"

inline tmap_sff_io_t *
tmap_sff_io_init(tmap_file_t *fp)
{
  tmap_sff_io_t *sffio = NULL;

  sffio = tmap_calloc(1, sizeof(tmap_sff_io_t), "sffio");

  sffio->fp = fp;
  sffio->gheader = tmap_sff_header_read(sffio->fp);
  sffio->n_read = 0;

  return sffio;
}

void
tmap_sff_io_destroy(tmap_sff_io_t *sffio)
{
  tmap_sff_header_destroy(sffio->gheader);
  free(sffio);
}

int32_t
tmap_sff_io_read(tmap_sff_io_t *sffio, tmap_sff_t *sff)
{
  // we have read them all
  if(sffio->gheader->n_reads <= sffio->n_read) return EOF;

  // destroy previous data, if any
  if(NULL != sff->rheader) {
      tmap_sff_read_header_destroy(sff->rheader);
      sff->rheader = NULL;
  }
  if(NULL != sff->read) {
      tmap_sff_read_destroy(sff->read);
      sff->read = NULL;
  }

  sff->gheader = sffio->gheader;
  sff->rheader = tmap_sff_read_header_read(sffio->fp);
  sff->read = tmap_sff_read_read(sffio->fp, sffio->gheader, sff->rheader);

  sffio->n_read++;

  return 1;
}

int32_t
tmap_sff_io_read_buffer(tmap_sff_io_t *sffio, tmap_sff_t **sff_buffer, int32_t buffer_length)
{
  int32_t n = 0;

  if(buffer_length <= 0) return 0;

  while(n < buffer_length && 0 <= tmap_sff_io_read(sffio, sff_buffer[n])) {
      n++;
  }

  return n;
}
