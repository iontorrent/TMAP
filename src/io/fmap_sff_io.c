#include <stdlib.h>
#include <stdio.h>

#include "../util/fmap_alloc.h"
#include "../util/fmap_sff.h"
#include "fmap_file.h"
#include "fmap_sff_io.h"

inline fmap_sff_io_t *
fmap_sff_io_init(fmap_file_t *fp)
{
  fmap_sff_io_t *sffio = NULL;

  sffio = fmap_calloc(1, sizeof(fmap_sff_io_t), "sffio");

  sffio->fp = fp;
  sffio->gheader = fmap_sff_header_read(sffio->fp);
  sffio->n_read = 0;

  return sffio;
}

void
fmap_sff_io_destroy(fmap_sff_io_t *sffio)
{
  fmap_sff_header_destroy(sffio->gheader);
  free(sffio);
}

int32_t
fmap_sff_io_read(fmap_sff_io_t *sffio, fmap_sff_t *sff)
{
  // we have read them all
  if(sffio->gheader->n_reads <= sffio->n_read) return EOF;

  // destroy previous data, if any
  if(NULL != sff->rheader) {
      fmap_sff_read_header_destroy(sff->rheader);
      sff->rheader = NULL;
  }
  if(NULL != sff->read) {
      fmap_sff_read_destroy(sff->read);
      sff->read = NULL;
  }

  sff->gheader = sffio->gheader;
  sff->rheader = fmap_sff_read_header_read(sffio->fp);
  sff->read = fmap_sff_read_read(sffio->fp, sffio->gheader, sff->rheader);

  sffio->n_read++;

  return 1;
}

int32_t
fmap_sff_io_read_buffer(fmap_sff_io_t *sffio, fmap_sff_t **sff_buffer, int32_t buffer_length)
{
  int32_t n = 0;

  if(buffer_length <= 0) return 0;

  while(n < buffer_length && 0 <= fmap_sff_io_read(sffio, sff_buffer[n])) {
      n++;
  }

  return n;
}
