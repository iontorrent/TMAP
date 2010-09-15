#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_fq.h"
#include "fmap_file.h"
#include "fmap_fq_io.h"

/* The MIT License

   Copyright (c) 2008, by Heng Li <lh3@sanger.ac.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

*/

/* Nils Homer - modified not to be macro-ized */

static inline fmap_stream_t *
fmap_stream_init(fmap_file_t *f, int32_t bufsize)
{
  fmap_stream_t *ks = fmap_calloc(1, sizeof(fmap_stream_t), "ks");
  ks->f = f;
  ks->bufsize = bufsize;
  ks->buf = fmap_malloc(sizeof(char)*ks->bufsize, "ks->buf");
  return ks;
}

static inline void 
fmap_stream_destroy(fmap_stream_t *ks)
{
  if(NULL != ks) {
      free(ks->buf);
      free(ks);
  }
}

#define fmap_stream_eof(ks) \
  ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define fmap_stream_rewind(ks) \
  ((ks)->is_eof = (ks)->begin = (ks)->end = 0)

static inline int 
fmap_stream_getc(fmap_stream_t *ks)
{
  if (ks->is_eof && ks->begin >= ks->end) return -1;
  if (ks->begin >= ks->end) {
      ks->begin = 0;
      ks->end = fmap_file_fread2(ks->f, ks->buf, ks->bufsize);
      if (ks->end < ks->bufsize) ks->is_eof = 1;
      if (ks->end == 0) return -1;
  }
  return (int)ks->buf[ks->begin++];
}

static int 
fmap_stream_getuntil(fmap_stream_t *ks, int delimiter, fmap_string_t *str, int *dret)
{
  if (dret) *dret = 0;
  str->l = 0;
  if (ks->begin >= ks->end && ks->is_eof) return -1;
  for (;;) {
      int i;
      if (ks->begin >= ks->end) {
          if (!ks->is_eof) {
              ks->begin = 0;
              ks->end = fmap_file_fread2(ks->f, ks->buf, ks->bufsize);
              if (ks->end < ks->bufsize) ks->is_eof = 1;
              if (ks->end == 0) break;
          } else break;
      }
      if (delimiter) {
          for (i = ks->begin; i < ks->end; ++i)
            if (ks->buf[i] == delimiter) break;
      } else {
          for (i = ks->begin; i < ks->end; ++i)
            if (isspace(ks->buf[i])) break;
      }
      if (str->m - str->l < i - ks->begin + 1) {
          str->m = str->l + (i - ks->begin) + 1;
          kroundup32(str->m);
          str->s = fmap_realloc(str->s, str->m, "str->s");
      }
      memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin);
      str->l = str->l + (i - ks->begin);
      ks->begin = i + 1;
      if (i < ks->end) {
          if (dret) *dret = ks->buf[i];
          break;
      }
  }
  str->s[str->l] = '\0';
  return str->l;
}

inline fmap_fq_io_t *
fmap_fq_io_init(fmap_file_t *fp)
{
  fmap_fq_io_t *s = fmap_calloc(1, sizeof(fmap_fq_io_t), "s");
  s->f = fmap_stream_init(fp, FMAP_STREAM_BUFFER_SIZE);
  return s;
}

static inline void 
fmap_fq_io_rewind(fmap_fq_io_t *fq)
{
  fq->last_char = 0;
  fq->f->is_eof = fq->f->begin = fq->f->end = 0;
}

inline void 
fmap_fq_io_destroy(fmap_fq_io_t *fqio)
{
  if(NULL == fqio) return;
  fmap_stream_destroy(fqio->f);
  free(fqio);
}

/* Return value:
   >=0  length of the fastq (normal)
   -1   end-of-file
   -2   truncated quality string
   */
int 
fmap_fq_io_read(fmap_fq_io_t *fqio, fmap_fq_t *fq)
{
  int c;
  fmap_stream_t *ks = fqio->f;
  if (fqio->last_char == 0) { /* then jump to the next header line */
      while ((c = fmap_stream_getc(ks)) != -1 && c != '>' && c != '@');
      if (c == -1) return -1; /* end of file */
      fqio->last_char = c;
  } /* the first header char has been read */
  fq->comment->l = fq->seq->l = fq->qual->l = 0;
  if (fmap_stream_getuntil(ks, 0, fq->name, &c) < 0) return -1;
  if (c != '\n') fmap_stream_getuntil(ks, '\n', fq->comment, 0);
  while ((c = fmap_stream_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
      if (isgraph(c)) { /* printable non-space character */
          if (fq->seq->l + 1 >= fq->seq->m) { /* double the memory */
              fq->seq->m = fq->seq->l + 2;
              kroundup32(fq->seq->m); /* rounded to next closest 2^k */
              fq->seq->s = fmap_realloc(fq->seq->s, fq->seq->m, "fq->seq->s");
          }
          fq->seq->s[fq->seq->l++] = (char)c;
      }
  }
  if (c == '>' || c == '@') fqio->last_char = c; /* the first header char has been read */
  fq->seq->s[fq->seq->l] = 0;	/* null terminated string */
  if (c != '+') return fq->seq->l; /* FASTA */
  if (fq->qual->m < fq->seq->m) {	/* allocate enough memory */
      fq->qual->m = fq->seq->m;
      fq->qual->s = fmap_realloc(fq->qual->s, fq->qual->m, "fq->qual->s");
  }
  while ((c = fmap_stream_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */
  if (c == -1) return -2; /* we should not stop here */
  while ((c = fmap_stream_getc(ks)) != -1 && fq->qual->l < fq->seq->l)
    if (c >= 33 && c <= 127) fq->qual->s[fq->qual->l++] = (unsigned char)c;
  fq->qual->s[fq->qual->l] = 0; /* null terminated string */
  fqio->last_char = 0;	/* we have not come to the next header line */
  if (fq->seq->l != fq->qual->l) return -2; /* qual string is shorter than fq string */
  return fq->seq->l;
}

int
fmap_fq_io_read_buffer(fmap_fq_io_t *fqio, fmap_fq_t **fq_buffer, int32_t buffer_length)
{
  int32_t n = 0;
  
  if(buffer_length <= 0) return 0;

  while(n < buffer_length && 0 <= fmap_fq_io_read(fqio, fq_buffer[n])) {
      n++;
  }

  return n;
}
