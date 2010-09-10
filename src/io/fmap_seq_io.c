#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_seq.h"
#include "fmap_io.h"
#include "fmap_seq_io.h"

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

inline fmap_seq_io_t *
fmap_seq_io_init(fmap_file_t *fp)
{
  fmap_seq_io_t *s = fmap_calloc(1, sizeof(fmap_seq_io_t), "s");
  s->f = fmap_stream_init(fp, FMAP_STREAM_BUFFER_SIZE);
  return s;
}

static inline void 
fmap_seq_io_rewind(fmap_seq_io_t *seq)
{
  seq->last_char = 0;
  seq->f->is_eof = seq->f->begin = seq->f->end = 0;
}

inline void 
fmap_seq_io_destroy(fmap_seq_io_t *seqio)
{
  if(NULL == seqio) return;
  fmap_stream_destroy(seqio->f);
  free(seqio);
}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
   */
int 
fmap_seq_io_read(fmap_seq_io_t *seqio, fmap_seq_t *seq)
{
  int c;
  fmap_stream_t *ks = seqio->f;
  if (seqio->last_char == 0) { /* then jump to the next header line */
      while ((c = fmap_stream_getc(ks)) != -1 && c != '>' && c != '@');
      if (c == -1) return -1; /* end of file */
      seqio->last_char = c;
  } /* the first header char has been read */
  seq->comment.l = seq->seq.l = seq->qual.l = 0;
  if (fmap_stream_getuntil(ks, 0, &seq->name, &c) < 0) return -1;
  if (c != '\n') fmap_stream_getuntil(ks, '\n', &seq->comment, 0);
  while ((c = fmap_stream_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
      if (isgraph(c)) { /* printable non-space character */
          if (seq->seq.l + 1 >= seq->seq.m) { /* double the memory */
              seq->seq.m = seq->seq.l + 2;
              kroundup32(seq->seq.m); /* rounded to next closest 2^k */
              seq->seq.s = fmap_realloc(seq->seq.s, seq->seq.m, "seq->seq.s");
          }
          seq->seq.s[seq->seq.l++] = (char)c;
      }
  }
  if (c == '>' || c == '@') seqio->last_char = c; /* the first header char has been read */
  seq->seq.s[seq->seq.l] = 0;	/* null terminated string */
  if (c != '+') return seq->seq.l; /* FASTA */
  if (seq->qual.m < seq->seq.m) {	/* allocate enough memory */
      seq->qual.m = seq->seq.m;
      seq->qual.s = fmap_realloc(seq->qual.s, seq->qual.m, "seq->qual.s");
  }
  while ((c = fmap_stream_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */
  if (c == -1) return -2; /* we should not stop here */
  while ((c = fmap_stream_getc(ks)) != -1 && seq->qual.l < seq->seq.l)
    if (c >= 33 && c <= 127) seq->qual.s[seq->qual.l++] = (unsigned char)c;
  seq->qual.s[seq->qual.l] = 0; /* null terminated string */
  seqio->last_char = 0;	/* we have not come to the next header line */
  if (seq->seq.l != seq->qual.l) return -2; /* qual string is shorter than seq string */
  return seq->seq.l;
}

int
fmap_seq_io_read_buffer(fmap_seq_io_t *seqio, fmap_seq_t **seq_buffer, int32_t buffer_length)
{
  int32_t n = 0;
  
  if(buffer_length <= 0) return 0;

  while(n < buffer_length && 0 <= fmap_seq_io_read(seqio, seq_buffer[n])) {
      n++;
  }

  return n;
}
