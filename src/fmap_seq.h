#ifndef FMAP_SEQ_H_
#define FMAP_SEQ_H_

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

/* Nils Homer - modified not to macro-ized */

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "fmap_io.h"

#define FMAP_STREAM_BUFFER_SIZE 4096

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/*! @typedef
  @abstract
  @field  buf       buffer for the stream
  @filed  begin
  @field  end
  @field  is_eof    reached EOF
  @field  buf_size  buffer size
  @field  f         file pointer
  */
typedef struct {
    char *buf;
    int begin, end, is_eof, bufsize;
    fmap_file_t *f;
} fmap_stream_t;

/*! @typedef
  @abstract
  @field  l  the length of the string
  @field  m  the memory allocated for this string
  @field  s  the pointer to the string
  */
typedef struct {
    size_t l, m;
    char *s;
} fmap_string_t;

/*! @typedef
  @abstract         structure for holding FASTA/FASTQ strings
  @field  name       the name string
  @field  comment    the comment string
  @field  seq        the sequence string
  @field  qual       the quality string
  @field  last_char  the last character read
  @field  f          pointer to the file structure
  */
typedef struct {
    fmap_string_t name, comment, seq, qual;
    int last_char;
    fmap_stream_t *f;
} fmap_seq_t;

/*! @function
  @abstract   initializes sequence reading structure
  @param  fp  a pointer to a file structure from which to read
  @return     pointer to the initialized memory for reading in sequences
  */
inline fmap_seq_t *
fmap_seq_init(fmap_file_t *fp);

/*! @function
  @abstract   initializes sequence reading structure
  @param  fp  a pointer to a file structure from which to read
  @return     pointer to the initialized memory for reading in sequences
  */
inline void 
fmap_seq_destroy(fmap_seq_t *seq);

/*! @function
  @abstract    reads in a reading structure
  @param  seq  a pointer to a previously initialized sequence structure
  @return      the length of the sequence read, -1 indicates an a EOF, -2 indicates a truncated quality string
  */
int 
fmap_seq_read(fmap_seq_t *seq);

#endif
