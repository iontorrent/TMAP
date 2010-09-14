#ifndef FMAP_SEQ_IO_H_
#define FMAP_SEQ_IO_H_

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

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "../util/fmap_seq.h"
#include "fmap_file.h"

#define FMAP_STREAM_BUFFER_SIZE 4096

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/*! @typedef
  @abstract         
  @field  buf      the character buffer
  @field  begin    the index of the next character in the buffer
  @field  end      the number of characters last read
  @field  is_eof   1 if the EOF marker has been reached, 0 otherwise
  @field  bufsize  the size of the character buffer
  */
typedef struct {
    char *buf;
    int32_t begin, end, is_eof;
    fmap_file_t *f;
    int32_t bufsize;
} fmap_stream_t; 

/*! @typedef
  @abstract         structure for reading FASTA/FASTQ strings
  @field  last_char  the last character read
  @field  f          pointer to the file structure
  */
typedef struct {
    int last_char;
    fmap_stream_t *f;
} fmap_seq_io_t;

/*! @function
  @abstract   initializes sequence reading structure
  @param  fp  a pointer to a file structure from which to read
  @return     pointer to the initialized memory for reading in sequences
  */
inline fmap_seq_io_t *
fmap_seq_io_init(fmap_file_t *fp);

/*! @function
  @abstract      destroys sequence reading structure
  @param  seqio  a pointer to the sequence structure
  */
inline void 
fmap_seq_io_destroy(fmap_seq_io_t *seqio);

/*! @function
  @abstract      reads in a reading structure
  @param  seqio  a pointer to a previously initialized sequence structure
  @param  seq    the sequence structure in which to store the data
  @return        the length of the sequence read, -1 indicates an a EOF, -2 indicates a truncated quality string
  */
int 
fmap_seq_io_read(fmap_seq_io_t *seqio, fmap_seq_t *seq);

/*! @function
  @abstract              reads sequences into a buffer
  @param  seqio          a pointer to a previously initialized sequence structure
  @param  seq_buffer     the sequence structure in which to store the data
  @param  buffer_length  the number of sequences to read
  @return                the number of sequences read
  */
int
fmap_seq_io_read_buffer(fmap_seq_io_t *seqio, fmap_seq_t **seq_buffer, int32_t buffer_length);

#endif
