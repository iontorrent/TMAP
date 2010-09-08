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

/* Nils Homer - modified not to be macro-ized */

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

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
  */
typedef struct {
    fmap_string_t name, comment, seq, qual;
} fmap_seq_t;

/*! @function
  @abstract   initializes sequence read structure
  @return     pointer to the initialized memory 
  */
inline fmap_seq_t *
fmap_seq_init();

/*! @function
  @abstract   
  @param  seq  a pointer to the sequence structure
  */
inline void 
fmap_seq_destroy(fmap_seq_t *seq);

/*! @function
  @abstract         reverse the seq and qual fields
  @param  seq       a pointer to a sequence structure
  @param  rev_comp  0 only to reverse, otherwise the compliment will also be taken
  */
void
fmap_seq_reverse(fmap_seq_t *seq, int32_t rev_comp);

#endif
