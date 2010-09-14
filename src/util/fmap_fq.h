#ifndef FMAP_FQ_H_
#define FMAP_FQ_H_

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
#include "fmap_definitions.h"

/*! @typedef
  @abstract         structure for holding FASTA/FASTQ strings
  @field  name       the name string
  @field  comment    the comment string
  @field  seq        the sequence string
  @field  qual       the quality string
  @field  is_int     1 if the sequence is in integer format, 0 otherwise 
  */
typedef struct {
    fmap_string_t name, comment, seq, qual;
    int32_t is_int;
} fmap_fq_t;

/*! @function
  @abstract   initializes sequence read structure
  @return     pointer to the initialized memory 
  */
inline fmap_fq_t *
fmap_fq_init();

/*! @function
  @abstract   
  @param  seq  a pointer to the sequence structure
  */
inline void 
fmap_fq_destroy(fmap_fq_t *seq);

/*! @function
  @abstract   clones the given sequence read structure
  @param  pointer to the sequence read structure to be copied  
  @return  pointer to the initialized memory 
  */
inline fmap_fq_t *
fmap_fq_clone(fmap_fq_t *seq);

/*! @function
  @abstract         reverse the seq and qual fields
  @param  seq       a pointer to a sequence structure
  @param  rev_comp  0 only to reverse, otherwise the compliment will also be taken
  */
void
fmap_fq_reverse(fmap_fq_t *seq, int32_t rev_comp);

/*! @function
  @abstract         converts bases to integer values
  @param  seq       a pointer to a sequence structure
  */
void
fmap_fq_to_int(fmap_fq_t *seq);

/*! @function
  @abstract         converts bases to character values
  @param  seq       a pointer to a sequence structure
  */
void
fmap_fq_to_char(fmap_fq_t *seq);

#endif
