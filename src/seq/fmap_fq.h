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

/* Nils Homer - modified not to be define-ized */

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "../util/fmap_string.h"
#include "../util/fmap_definitions.h"

/*! 
  A Library for FASTQ data
  */

/*! 
  structure for holding FASTA/FASTQ strings
  */
typedef struct {
    fmap_string_t *name;  /*!< the name string */
    fmap_string_t *comment;  /*!< the comment string */
    fmap_string_t *seq;  /*!< the sequence string */
    fmap_string_t *qual;  /*!< the quality string */
    int32_t is_int;  /*!< 1 if the sequence is in integer format, 0 otherwise  */
} fmap_fq_t;

/*! 
  initializes sequence read structure
  @return     pointer to the initialized memory 
  */
inline fmap_fq_t *
fmap_fq_init();

/*! 
  @param  fq  a pointer to the sequence structure
  */
inline void 
fmap_fq_destroy(fmap_fq_t *fq);

/*! 
  clones the given sequence read structure
  @param  fq  pointer to the sequence read structure to be copied  
  @return     pointer to the initialized memory 
  */
inline fmap_fq_t *
fmap_fq_clone(fmap_fq_t *fq);

/*! 
  reverse compliments the sequence and reverse the qualities
  @param  fq  a pointer to a sequence structure
  */
void
fmap_fq_reverse_compliment(fmap_fq_t *fq);

/*! 
  converts bases to integer values
  @param  fq  a pointer to a sequence structure
  */
void
fmap_fq_to_int(fmap_fq_t *fq);

/*! 
  converts bases to character values
  @param  fq  a pointer to a sequence structure
  */
void
fmap_fq_to_char(fmap_fq_t *fq);

/*!
  gets the read's bases
  @param  fq  a pointer to a sequence structure
 */
inline fmap_string_t *
fmap_fq_get_bases(fmap_fq_t *fq);

/*!
  gets the read's qualities
  @param  fq  a pointer to a sequence structure
 */
inline fmap_string_t *
fmap_fq_get_qualities(fmap_fq_t *fq);

#endif
