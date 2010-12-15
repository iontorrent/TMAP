/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_FQ_IO_H_
#define TMAP_FQ_IO_H_

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

#include "../seq/tmap_fq.h"
#include "tmap_file.h"

/*! 
  A FASTQ Reading Library
  */

#define TMAP_STREAM_BUFFER_SIZE 4096

/*! 
  */
typedef struct {
    char *buf;  /*!< the character buffer */
    int32_t begin;  /*!< the index of the next character in the buffer */
    int32_t end;  /*!< the number of characters last read */
    int32_t is_eof;  /*!< 1 if the EOF marker has been reached, 0 otherwise */
    tmap_file_t *f;  /*!< the file pointer associated with this stream */
    int32_t bufsize;  /*!< the size of the character buffer */
} tmap_stream_t; 

/*! 
  structure for reading FASTA/FASTQ strings
  */
typedef struct {
    int last_char;  /*!< the last character read */
    tmap_stream_t *f;  /*!< pointer to the file structure */
} tmap_fq_io_t;

/*! 
  initializes fastq reading structure
  @param  fp  a pointer to a file structure from which to read
  @return     pointer to the initialized memory for reading in fastqs
  */
inline tmap_fq_io_t *
tmap_fq_io_init(tmap_file_t *fp);

/*! 
  destroys fastq reading structure
  @param  fqio  a pointer to the fastq structure
  */
inline void 
tmap_fq_io_destroy(tmap_fq_io_t *fqio);

/*! 
  reads in a reading structure
  @param  fqio  a pointer to a previously initialized fastq structure
  @param  fq    the fastq structure in which to store the data
  @return       the length of the fastq read, -1 indicates an a EOF, -2 indicates a truncated quality string
  */
int 
tmap_fq_io_read(tmap_fq_io_t *fqio, tmap_fq_t *fq);

/*! 
  reads fastqs into a buffer
  @param  fqio           a pointer to a previously initialized fastq structure
  @param  fq_buffer      the fastq structure in which to store the data
  @param  buffer_length  the number of fastqs to read
  @return                the number of fastqs read
  */
int
tmap_fq_io_read_buffer(tmap_fq_io_t *fqio, tmap_fq_t **fq_buffer, int32_t buffer_length);

#endif
