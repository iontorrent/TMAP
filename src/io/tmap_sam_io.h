/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SAM_IO_H
#define TMAP_SAM_IO_H

#include <config.h>

#include "../samtools/bam.h"
#include "../samtools/sam.h"

/*! 
  A SAM/BAM Reading Library
  */

/*! 
*/
typedef struct _tmap_sam_io_t {
    samfile_t *fp;  /*!< the file pointer to the SAM/BAM file */
} tmap_sam_io_t;

#include "../seq/tmap_sam.h"

/*! 
  initializes SAM reading structure
  @param  fn  the input file name, or "-" for stdin
  @return     pointer to the initialized memory for reading in SAMs/BAMs
  */
inline tmap_sam_io_t *
tmap_sam_io_init(const char *fn);

/*! 
  initializes BAM reading structure
  @param  fn  the input file name, or "-" for stdin
  @return     pointer to the initialized memory for reading in SAMs/BAMs
  */
inline tmap_sam_io_t *
tmap_bam_io_init(const char *fn);

/*!
  initialized SAM/BAM writing structure
  @param  fn  the output file name, or "-" for stdout
  @param  mode  the mode; must be one of "wh", "wb", or "wbu"
  @param  header  the output BAM Header
  @return  a pointer to the initialized memory for writing SAMs/BAMs
  */
inline tmap_sam_io_t *
tmap_sam_io_init2(const char *fn, const char *mode,
                  bam_header_t *header);

/*! 
  destroys SAM/BAM reading structure
  @param  samio  a pointer to the SAM/BAM structure
  */
void
tmap_sam_io_destroy(tmap_sam_io_t *samio);

/*! 
  reads in a reading structure
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @param  sam    the SAM/BAM structure in which to store the data
  @return        1 if the SAM/BAM record was read in correctly, otherwise -1 indicates an a EOF
  */
int32_t
tmap_sam_io_read(tmap_sam_io_t *samio, tmap_sam_t *sam);

/*! 
  reads SAMs/BAMs into a buffer
  @param  samio           a pointer to a previously initialized SAM/BAM structure
  @param  sam_buffer      the SAM/BAM structure in which to store the data
  @param  buffer_length   the number of SAMs/BAMs to read
  @return                 the number of SAMs/BAMs read
  */
int32_t
tmap_sam_io_read_buffer(tmap_sam_io_t *samio, tmap_sam_t **sam_buffer, int32_t buffer_length);

/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return   the SAM header structure 
 */
sam_header_t*
tmap_sam_io_get_sam_header(tmap_sam_io_t *samio);

#endif
