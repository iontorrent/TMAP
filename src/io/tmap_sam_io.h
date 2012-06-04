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
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the ID tag from the SAM Header
  */
#define tmap_sam_io_get_rg_id(samio) (samio->rg_tags[0]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the CN tag from the SAM Header
  */
#define tmap_sam_io_get_rg_cn(samio) (samio->rg_tags[1]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the DS tag from the SAM Header
  */
#define tmap_sam_io_get_rg_ds(samio) (samio->rg_tags[2]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the DT tag from the SAM Header
  */
#define tmap_sam_io_get_rg_dt(samio) (samio->rg_tags[3]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the FO tag from the SAM Header
  */
#define tmap_sam_io_get_rg_fo(samio) (samio->rg_tags[4]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the KS tag from the SAM Header
  */
#define tmap_sam_io_get_rg_ks(samio) (samio->rg_tags[5]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the LB tag from the SAM Header
  */
#define tmap_sam_io_get_rg_lb(samio) (samio->rg_tags[6]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the PG tag from the SAM Header
  */
#define tmap_sam_io_get_rg_pg(samio) (samio->rg_tags[7]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the PI tag from the SAM Header
  */
#define tmap_sam_io_get_rg_pi(samio) (samio->rg_tags[8]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the PL tag from the SAM Header
  */
#define tmap_sam_io_get_rg_pl(samio) (samio->rg_tags[9]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the PU tag from the SAM Header
  */
#define tmap_sam_io_get_rg_pu(samio) (samio->rg_tags[10]) 
/*!
  @param  samio  a pointer to a previously initialized SAM/BAM structure
  @return the SM tag from the SAM Header
  */
#define tmap_sam_io_get_rg_sm(samio) (samio->rg_tags[11]) 

/*! 
  initializes SAM reading structure
  @param  fp  a pointer to a file structure from which to read
  @return     pointer to the initialized memory for reading in SAMs/BAMs
  */
inline tmap_sam_io_t *
tmap_sam_io_init(const char *fn);

/*! 
  initializes BAM reading structure
  @param  fp  a pointer to a file structure from which to read
  @return     pointer to the initialized memory for reading in SAMs/BAMs
  */
inline tmap_sam_io_t *
tmap_bam_io_init(const char *fn);

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
