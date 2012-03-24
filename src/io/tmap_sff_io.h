/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SFF_IO_H
#define TMAP_SFF_IO_H

#include <stdint.h>
#include "../seq/tmap_sff.h"
#include "tmap_file.h"

/*! 
  A SFF Reading Library
  */

/*! 
  structure for reading SFFs
  */
typedef struct {
    tmap_file_t *fp;  /*!< pointer to the file structure */
    tmap_sff_header_t *gheader;  /*!< pointer to the global SFF header */
    uint32_t n_read;  /*!< the number of SFF reads read */
    int32_t early_eof_ok;  /*!< 0 if the number of reads to read in should match the SFF header, 0 otherwise */
} tmap_sff_io_t;

/*!
  @param  sffio  a pointer to a previously initialized SAM/BAM structure
  @return the FO tag from the SAM Header
  */
#define tmap_sff_io_get_rg_fo(sffio) (sffio->gheader->flow->s)
/*!
  @param  sffio  a pointer to a previously initialized SAM/BAM structure
  @return the KS tag from the SAM Header
  */
#define tmap_sff_io_get_rg_ks(sffio) (sffio->gheader->key->s)

/*! 
  initializes sff reading structure
  @param  fp  a pointer to a file structure from which to read
  @return     pointer to the initialized memory for reading in sffs
  */
inline tmap_sff_io_t *
tmap_sff_io_init(tmap_file_t *fp);

/*! 
  initializes sff reading structure
  @param  fp  a pointer to a file structure from which to read
  @param  early_eof_ok  set this to 1 if the number of reads to be read in disagrees with the SFF header
  @return     pointer to the initialized memory for reading in sffs
  */
inline tmap_sff_io_t *
tmap_sff_io_init2(tmap_file_t *fp, int32_t early_eof_ok);

/*! 
  destroys sff reading structure
  @param  sffio  a pointer to the sff structure
  */
inline void 
tmap_sff_io_destroy(tmap_sff_io_t *sffio);

/*! 
  reads in a sff structure
  @param  sffio  a pointer to a previously initialized sff structure
  @param  sff   the sff structure in which to store the data
  @return       1 if successful, -1 if unsuccessful (EOF)
  */
int 
tmap_sff_io_read(tmap_sff_io_t *sffio, tmap_sff_t *sff);

/*! 
  reads sffs into a buffer
  @param  sffio           a pointer to a previously initialized sff structure
  @param  sff_buffer     the sff structure in which to store the data
  @param  buffer_length  the number of sffs to read
  @return                the number of sffs read
  */
int
tmap_sff_io_read_buffer(tmap_sff_io_t *sffio, tmap_sff_t **sff_buffer, int32_t buffer_length);

#endif
