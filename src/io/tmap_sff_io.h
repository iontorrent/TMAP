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
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the ID tag for the SAM Header
  */
#define tmap_sff_io_get_rg_id(sffio) (NULL) 
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the CN tag for the SAM Header
  */
#define tmap_sff_io_get_rg_cn(sffio) (NULL) 
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the DS tag for the SAM Header
  */
#define tmap_sff_io_get_rg_ds(sffio) (NULL) 
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the DT tag for the SAM Header
  */
#define tmap_sff_io_get_rg_dt(sffio) (NULL) 
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the FO tag for the SAM Header
  */
#define tmap_sff_io_get_rg_fo(sffio) (sffio->gheader->flow->s)
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the KS tag for the SAM Header
  */
#define tmap_sff_io_get_rg_ks(sffio) (sffio->gheader->key->s)
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the LB tag for the SAM Header
  */
#define tmap_sff_io_get_rg_lb(sffio) (NULL) 
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the PG tag for the SAM Header
  */
#define tmap_sff_io_get_rg_pg(sffio) (NULL) 
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the PI tag for the SAM Header
  */
#define tmap_sff_io_get_rg_pi(sffio) (NULL) 
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the PL tag for the SAM Header
  */
#define tmap_sff_io_get_rg_pl(sffio) (NULL) 
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the PU tag for the SAM Header
  */
#define tmap_sff_io_get_rg_pu(sffio) (NULL) 
/*!
  @param  sffio  a pointer to a previously initialized SFF structure
  @return the SM tag for the SAM Header
  */
#define tmap_sff_io_get_rg_sm(sffio) (NULL) 

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

/*!
  @param  sffio  a pointer to a previously initialized sff structure
  @param  n     stores the number of rg ids 
  @return   the header structure (rg-ids x rg tags)
 */
char***
tmap_sff_io_get_rg_header(tmap_sff_io_t *sffio, int32_t *n);

#endif
