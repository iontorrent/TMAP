#ifndef TMAP_SFF_IO_H_
#define TMAP_SFF_IO_H_

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
    int32_t n_read;  /*!< the number of SFF reads read */
} tmap_sff_io_t;

/*! 
  initializes sff reading structure
  @param  fp  a pointer to a file structure from which to read
  @return     pointer to the initialized memory for reading in sffs
  */
inline tmap_sff_io_t *
tmap_sff_io_init(tmap_file_t *fp);

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
