#ifndef FMAP_SFF_IO_H_
#define FMAP_SFF_IO_H_

#include <stdint.h>
#include "../seq/fmap_sff.h"
#include "fmap_file.h"

/*! 
  A SFF Reading Library
  */

/*! 
  structure for reading SFFs
  */
typedef struct {
    fmap_file_t *fp;  /*!< pointer to the file structure */
    fmap_sff_header_t *gheader;  /*!< pointer to the global SFF header */
    int32_t n_read;  /*!< the number of SFF reads read */
} fmap_sff_io_t;

/*! 
  initializes sff reading structure
  @param  fp  a pointer to a file structure from which to read
  @return     pointer to the initialized memory for reading in sffs
  */
inline fmap_sff_io_t *
fmap_sff_io_init(fmap_file_t *fp);

/*! 
  destroys sff reading structure
  @param  sffio  a pointer to the sff structure
  */
inline void 
fmap_sff_io_destroy(fmap_sff_io_t *sffio);

/*! 
  reads in a sff structure
  @param  sffio  a pointer to a previously initialized sff structure
  @param  sff   the sff structure in which to store the data
  @return       1 if successful, -1 if unsuccessful (EOF)
  */
int 
fmap_sff_io_read(fmap_sff_io_t *sffio, fmap_sff_t *sff);

/*! 
  reads sffs into a buffer
  @param  sffio           a pointer to a previously initialized sff structure
  @param  sff_buffer     the sff structure in which to store the data
  @param  buffer_length  the number of sffs to read
  @return                the number of sffs read
  */
int
fmap_sff_io_read_buffer(fmap_sff_io_t *sffio, fmap_sff_t **sff_buffer, int32_t buffer_length);

#endif
