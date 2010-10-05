#ifndef FMAP_FILE_H_
#define FMAP_FILE_H_

#include <stdio.h>
#include <zlib.h>
#ifndef DISABLE_BZ2
#include <bzlib.h>
#endif 
#include <config.h>
#include <stdarg.h>

/*! 
  File handling routines analgous to those in stdio.h
  */

/*! 
  @details  the various supported file compression types
  */
enum {
    FMAP_FILE_NO_COMPRESSION=0,  /*!< no compression */
    FMAP_FILE_BZ2_COMPRESSION,  /*!< bzip2 compression */
    FMAP_FILE_GZ_COMPRESSION  /*!< gzip compression */
};

/*! 
  @details  the type of bzip2 stream (read/write)
  */
enum {
    FMAP_FILE_BZ2_READ=0,  /*!< a reading bzip2 stream  */
    FMAP_FILE_BZ2_WRITE  /*!< a writing bzip2 stream */
};

/*! 
  structure aggregating common file structures
  */
typedef struct {
    FILE *fp;  /*!< stdio file pointer */
    gzFile gz;  /*!< gzip file pointer */
    BZFILE *bz2;  /*!< bz2 file pointer */
    int32_t c;  /*!< compression type */
    char unused[BZ_MAX_UNUSED];  /*!< for bz2 function 'BZ2_bzReadOpen' */
    int32_t n_unused;  /*!< for bz2 function 'BZ2_bzReadGetUnused' */
    int32_t bzerror;  /*!< stores the last BZ2 error */
    int32_t open_type;  /*!< the type of bzip2 stream */
#ifndef DISABLE_BZ2
#endif
#ifndef DISABLE_BZ2
#endif
} fmap_file_t;

extern fmap_file_t *fmap_file_stdout; // to use, initialize this in your main
extern fmap_file_t *fmap_file_stderr; // to use, initialize this in your main

/*! 
  emulates fopen from stdio.h
  @param  path         filename to open
  @param  mode         access modes
  @param  compression  compression type
  @return              a pointer to the initialized file structure
  */
fmap_file_t *
fmap_file_fopen(const char* path, const char *mode, int32_t compression);

/*! 
  emulates fdopen from stdio.h
  @param  filedes       file descriptor to open
  @param  mode         access modes
  @param  compression  compression type
  @return              a pointer to the initialized file structure
  */
fmap_file_t *
fmap_file_fdopen(int filedes, const char *mode, int32_t compression);

/*! 
  closes the file associated with the file pointer
  @param  fp  pointer to the file structure to close
  @return     a pointer to the initialized file structure
  */
void 
fmap_file_fclose(fmap_file_t *fp); 

/*! 
  emulates fread from stdio.h
  @param  ptr    pointer to a block of memory with a minimum size of (size*count) bytes
  @param  size   size in bytes of each element to be read
  @param  count  number of elements, each one with a size of size bytes
  @param  fp     pointer to the file structure from which to read
  @return        the number of elements read (should equal count)
  */
size_t 
fmap_file_fread(void *ptr, size_t size, size_t count, fmap_file_t *fp);

/*! 
  emulates gzread from zlib.h
  @param  fp   pointer to the file structure from which to read
  @param  ptr  pointer to a block of memory with a minimum size of len bytes
  @param  len  number bytes to read
  @return      the number of elements read (should equal count)
  */
int 
fmap_file_fread2(fmap_file_t *fp, void *ptr, unsigned int len);

/*! 
  emulates fgetc from stdio.h
  @param  fp  pointer to the file structure from which to read
  @return     the character read is returned as an int value 
  */
int 
fmap_file_fgetc(fmap_file_t *fp);

/*! 
  emulates fwrite from stdio.h
  @param  ptr    pointer to the array of elements to be written
  @param  size   size in bytes of each element to be written
  @param  count  number of elements, each one with a size of size bytes
  @param  fp     pointer to the file structure to which to write
  @return        the number of elements written (should equal count)
  */
size_t 
fmap_file_fwrite(void *ptr, size_t size, size_t count, fmap_file_t *fp); 

/*! 
  emulates vfprintf from stdio.h
  @param  fp      pointer to the file structure to which to write
  @param  format  the text and format of what to print
  @param  ap      depending on the format string, the function may expect a sequence of additional arguments
  @return         the number of characters written
  */
int32_t
fmap_file_vfprintf(fmap_file_t *fp, const char *format, va_list ap);

/*! 
  emulates fprintf from stdio.h
  @param  fp      pointer to the file structure to which to write
  @param  format  the text and format of what to print
  @param  ...     depending on the format string, the function may expect a sequence of additional arguments
  @return         the number of characters written
  */
int32_t
fmap_file_fprintf(fmap_file_t *fp, const char *format, ...);

/*!
  emulates fflush from stdio.h
  @param  fp       pointer to the file structure to which to flush (the file should have been opened for writing)
  @param  gz_flush  if the file is writing gzip compressed output, 0 will flush the compressed buffer and 1 will flush all data
  @details         this will have no effect if the file is writing bzip2 data. 
  */
int32_t
fmap_file_fflush(fmap_file_t *fp, int32_t gz_flush);

#endif
