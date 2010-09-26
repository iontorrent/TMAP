#ifndef FMAP_FILE_H_
#define FMAP_FILE_H_

#include <stdio.h>
#include <zlib.h>
#ifndef DISABLE_BZ2
#include <bzlib.h>
#endif 
#include <config.h>
#include <stdarg.h>

/*! @header
  @abstract  File handling routines analgous to those in stdio.h
  */

/*! 
  @discussion  the various supported file compression types
  */
enum {
    FMAP_FILE_NO_COMPRESSION=0,  /*!< no compression */
    FMAP_FILE_BZ2_COMPRESSION,  /*!< bzip2 compression */
    FMAP_FILE_GZ_COMPRESSION  /*!< gzip compression */
};

/*! 
  @discussion  the type of bzip2 stream (read/write)
  */
enum {
    FMAP_FILE_BZ2_READ=0,  /*!< a reading bzip2 stream  */
    FMAP_FILE_BZ2_WRITE  /*!< a writing bzip2 stream */
};

/*! @typedef
  @abstract         structure aggregating common file structures
  @field  fp         stdio file pointer
  @field  gz         gzip file pointer
  @field  c          compression type
  @field  bz2        bz2 file pointer
  @field  unused     for bz2 function 'BZ2_bzReadOpen'
  @field  n_unused   for bz2 function 'BZ2_bzReadGetUnused'
  @field  bzerror    stores the last BZ2 error
  @field  open_type  the type of bzip2 stream
  */
typedef struct {
    FILE *fp;
    gzFile gz;
#ifndef DISABLE_BZ2
    BZFILE *bz2;
#endif
    int32_t c;

#ifndef DISABLE_BZ2
    char unused[BZ_MAX_UNUSED];
    int32_t n_unused;
    int32_t bzerror;
    int32_t open_type;
#endif
} fmap_file_t;

extern fmap_file_t *fmap_file_stdout; // to use, initialize this in your main
extern fmap_file_t *fmap_file_stderr; // to use, initialize this in your main

/*! @function
  @abstract            emulates fopen from stdio.h
  @param  path         filename to open
  @param  mode         access modes
  @param  compression  compression type
  @return              a pointer to the initialized file structure
  */
fmap_file_t *
fmap_file_fopen(const char* path, const char *mode, int32_t compression);

/*! @function
  @abstract            emulates fdopen from stdio.h
  @param  filedes       file descriptor to open
  @param  mode         access modes
  @param  compression  compression type
  @return              a pointer to the initialized file structure
  */
fmap_file_t *
fmap_file_fdopen(int filedes, const char *mode, int32_t compression);

/*! @function
  @abstract   closes the file associated with the file pointer
  @param  fp  pointer to the file structure to close
  @return     a pointer to the initialized file structure
  */
void 
fmap_file_fclose(fmap_file_t *fp); 

/*! @function
  @abstract      emulates fread from stdio.h
  @param  ptr    pointer to a block of memory with a minimum size of (size*count) bytes
  @param  size   size in bytes of each element to be read
  @param  count  number of elements, each one with a size of size bytes
  @param  fp     pointer to the file structure from which to read
  @return        the number of elements read (should equal count)
  */
size_t 
fmap_file_fread(void *ptr, size_t size, size_t count, fmap_file_t *fp);

/*! @function
  @abstract    emulates gzread from zlib.h
  @param  fp   pointer to the file structure from which to read
  @param  ptr  pointer to a block of memory with a minimum size of len bytes
  @param  len  number bytes to read
  @return      the number of elements read (should equal count)
  */
int 
fmap_file_fread2(fmap_file_t *fp, void *ptr, unsigned int len);

/*! @function
  @abstract   emulates fgetc from stdio.h
  @param  fp  pointer to the file structure from which to read
  @return     the character read is returned as an int value 
  */
int 
fmap_file_fgetc(fmap_file_t *fp);

/*! @function
  @abstract      emulates fwrite from stdio.h
  @param  ptr    pointer to the array of elements to be written
  @param  size   size in bytes of each element to be written
  @param  count  number of elements, each one with a size of size bytes
  @param  fp     pointer to the file structure to which to write
  @return        the number of elements written (should equal count)
  */
size_t 
fmap_file_fwrite(void *ptr, size_t size, size_t count, fmap_file_t *fp); 

/*! @function
  @abstract       emulates vfprintf from stdio.h
  @param  fp      pointer to the file structure to which to write
  @param  format  the text and format of what to print
  @param  ap      depending on the format string, the function may expect a sequence of additional arguments
  @return         the number of characters written
  */
int32_t
fmap_file_vfprintf(fmap_file_t *fp, const char *format, va_list ap);

/*! @function
  @abstract       emulates fprintf from stdio.h
  @param  fp      pointer to the file structure to which to write
  @param  format  the text and format of what to print
  @param  ...     depending on the format string, the function may expect a sequence of additional arguments
  @return         the number of characters written
  */
int32_t
fmap_file_fprintf(fmap_file_t *fp, const char *format, ...);
#endif
