#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <bzlib.h>
#include <zlib.h>
#include <ctype.h>
#include <config.h>

#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_io.h"

fmap_file_t *
fmap_file_fopen(const char* path, const char *mode, int32_t compression) 
{
  fmap_file_t *fp = NULL;
  int32_t open_ok = 1;

  if(NULL == strchr(mode, 'r') && NULL == strchr(mode, 'w')) {
      fmap_error("mode", Exit, OutOfRange);
  }

  fp = fmap_calloc(1, sizeof(fmap_file_t), "fp");
  fp->fp=NULL;
#ifndef DISABLE_BZ2 
  fp->bz2=NULL;
  fp->n_unused=0;
#endif
  fp->gz=NULL;
  fp->c=compression;

  switch(fp->c) {
    case FMAP_FILE_NO_COMPRESSION:
      fp->fp = fopen(path, mode);
      if(NULL == fp->fp) {
          free(fp); 
          open_ok = 0;
          break;
      }
      break;
#ifndef DISABLE_BZ2 
    case FMAP_FILE_BZ2_COMPRESSION:
      fp->fp = fopen(path, mode);
      if(NULL == fp->fp) {
          free(fp); 
          open_ok = 0;
          break;
      }
      if(NULL != strchr(mode, 'r')) {
          fp->open_type = FMAP_FILE_BZ2_READ;
          fp->bz2 = BZ2_bzReadOpen(&fp->bzerror, fp->fp, 0, 0, fp->unused, fp->n_unused);
          if(NULL == fp->bz2) {
              free(fp); 
              open_ok = 0;
              break;
          }
      }
      else {
          fp->open_type = FMAP_FILE_BZ2_WRITE;
          // 900K blockSize100k
          // 30 workFactor
          fp->bz2 = BZ2_bzWriteOpen(&fp->bzerror, fp->fp, 9, 0, 30); 
      }
      break;
#endif
    case FMAP_FILE_GZ_COMPRESSION:
      fp->gz = gzopen(path, mode);
      if(NULL == fp->gz) {
          free(fp); 
          open_ok = 0;
          break;
      }
      break;
    default:
      fmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(0 == open_ok) {
      fmap_error(path, Exit, OpenFileError);
  }

  return fp;
}

fmap_file_t *
fmap_file_fdopen(int filedes, const char *mode, int32_t compression) 
{
  fmap_file_t *fp = NULL;
  int32_t open_ok = 1;

  if(NULL == strchr(mode, 'r') && NULL == strchr(mode, 'w')) {
      fmap_error("mode", Exit, OutOfRange);
  }

  fp = fmap_calloc(1, sizeof(fmap_file_t), "fp");
  fp->fp=NULL;
#ifndef DISABLE_BZ2 
  fp->bz2=NULL;
  fp->n_unused=0;
#endif
  fp->gz=NULL;
  fp->c=compression;

  switch(fp->c) {
    case FMAP_FILE_NO_COMPRESSION:
      fp->fp = fdopen(filedes, mode);
      if(NULL == fp->fp) {
          free(fp); 
          open_ok = 0;
      }
      break;
#ifndef DISABLE_BZ2 
    case FMAP_FILE_BZ2_COMPRESSION:
      fp->fp = fdopen(filedes, mode);
      if(NULL == fp->fp) {
          free(fp); 
          open_ok = 0;
      }
      if(NULL != strchr(mode, 'r')) {
          fp->open_type = FMAP_FILE_BZ2_READ;
          fp->bz2 = BZ2_bzReadOpen(&fp->bzerror, fp->fp, 0, 0, fp->unused, fp->n_unused);
          if(NULL == fp->bz2) {
              free(fp); 
              open_ok = 0;
              break;
          }
      }
      else {
          fp->open_type = FMAP_FILE_BZ2_WRITE;
          // 900K blockSize100k
          // 30 workFactor
          fp->bz2 = BZ2_bzWriteOpen(&fp->bzerror, fp->fp, 9, 0, 30); 
      }
#endif
    case FMAP_FILE_GZ_COMPRESSION:
      fp->gz = gzdopen(filedes, mode);
      if(NULL == fp->gz) {
          free(fp); 
          open_ok = 0;
          break;
      }
      break;
    default:
      fmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(0 == open_ok) {
      fmap_error(NULL, Exit, OpenFileError);
  }

  return fp;
}

void 
fmap_file_fclose(fmap_file_t *fp) 
{
  int closed_ok = 1;
  switch(fp->c) {
    case FMAP_FILE_NO_COMPRESSION:
#ifndef DISABLE_BZ2 
    case FMAP_FILE_BZ2_COMPRESSION:
      if(FMAP_FILE_BZ2_WRITE == fp->open_type) {
          BZ2_bzWriteClose(&fp->bzerror, fp->bz2, 0, NULL, NULL);
          if(fp->bzerror == BZ_IO_ERROR) {
              closed_ok = 0;
              break;
          }
      }
      if(EOF == fclose(fp->fp) ) {
          closed_ok = 0;
          break;
      }
      break;
#endif
    case FMAP_FILE_GZ_COMPRESSION:
      if(EOF == gzclose(fp->gz)) {
          closed_ok = 0;
          break;
      }
      break;
    default:
      fmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(closed_ok == 0) {
      fmap_error(NULL, Exit, OpenFileError);
  }

  free(fp);
}

size_t 
fmap_file_fread(void *ptr, size_t size, size_t count, fmap_file_t *fp) 
{
  size_t num_read = 0;
  int error;
#ifndef DISABLE_BZ2 
  int32_t nbuf=0, i;
  void *unused_tmp_void=NULL;
  char *unused_tmp=NULL;
#endif

  switch(fp->c) {
    case FMAP_FILE_NO_COMPRESSION:
      num_read = fread(ptr, size, count, fp->fp);
      break;
#ifndef DISABLE_BZ2 
    case FMAP_FILE_BZ2_COMPRESSION:
      while(0 == nbuf && 
            !(BZ_STREAM_END == fp->bzerror && 0 == fp->n_unused && feof(fp->fp))) {
          nbuf = BZ2_bzRead(&fp->bzerror, fp->bz2, ptr, size*count);
          if(BZ_OK == fp->bzerror) {
              break;
          }
          else if(BZ_STREAM_END == fp->bzerror) {
              // Get unused
              BZ2_bzReadGetUnused(&fp->bzerror, fp->bz2, &unused_tmp_void, &fp->n_unused);
              if(BZ_OK != fp->bzerror) fmap_error("BZ2_bzReadGetUnused", Exit, OutOfRange); 
              unused_tmp = (char*)unused_tmp_void;
              for(i=0;i<fp->n_unused;i++) {
                  fp->unused[i] = unused_tmp[i];
              }
              // Close
              BZ2_bzReadClose(&fp->bzerror, fp->bz2);
              if(BZ_OK != fp->bzerror) fmap_error("BZ2_bzReadClose", Exit, OutOfRange); 
              fp->bzerror = BZ_STREAM_END; // set to the stream end for next call to this function
              // Open again if necessary
              if(0 == fp->n_unused && feof(fp->fp)) {
                  break;
              }
              else {
                  fp->bz2 = BZ2_bzReadOpen(&fp->bzerror, fp->fp, 0, 0, fp->unused, fp->n_unused);
                  if(NULL == fp->bz2) fmap_error("fp->bz2", Exit, OpenFileError);
              }
          }
          else {
              fmap_error("fp->c", Exit, ReadFileError);
          }
      }
      num_read = nbuf / size;
      break;
#endif
    case FMAP_FILE_GZ_COMPRESSION:
      num_read = gzread(fp->gz, ptr, size*count) / size;
      break;
    default:
      fmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(num_read != count) {
      // check if it was an end of file/stream
      switch(fp->c) {
        case FMAP_FILE_NO_COMPRESSION:
          if(0 != feof(fp->fp)) {
              return num_read;
          }
          break;
#ifndef DISABLE_BZ2 
        case FMAP_FILE_BZ2_COMPRESSION:
          if(BZ_STREAM_END == fp->bzerror) {
              return num_read;
          }
          break;
#endif
        case FMAP_FILE_GZ_COMPRESSION:
          gzerror(fp->gz, &error);
          if(Z_STREAM_END == error) {
              return num_read;
          }
          break;
        default:
          fmap_error("fp->c", Exit, OutOfRange);
          break;
      }
      fmap_error(NULL, Exit, ReadFileError);
  }

  return num_read;
}

int 
fmap_file_fread2(fmap_file_t *fp, void *ptr, unsigned int len)
{
  return fmap_file_fread(ptr, sizeof(char), (size_t)len, fp);
}

int 
fmap_file_fgetc(fmap_file_t *fp)
{
  char c;
  fmap_file_fread(&c, sizeof(char), 1, fp);
  return (int)c;
}

size_t 
fmap_file_fwrite(void *ptr, size_t size, size_t count, fmap_file_t *fp) 
{
  size_t num_written = 0;

  switch(fp->c) {
    case FMAP_FILE_NO_COMPRESSION:
      num_written = fwrite(ptr, size, count, fp->fp);
      break;
#ifndef DISABLE_BZ2 
    case FMAP_FILE_BZ2_COMPRESSION:
      num_written = BZ2_bzwrite(fp->bz2, ptr, size*count) / size;
      break;
#endif
    case FMAP_FILE_GZ_COMPRESSION:
      num_written = gzwrite(fp->gz, ptr, size*count) / size;
      break;
    default:
      fmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(num_written != count) {
      fmap_error(NULL, Exit, WriteFileError);
  }

  return num_written;
}
