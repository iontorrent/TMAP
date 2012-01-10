/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <bzlib.h>
#include <zlib.h>
#include <ctype.h>
#include <config.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "tmap_file.h"

tmap_file_t *
tmap_file_fopen(const char* path, const char *mode, int32_t compression) 
{
  tmap_file_t *fp = NULL;
  int32_t open_ok = 1;

  if(NULL == strchr(mode, 'r') && NULL == strchr(mode, 'w')) {
      tmap_error("mode", Exit, OutOfRange);
  }

  fp = tmap_calloc(1, sizeof(tmap_file_t), "fp");
  fp->fp=NULL;
#ifndef DISABLE_BZ2 
  fp->bz2=NULL;
  fp->n_unused=0;
#endif
  fp->gz=NULL;
  fp->c=compression;

  switch(fp->c) {
    case TMAP_FILE_NO_COMPRESSION:
      fp->fp = fopen(path, mode);
      if(NULL == fp->fp) {
          free(fp); 
          open_ok = 0;
          break;
      }
      break;
#ifndef DISABLE_BZ2 
    case TMAP_FILE_BZ2_COMPRESSION:
      fp->fp = fopen(path, mode);
      if(NULL == fp->fp) {
          free(fp); 
          open_ok = 0;
          break;
      }
      if(NULL != strchr(mode, 'r')) {
          fp->open_type = TMAP_FILE_BZ2_READ;
          fp->bz2 = BZ2_bzReadOpen(&fp->bzerror, fp->fp, 0, 0, fp->unused, fp->n_unused);
          if(NULL == fp->bz2) {
              free(fp); 
              open_ok = 0;
              break;
          }
      }
      else {
          fp->open_type = TMAP_FILE_BZ2_WRITE;
          // 900K blockSize100k
          // 30 workFactor
          fp->bz2 = BZ2_bzWriteOpen(&fp->bzerror, fp->fp, 9, 0, 30); 
      }
      break;
#endif
    case TMAP_FILE_GZ_COMPRESSION:
      fp->gz = gzopen(path, mode);
      if(NULL == fp->gz) {
          free(fp); 
          open_ok = 0;
          break;
      }
      break;
    default:
      tmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(0 == open_ok) {
      tmap_error(path, Exit, OpenFileError);
  }

  return fp;
}

tmap_file_t *
tmap_file_fdopen(int filedes, const char *mode, int32_t compression) 
{
  tmap_file_t *fp = NULL;
  int32_t open_ok = 1;

  if(NULL == strchr(mode, 'r') && NULL == strchr(mode, 'w')) {
      tmap_error("mode", Exit, OutOfRange);
  }

  fp = tmap_calloc(1, sizeof(tmap_file_t), "fp");
  fp->fp=NULL;
#ifndef DISABLE_BZ2 
  fp->bz2=NULL;
  fp->n_unused=0;
#endif
  fp->gz=NULL;
  fp->c=compression;

  switch(fp->c) {
    case TMAP_FILE_NO_COMPRESSION:
      fp->fp = fdopen(filedes, mode);
      if(NULL == fp->fp) {
          free(fp); 
          open_ok = 0;
      }
      break;
#ifndef DISABLE_BZ2 
    case TMAP_FILE_BZ2_COMPRESSION:
      fp->fp = fdopen(filedes, mode);
      if(NULL == fp->fp) {
          free(fp); 
          open_ok = 0;
      }
      if(NULL != strchr(mode, 'r')) {
          fp->open_type = TMAP_FILE_BZ2_READ;
          fp->bz2 = BZ2_bzReadOpen(&fp->bzerror, fp->fp, 0, 0, fp->unused, fp->n_unused);
          if(NULL == fp->bz2) {
              free(fp); 
              open_ok = 0;
              break;
          }
      }
      else {
          fp->open_type = TMAP_FILE_BZ2_WRITE;
          // 900K blockSize100k
          // 30 workFactor
          fp->bz2 = BZ2_bzWriteOpen(&fp->bzerror, fp->fp, 9, 0, 30); 
      }
#endif
    case TMAP_FILE_GZ_COMPRESSION:
      fp->gz = gzdopen(filedes, mode);
      if(NULL == fp->gz) {
          free(fp); 
          open_ok = 0;
          break;
      }
      break;
    default:
      tmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(0 == open_ok) {
      tmap_error(NULL, Exit, OpenFileError);
  }

  return fp;
}

void 
tmap_file_fclose1(tmap_file_t *fp, int32_t close_underlyingfp) 
{
  int closed_ok = 1;
  switch(fp->c) {
    case TMAP_FILE_NO_COMPRESSION:
      if(1 == close_underlyingfp) {
          if(EOF == fclose(fp->fp) ) {
              closed_ok = 0;
              break;
          }
      }
      break;
#ifndef DISABLE_BZ2 
    case TMAP_FILE_BZ2_COMPRESSION:
      if(TMAP_FILE_BZ2_WRITE == fp->open_type) {
          BZ2_bzWriteClose(&fp->bzerror, fp->bz2, 0, NULL, NULL);
          if(fp->bzerror == BZ_IO_ERROR) {
              closed_ok = 0;
              break;
          }
      }
      if(1 == close_underlyingfp) {
          if(EOF == fclose(fp->fp) ) {
              closed_ok = 0;
              break;
          }
      }
      break;
#endif
    case TMAP_FILE_GZ_COMPRESSION:
      if(EOF == gzclose(fp->gz)) {
          closed_ok = 0;
          break;
      }
      break;
    default:
      tmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(closed_ok == 0) {
      tmap_error(NULL, Exit, OpenFileError);
  }

  free(fp);
}

void 
tmap_file_fclose(tmap_file_t *fp) 
{
  tmap_file_fclose1(fp, 1);
}

size_t 
tmap_file_fread(void *ptr, size_t size, size_t count, tmap_file_t *fp) 
{
  size_t num_read = 0, to_read, cur_read;
  int error;
#ifndef DISABLE_BZ2 
  int32_t nbuf=0, i;
  void *unused_tmp_void=NULL;
  char *unused_tmp=NULL;
#endif

  switch(fp->c) {
    case TMAP_FILE_NO_COMPRESSION:
      //cur_read = fread(ptr, size, count, fp->fp);
      // NB: some fread's are buggy when reading in large files, 
      // like apple's fread, so use a buffered approach
      while(num_read < count) {
          to_read = count - num_read;
          if(TMAP_1GB < to_read * size) {
              to_read = (size_t)(TMAP_1GB / size);
          }
          cur_read = fread(ptr + (num_read * size), size, to_read, fp->fp);
          num_read += cur_read;
          if(cur_read != to_read) break; // error 
      }
      break;
#ifndef DISABLE_BZ2 
    case TMAP_FILE_BZ2_COMPRESSION:
      while(0 == nbuf && 
            !(BZ_STREAM_END == fp->bzerror && 0 == fp->n_unused && feof(fp->fp))) {
          nbuf = BZ2_bzRead(&fp->bzerror, fp->bz2, ptr, size*count);
          if(BZ_OK == fp->bzerror) {
              break;
          }
          else if(BZ_STREAM_END == fp->bzerror) {
              // Get unused
              BZ2_bzReadGetUnused(&fp->bzerror, fp->bz2, &unused_tmp_void, &fp->n_unused);
              if(BZ_OK != fp->bzerror) tmap_error("BZ2_bzReadGetUnused", Exit, OutOfRange); 
              unused_tmp = (char*)unused_tmp_void;
              for(i=0;i<fp->n_unused;i++) {
                  fp->unused[i] = unused_tmp[i];
              }
              // Close
              BZ2_bzReadClose(&fp->bzerror, fp->bz2);
              if(BZ_OK != fp->bzerror) tmap_error("BZ2_bzReadClose", Exit, OutOfRange); 
              fp->bzerror = BZ_STREAM_END; // set to the stream end for next call to this function
              // Open again if necessary
              if(0 == fp->n_unused && feof(fp->fp)) {
                  break;
              }
              else {
                  fp->bz2 = BZ2_bzReadOpen(&fp->bzerror, fp->fp, 0, 0, fp->unused, fp->n_unused);
                  if(NULL == fp->bz2) tmap_error("fp->bz2", Exit, OpenFileError);
              }
          }
          else {
              tmap_error("fp->c", Exit, ReadFileError);
          }
      }
      num_read = nbuf / size;
      break;
#endif
    case TMAP_FILE_GZ_COMPRESSION:
      num_read = gzread(fp->gz, ptr, size*count) / size;
      break;
    default:
      tmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(num_read != count) {
      // check if it was an end of file/stream
      switch(fp->c) {
        case TMAP_FILE_NO_COMPRESSION:
          if(0 != feof(fp->fp)) {
              return num_read;
          }
          break;
#ifndef DISABLE_BZ2 
        case TMAP_FILE_BZ2_COMPRESSION:
          if(BZ_STREAM_END == fp->bzerror) {
              return num_read;
          }
          break;
#endif
        case TMAP_FILE_GZ_COMPRESSION:
          gzerror(fp->gz, &error);
          if(Z_STREAM_END == error) {
              return num_read;
          }
          break;
        default:
          tmap_error("fp->c", Exit, OutOfRange);
          break;
      }
      tmap_error(NULL, Exit, ReadFileError);
  }

  return num_read;
}

int 
tmap_file_fread2(tmap_file_t *fp, void *ptr, unsigned int len)
{
  return tmap_file_fread(ptr, sizeof(char), (size_t)len, fp);
}

int 
tmap_file_fgetc(tmap_file_t *fp)
{
  char c;
  tmap_file_fread(&c, sizeof(char), 1, fp);
  return (int)c;
}

size_t 
tmap_file_fwrite(void *ptr, size_t size, size_t count, tmap_file_t *fp) 
{
  size_t num_written = 0;

  switch(fp->c) {
    case TMAP_FILE_NO_COMPRESSION:
      num_written = fwrite(ptr, size, count, fp->fp);
      break;
#ifndef DISABLE_BZ2 
    case TMAP_FILE_BZ2_COMPRESSION:
      num_written = BZ2_bzwrite(fp->bz2, ptr, size*count) / size;
      break;
#endif
    case TMAP_FILE_GZ_COMPRESSION:
      num_written = gzwrite(fp->gz, ptr, size*count) / size;
      break;
    default:
      tmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(num_written != count) {
      tmap_error(NULL, Exit, WriteFileError);
  }

  return num_written;
}

int32_t 
tmap_file_vfprintf(tmap_file_t *fp, const char *format, va_list ap)
{
  int32_t n;

  if(TMAP_FILE_NO_COMPRESSION != fp->c) {
      tmap_error("compression not supported", Exit, OutOfRange);
  }

  n = vfprintf(fp->fp, format, ap);

  if(n < 0) {
      tmap_error("vfprintf failed", Exit, WriteFileError);
  }

  return n;
} 

int32_t 
tmap_file_fprintf(tmap_file_t *fp, const char *format, ...)
{
  int32_t n;
  va_list ap;

  if(NULL == fp) tmap_error("input file pointer was null", Exit, WriteFileError);

  va_start(ap, format);
  n = vfprintf(fp->fp, format, ap);
  va_end(ap);

  return n;
} 

int32_t 
tmap_file_printf(const char *format, ...)
{
  int32_t n;
  va_list ap;

  if(NULL == tmap_file_stdout) tmap_error("stdout file pointer was null", Exit, WriteFileError);

  va_start(ap, format);
  n = vfprintf(tmap_file_stdout->fp, format, ap);
  va_end(ap);

  return n;
} 

int32_t
tmap_file_fflush(tmap_file_t *fp, int32_t gz_flush)
{
  int32_t ret = EOF;

  switch(fp->c) {
    case TMAP_FILE_NO_COMPRESSION:
      ret = fflush(fp->fp);
      break;
#ifndef DISABLE_BZ2 
    case TMAP_FILE_BZ2_COMPRESSION:
      // do nothign
      ret = 0;
      break;
#endif
    case TMAP_FILE_GZ_COMPRESSION:
      if(Z_OK == gzflush(fp->gz, (0 == gz_flush) ? Z_SYNC_FLUSH : Z_FULL_FLUSH)) {
          ret = 0;
      }
      break;
    default:
      tmap_error("fp->c", Exit, OutOfRange);
      break;
  }

  if(0 != ret) {
      tmap_error(NULL, Exit, WriteFileError);
  }

  return ret;
}
