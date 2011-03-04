/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <config.h>

#include "tmap_error.h"
#include "../io/tmap_file.h"
#include "tmap_progress.h"

static char tmap_progress_command[1024]="(null)";
static clock_t tmap_progress_start_time=0;
static int32_t tmap_progress_verbosity=0;

void 
tmap_progress_set_command(const char *command)
{
  strcpy(tmap_progress_command, PACKAGE);
  strcat(tmap_progress_command, " ");
  strcat(tmap_progress_command, command);
}

void 
tmap_progress_set_start_time(clock_t start_time)
{
  tmap_progress_start_time = start_time;
}

void
tmap_progress_set_verbosity(int32_t verbosity)
{
  tmap_progress_verbosity = verbosity;
}

static void 
tmap_progress_vprint1(const char *format, clock_t start_time, va_list ap)
{
  static char tmap_progress_format[2048]="\0";
  
  if(0 != tmap_progress_verbosity) {
      if(0 <= (float)start_time) {
          if(sprintf(tmap_progress_format, "[%s] %.2f sec: ", tmap_progress_command, (float)(clock()-start_time) / CLOCKS_PER_SEC) < 0) {
              tmap_error(NULL, Exit, OutOfRange);
          }
      }
      else {
          if(sprintf(tmap_progress_format, "[%s] ", tmap_progress_command) < 0) {
              tmap_error(NULL, Exit, OutOfRange);
          }
      }
      strcat(tmap_progress_format, format);
      strcat(tmap_progress_format, ".\n");

      tmap_file_vfprintf(tmap_file_stderr, tmap_progress_format, ap);
  }
}

void 
tmap_progress_print(const char *format, ...)
{
  va_list ap;
  if(0 != tmap_progress_verbosity) {
      va_start(ap, format);
      tmap_progress_vprint1(format, -1, ap);
      va_end(ap);
  }
}

void 
tmap_progress_print1(const char *format, clock_t start_time, ...)
{
  va_list ap;
  if(0 != tmap_progress_verbosity) {
      va_start(ap, start_time);
      tmap_progress_vprint1(format, start_time, ap);
      va_end(ap);
  }
}

inline void 
tmap_progress_print2(const char *format, ...)
{
  va_list ap;
  if(0 != tmap_progress_verbosity) {
      va_start(ap, format);
      tmap_progress_vprint1(format, tmap_progress_start_time, ap);
      va_end(ap);
  }
}
