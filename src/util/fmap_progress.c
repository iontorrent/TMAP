#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <config.h>

#include "fmap_error.h"
#include "../io/fmap_file.h"
#include "fmap_progress.h"

static char fmap_progress_command[1024]="(null)";
static clock_t fmap_progress_start_time=0;
static int32_t fmap_progress_verbosity=0;

void 
fmap_progress_set_command(const char *command)
{
  strcpy(fmap_progress_command, PACKAGE);
  strcat(fmap_progress_command, " ");
  strcat(fmap_progress_command, command);
}

void 
fmap_progress_set_start_time(clock_t start_time)
{
  fmap_progress_start_time = start_time;
}

void
fmap_progress_set_verbosity(int32_t verbosity)
{
  fmap_progress_verbosity = verbosity;
}

static void 
fmap_progress_vprint1(const char *format, clock_t start_time, va_list ap)
{
  static char fmap_progress_format[2048]="\0";
  
  if(0 != fmap_progress_verbosity) {
      if(0 < (float)start_time) {
          if(sprintf(fmap_progress_format, "[%s] %.2f sec: ", fmap_progress_command, (float)(clock()-start_time) / CLOCKS_PER_SEC) < 0) {
              fmap_error(NULL, Exit, OutOfRange);
          }
      }
      else {
          if(sprintf(fmap_progress_format, "[%s] ", fmap_progress_command) < 0) {
              fmap_error(NULL, Exit, OutOfRange);
          }
      }
      strcat(fmap_progress_format, format);
      strcat(fmap_progress_format, ".\n");

      fmap_file_vfprintf(fmap_file_stderr, fmap_progress_format, ap);
  }
}

void 
fmap_progress_print(const char *format, ...)
{
  va_list ap;
  if(0 != fmap_progress_verbosity) {
      va_start(ap, format);
      fmap_progress_vprint1(format, 0, ap);
      va_end(ap);
  }
}

void 
fmap_progress_print1(const char *format, clock_t start_time, ...)
{
  va_list ap;
  if(0 != fmap_progress_verbosity) {
      va_start(ap, start_time);
      fmap_progress_vprint1(format, start_time, ap);
      va_end(ap);
  }
}

inline void 
fmap_progress_print2(const char *format, ...)
{
  va_list ap;
  if(0 != fmap_progress_verbosity) {
      va_start(ap, format);
      fmap_progress_vprint1(format, fmap_progress_start_time, ap);
      va_end(ap);
  }
}
