/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_PROGRESS_H
#define TMAP_PROGRESS_H

#include <time.h>

/*! 
  Progress Messaging 
  */

/*! 
  @param  command  the command string
  */
void
tmap_progress_set_command(const char *command);

/*! 
  Set the start time of the clock
  */
void
tmap_progress_set_start_time();

/*! 
  @param  verbosity  the verbosity level
  */
void
tmap_progress_set_verbosity(int32_t verbosity);

/*! 
  @return  the verbosity level
  */
int32_t
tmap_progress_get_verbosity();

/*! 
  @param  format      the format for the message
  @param  ...         the arguments to fill in the format
  */
void
tmap_progress_print(const char *format, ...);

/*! 
  @param  format      the format for the message
  @param  start_time  the execution start time 
  @param  ...         the arguments to fill in the format
  */
void
tmap_progress_print1(const char *format, clock_t start_time, ...);

/*! 
  @param  format   the format for the message
  @param  ...      the arguments to fill in the format
  */
void
tmap_progress_print2(const char *format, ...);

#endif
