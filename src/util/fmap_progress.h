#ifndef FMAP_PROGRESS_H_
#define FMAP_PROGRESS_H_

#include <time.h>

/*! @header
  @abstract  Progress Messaging 
  */

/*! @function
  @abstract
  @param  command  the command string
  */
void
fmap_progress_set_command(const char *command);

/*! @function
  @abstract
  @param  start_time  the start time of the clock
  */
void
fmap_progress_set_start_time(clock_t start_time);

/*! @function
  @abstract
  @param  verbosity  the verbosity level
  */
void
fmap_progress_set_verbosity(int32_t verbosity);

/*! @function
  @abstract
  @param  format      the format for the message
  @param  ...         the arguments to fill in the format
  */
void
fmap_progress_print(const char *format, ...);

/*! @function
  @abstract
  @param  format      the format for the message
  @param  start_time  the execution start time 
  @param  ...         the arguments to fill in the format
  */
void
fmap_progress_print1(const char *format, clock_t start_time, ...);

/*! @function
  @abstract
  @param  format   the format for the message
  @param  ...      the arguments to fill in the format
  */
void
fmap_progress_print2(const char *format, ...);

#endif
