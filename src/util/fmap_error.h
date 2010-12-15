#ifndef TMAP_ERROR_H_
#define TMAP_ERROR_H_

#define BREAK_LINE "************************************************************\n"

#include <stdint.h>

/*! 
  Error handling routines.
 */

/*! 
  the type of action to be taken
  @details  the type of action to take upon the detection of an error
  */
enum {
    Exit,  /*!< exit the program */
    Warn,  /*!< print a warning  */
    LastActionType /*!< dummy action type */
};

/*! 
  the type of error
  @details  the type of error detected
  */
enum {
    OutOfRange=0,  /*!< value was out of range */
    CommandLineArgument, /*!< improper command line argument */
    ReallocMemory, /*!< memory re-allocation failure */
    MallocMemory, /*!< memory allocation failure */
    OpenFileError, /*!< could not open a file */
    CloseFileError, /*!< could not close a file */
    ReadFileError, /*!< could not read from a file */
    WriteFileError, /*!< could not write from a file */
    EndOfFile, /*!< reached the end-of-file prematurely */
    ThreadError, /*!< error starting/joining threads */
    SigInt, /*!< SIGINT signal caught */
    SharedMemoryGet, /*!< could not get the shared memory */
    SharedMemoryAttach, /*!< could not attach the shared memory */
    SharedMemoryControl, /*!< could not control the shared memory */
    SharedMemoryDetach, /*!< could not detach the shared memory */
    SharedMemoryListing, /*!< could not find the listing in shared memory */
    LastErrorType, /*!< dummy error type  */
};

/*! 
  checks if the value falls within the bounds
  @param  val     the value to be checked
  @param  lower   the lower integer value (inclusive)
  @param  upper   the upper integer value (inclusive)
  @param  option  the option being checked 
  @details        throws a command line argument error if the value is not within the bounds
  */
void
tmap_error_cmd_check_int(int32_t val, int32_t lower, int32_t upper, char *option);

/*! 
  process an error based on the given action
  @param  variable_name  the variable name or value associated with the error
  @param  action_type    the action to be taken
  @param  error_type     the error type 
  */
#define tmap_error(variable_name, action_type, error_type) \
  (tmap_error_full(__FILE__, __LINE__, __func__, variable_name, action_type, error_type))

/*! 
  process an error based on the given action
  @param  function_name  the function name reporting the error
  @param  variable_name  the variable name or value associated with the error
  @param  action_type    the action to be taken
  @param  error_type     the error type 
  */
#define tmap_error1(function_name, variable_name, action_type, error_type) \
  (tmap_error_full(__FILE__, __LINE__, function_name, variable_name, action_type, error_type))

/*! 
  process an error based on the given action
  @param  file            the calling file (use __FILE__)
  @param  line           the line number in the calling file (use __LINE__)
  @param  function_name  the function name reporting the error
  @param  variable_name  the variable name or value associated with the error
  @param  action_type    the action to be taken
  @param  error_type     the error type 
  */
void 
tmap_error_full(const char *file, const unsigned int line, const char *function_name, const char *variable_name, int action_type, int error_type);

#endif
