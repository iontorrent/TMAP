#ifndef FMAP_ERROR_H_
#define FMAP_ERROR_H_

#define BREAK_LINE "************************************************************\n"

#include <stdint.h>

/*! @header
  @abstract  Error handling routines.
 */

/*! 
  @abstract                  the type of action to be taken
  @constant  Exit            exit the program
  @constant  Warn            print a warning 
  @constant  LastActionType  dummy action type
  @discussion  the type of action to take upon the detection of an error
  */
enum {
    Exit, 
    Warn, 
    LastActionType
};

/*! 
  @abstract                   the type of error
  @constant  OutOfRange           value was out of range
  @constant  CommandLineArgument  improper command line argument
  @constant  ReallocMemory        memory re-allocation failure
  @constant  MallocMemory         memory allocation failure
  @constant  OpenFileError        could not open a file
  @constant  CloseFileError       could not close a file
  @constant  ReadFileError        could not read from a file
  @constant  WriteFileError       could not write from a file
  @constant  EndOfFile            reached the end-of-file prematurely
  @constant  ThreadError          error starting/joining threads
  @constant  SigInt               SIGINT signal caught
  @constant  SharedMemoryGet      could not get the shared memory
  @constant  SharedMemoryAttach   could not attach the shared memory
  @constant  SharedMemoryControl  could not control the shared memory
  @constant  SharedMemoryDetach   could not detach the shared memory
  @constant  SharedMemoryListing  could not find the listing in shared memory
  @constant  LastErrorType        dummy error type 
  @discussion  the type of error detected
  */
enum {
    OutOfRange=0, 
    CommandLineArgument,
    ReallocMemory,
    MallocMemory,
    OpenFileError,
    CloseFileError,
    ReadFileError,
    WriteFileError,
    EndOfFile,
    ThreadError,
    SigInt,
    SharedMemoryGet,
    SharedMemoryAttach,
    SharedMemoryControl,
    SharedMemoryDetach,
    SharedMemoryListing,
    LastErrorType,
};

/*! @function
  @abstract      checks if the value falls within the bounds
  @param  val    the value to be checked
  @param  lower  the lower integer value (inclusive)
  @param  upper  the upper integer value (inclusive)
  @discussion    throws a command line argument error if the value is not within the bounds
  */
void
fmap_error_cmd_check_int(int32_t val, int32_t lower, int32_t upper, char *option);

/*! @define
  @abstract              process an error based on the given action
  @param  variable_name  the variable name or value associated with the error
  @param  action_type    the action to be taken
  @param  error_type     the error type 
  */
#define fmap_error(variable_name, action_type, error_type) \
  (fmap_error_full(__FILE__, __LINE__, __func__, variable_name, action_type, error_type))

/*! @define
  @abstract              process an error based on the given action
  @param  function_name  the function name reporting the error
  @param  variable_name  the variable name or value associated with the error
  @param  action_type    the action to be taken
  @param  error_type     the error type 
  */
#define fmap_error1(function_name, variable_name, action_type, error_type) \
  (fmap_error_full(__FILE__, __LINE__, function_name, variable_name, action_type, error_type))

/*! @function
  @abstract              process an error based on the given action
  @param  function_name  the function name reporting the error
  @param  variable_name  the variable name or value associated with the error
  @param  action_type    the action to be taken
  @param  error_type     the error type 
  */
void 
fmap_error_full(const char *file, const unsigned int line, const char *function_name, const char *variable_name, int action_type, int error_type);

#endif
