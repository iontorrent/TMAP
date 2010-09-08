#ifndef FMAP_ERROR_H_
#define FMAP_ERROR_H_

#define BREAK_LINE "************************************************************\n"

/*! @enum
  @abstract              the type of action to be taken
  @field  Exit            exit the program
  @field  Warn            print a warning 
  @field  LastActionType  dummy action type
  */
enum {Exit, Warn, LastActionType};

/*! @enum
  @abstract                   the type of error
  @field  OutOfRange           value was out of range
  @field  CommandLineArgument  improper command line argument
  @field  ReallocMemory        memory re-allocation failure
  @field  MallocMemory         memory allocation failure
  @field  OpenFileError        could not open a file
  @field  CloseFileError       could not close a file
  @field  ReadFileError        could not read from a file
  @field  WriteFileError       could not write from a file
  @field  EndOfFile            reached the end-of-file prematurely
  @field  ThreadError          error starting/joining threads
  @field  LastErrorType        dummy error type 
  */
enum {
    OutOfRange=0, 
    CommandLineArgument,
    ReallocMemory,
    MallocMemory,
    OpenFileError,
    CloseFileEror,
    ReadFileError,
    WriteFileError,
    EndOfFile,
    ThreadError,
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

/*! @function
  @abstract              process an error based on the given action
  @param  variable_name  the variable name or value associated with the error
  @param  action_type    the action to be taken
  @param  error_type     the error type 
  */
#define fmap_error(_variable_name, _action_type, _error_type) \
  (fmap_error_full(__FILE__, __LINE__, __func__, _variable_name, _action_type, _error_type))

/*! @function
  @abstract              process an error based on the given action
  @param  function_name  the function name reporting the error
  @param  variable_name  the variable name or value associated with the error
  @param  action_type    the action to be taken
  @param  error_type     the error type 
  */
#define fmap_error1(_function_name, _variable_name, _action_type, _error_type) \
  (fmap_error_full(__FILE__, __LINE__, _function_name, _variable_name, _action_type, _error_type))

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
