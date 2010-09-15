#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fmap_error.h"

static char error_string[][64] =
{ 
  "the value is out of range",
  "command line argument",
  "could not re-allocate memory",
  "could not allocate memory",
  "could not open the file",
  "could not close the file",
  "could not read from the file",
  "could not write to the file",
  "encountered early end-of-file",
  "error running the threads",
  "last error type"
};	   

static char action_string[][20] =
{"fatal error", "warning", "LastActionType"};

void
fmap_error_cmd_check_int(int32_t val, int32_t lower, int32_t upper, char *option)
{
  if(val < lower || upper < val) {
      char str[1024] = "\0";
      strcpy(str, "option ");
      strcat(str, option);
      strcat(str, " out of range");
      fmap_error(str, Exit, CommandLineArgument);
  }
} 

void 
fmap_error_full(const char *file, const unsigned int line, const char *function_name, const char *variable_name, int action_type, int error_type) 
{
  // Note: using regular-old fprintf since fmap_file_fprintf will call this
  // function (infinite recursion)

  if(NULL == variable_name) {
      fprintf(stderr, "\n%s:%u: in function \"%s\"\n%s: %s\n",
              file, line, function_name, action_string[action_type], error_string[error_type]);
  }
  else {
      fprintf(stderr, "\n%s:%u: in function \"%s\"\n%s\n%s: %s\n",
              file, line, function_name, variable_name, action_string[action_type], error_string[error_type]);
  }

  if(error_type == ReadFileError || 
     error_type == OpenFileError || 
     error_type == WriteFileError ||
     error_type == EndOfFile) {
      perror("The file stream error was");
  }

  switch(action_type) {
    case Exit: 
#ifdef PACKAGE_BUGREPORT
      fprintf(stderr, "Please report bugs to %s.\n", PACKAGE_BUGREPORT);
#endif
      exit(EXIT_FAILURE); 
      break; /* Not necessary actually! */
    case Warn:
      fprintf(stderr, "\n");
      break;
    default:
      exit(EXIT_FAILURE); 
      break;
  }
}
