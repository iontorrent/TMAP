#ifndef FMAP_SERVER_H_
#define FMAP_SERVER_H_

/*! @header
  @abstract  A Server for Storing Reference Data
 */

/*! @define
  @abstract  How long the server should sleep
 */
#define FMAP_SERVER_SLEEP 5
//#define FMAP_SERVER_SLEEP 15

/*! @enum  The Server Command
  @constant  FMAP_SERVER_UNKNOWN  Unknown command
  @constant  FMAP_SERVER_START    Start command
  @constant  FMAP_SERVER_STOP     Stop command
  @constant  FMAP_SERVER_KILL     Kill command
  @discussion  The command to the server.
 */
enum {
    FMAP_SERVER_UNKNOWN = -1,
    FMAP_SERVER_START   = 0,
    FMAP_SERVER_STOP    = 1,
    FMAP_SERVER_KILL    = 2
};

/*! @function
  @abstract
  @param  signal  the signal to catch
  @discussion  use this function to catch a SIGINT signal
*/
void
fmap_server_sigint(int signal);

/*! @function
  @abstract
  @param  shm  the current shared memory
  */
void
fmap_server_set_sigint(fmap_shm_t *shm);

/*! @function
  @abstract         starts a server and loads in reference data
  @param  fn_fasta  the file name of the fasta file
  @param  key       the server key
  @param  listing   the server listing to load
  */
void
fmap_server_start(char *fn_fasta, key_t key, uint32_t listing);

/*! @function
  @abstract    stops the server with the given key
  @param  key  the server key
  */
void
fmap_server_stop(key_t key);

/*! @function
  @abstract    kills the server with the given key
  @param  key  the server key
  @discussion  the shared memory will be destroyed by the calling process 
  */
void
fmap_server_kill(key_t key);

/*! @function
  @abstract       parses the server command
  @param  optarg  the command to parse
  @return         the integer server command
  */
int32_t
fmap_server_get_command_int(char *optarg);

/*! @function
  @abstract     main-like function for 'fmap server'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
*/
int
fmap_server_main(int argc, char *argv[]);

#endif
