#ifndef TMAP_SERVER_H_
#define TMAP_SERVER_H_

/*! 
  A Server for Storing Reference Data
 */

/*! 
  How long the server should sleep
 */
#define TMAP_SERVER_SLEEP 5
//#define TMAP_SERVER_SLEEP 15

/*! 
  @details  The command to the server.
 */
enum {
    TMAP_SERVER_UNKNOWN = -1, /*!< Unknown command */
    TMAP_SERVER_START   = 0, /*!< Start command */
    TMAP_SERVER_STOP    = 1, /*!< Stop command */
    TMAP_SERVER_KILL    = 2 /*!< Kill command */
};

/*! 
  @param  signal  the signal to catch
  @details        use this function to catch a SIGINT signal
*/
void
tmap_server_sigint(int signal);
/*! 
  @param  shm  the current shared memory
  */
void
tmap_server_set_sigint(tmap_shm_t *shm);

/*! 
  starts a server and loads in reference data
  @param  fn_fasta  the file name of the fasta file
  @param  key       the server key
  @param  listing   the server listing to load
  */
void
tmap_server_start(char *fn_fasta, key_t key, uint32_t listing);

/*! 
  stops the server with the given key
  @param  key  the server key
  */
void
tmap_server_stop(key_t key);

/*! 
  kills the server with the given key
  @param  key  the server key
  @details     the shared memory will be destroyed by the calling process 
  */
void
tmap_server_kill(key_t key);

/*! 
  parses the server command
  @param  optarg  the command to parse
  @return         the integer server command
  */
int32_t
tmap_server_get_command_int(char *optarg);

/*! 
  main-like function for 'tmap server'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
*/
int
tmap_server_main(int argc, char *argv[]);
#endif
