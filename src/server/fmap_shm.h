#ifndef FMAP_SHM_H_
#define FMAP_SHM_H_

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <stdint.h>

#define FMAP_SHM_NOT_READY ~0xffaa6161
#define FMAP_SHM_READY 0xffaa6161
#define FMAP_SHM_DEAD  0xaabbccdd
#define FMAP_SHMGET_SLEEP 10
#define FMAP_SHMGET_RETRIES 10

/*! @header
  @abstract  Shared Memory Library
  */

/*! 
  @constant  FMAP_SHM_LISTING_REFSEQ      the forward packed reference sequence 
  @constant  FMAP_SHM_LISTING_REV_REFSEQ  the forward packed reference sequence 
  @constant  FMAP_SHM_LISTING_BWT         the forward BWT string
  @constant  FMAP_SHM_LISTING_REV_BWT     the reverse BWT string
  @constant  FMAP_SHM_LISTING_SA          the forward SA string
  @constant  FMAP_SHM_LISTING_REV_SA      the reverse SA string
  @discussion  The server listings are stored as 0/1 bits, with 32-bits currently supported.
  They list if a given data structure is loaded into memory.
  */
enum {
    FMAP_SHM_LISTING_REFSEQ     = 0x1,
    FMAP_SHM_LISTING_REV_REFSEQ = 0x2,
    FMAP_SHM_LISTING_BWT        = 0x4,
    FMAP_SHM_LISTING_REV_BWT    = 0x8,
    FMAP_SHM_LISTING_SA         = 0x10,
    FMAP_SHM_LISTING_REV_SA     = 0x20
};

/*! @typedef
  @abstract
  @field  key    the key of the shared memory 
  @field  size   the size of the shared memory
  @field  shmid  the id of the shared memory
  @field  ptr    pointer to the first byte of the shared memory
  @field  buf    pointer to the first byte of the data stored in the shared memory
  @discussion   four bytes begin the shared memory for lazy synchronization
  */
typedef struct {
    key_t key;
    size_t size;
    int32_t shmid;
    void *ptr;
    void *buf;
} fmap_shm_t;

/*! @function
  @abstract
  @param  shm  pointer to the shared memory structure
  @return      the state of the shared memory
  */
inline uint32_t
fmap_shm_get_state(fmap_shm_t *shm);

/*! @function
  @abstract
  @param  shm      pointer to the shared memory structure
  @param  listing  the listing to check
  @return          1 if the listing exists, 0 otherwise
  */
inline uint32_t
fmap_shm_listing_exists(fmap_shm_t *shm, uint32_t listing);

/*! @function
  @abstract
  @param  shm      pointer to the shared memory structure
  @param  listing  the listing of the shared memory
  @param  size     size of the listing in bytes
  */
inline void
fmap_shm_add_listing(fmap_shm_t *shm, uint32_t listing, size_t size);

/*! @function
  @abstract
  @param  shm      pointer to the shared memory structure
  @param  listing  the listing of the shared memory
  @return          the number of bytes from the start of the buffer, or SIZE_MAX if the listing does not exist
  */
inline size_t
fmap_shm_get_listing_bytes(fmap_shm_t *shm, uint32_t listing);

/*! @function
  @abstract
  @param  shm      pointer to the shared memory structure
  @param  listing  the listing of the shared memory
  @return          pointer to the beginning of the packed listing
  */
inline uint8_t *
fmap_shm_get_buffer(fmap_shm_t *shm, uint32_t listing);

/*! @function
  @abstract
  @param  shm  pointer to the shared memory structure
  */
inline void
fmap_shm_set_not_ready(fmap_shm_t *shm);
/*! @function
  @abstract
  @param  shm  pointer to the shared memory structure
  */
inline void
fmap_shm_set_ready(fmap_shm_t *shm);

/*! @function
  @abstract
  @param  shm  pointer to the shared memory structure
  */
inline void
fmap_shm_set_dead(fmap_shm_t *shm);

/*! @function
  @abstract
  @param  key         the key of the shared memory 
  @param  size        the total size of the shared memory
  @param  create      1 if the process is to create the shared memory, 0 otherwise
  @return             a pointer to the initialized shared memory
  */
fmap_shm_t *
fmap_shm_init(key_t key, size_t size, int32_t create);

/*! @function
  @abstract
  @param  shm    pointer to the shared memory structure
  @param  force  forces the shared memory to be destroyed
  */
void
fmap_shm_destroy(fmap_shm_t *shm, int32_t force);

#endif
