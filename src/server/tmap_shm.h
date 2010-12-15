#ifndef TMAP_SHM_H_
#define TMAP_SHM_H_

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <stdint.h>

#define TMAP_SHM_NOT_READY ~0xffaa6161
#define TMAP_SHM_READY 0xffaa6161
#define TMAP_SHM_DEAD  0xaabbccdd
#define TMAP_SHMGET_SLEEP 10
#define TMAP_SHMGET_RETRIES 10

/*! 
  Shared Memory Library
  */

/*! 
  @details  The server listings are stored as 0/1 bits, with 32-bits currently supported.
  They list if a given data structure is loaded into memory.
  */
enum {
    TMAP_SHM_LISTING_REFSEQ     = 0x1, /*!< the forward packed reference sequence  */
    TMAP_SHM_LISTING_REV_REFSEQ = 0x2, /*!< the forward packed reference sequence  */
    TMAP_SHM_LISTING_BWT        = 0x4, /*!< the forward BWT string */
    TMAP_SHM_LISTING_REV_BWT    = 0x8, /*!< the reverse BWT string */
    TMAP_SHM_LISTING_SA         = 0x10, /*!< the forward SA string */
    TMAP_SHM_LISTING_REV_SA     = 0x20 /*!< the reverse SA string */
};

/*! 
  @details  four bytes begin the shared memory for lazy synchronization
  */
typedef struct {
    key_t key;  /*!< the key of the shared memory  */
    size_t size;  /*!< the size of the shared memory */
    int32_t shmid;  /*!< the id of the shared memory */
    void *ptr;  /*!< pointer to the first byte of the shared memory */
    void *buf;  /*!< pointer to the first byte of the data stored in the shared memory */
    int32_t creator; /*!< 1 if this structure created the shared memory, 0 otherwise */
} tmap_shm_t;

/*! 
  @param  shm  pointer to the shared memory structure
  @return      the state of the shared memory
  */
inline uint32_t
tmap_shm_get_state(tmap_shm_t *shm);

/*! 
  @param  shm      pointer to the shared memory structure
  @param  listing  the listing to check
  @return          1 if the listing exists, 0 otherwise
  */
inline uint32_t
tmap_shm_listing_exists(tmap_shm_t *shm, uint32_t listing);

/*! 
  @param  shm      pointer to the shared memory structure
  @param  listing  the listing of the shared memory
  @param  size     size of the listing in bytes
  */
inline void
tmap_shm_add_listing(tmap_shm_t *shm, uint32_t listing, size_t size);

/*! 
  @param  shm      pointer to the shared memory structure
  @param  listing  the listing of the shared memory
  @return          the number of bytes from the start of the buffer, or SIZE_MAX if the listing does not exist
  */
inline size_t
tmap_shm_get_listing_bytes(tmap_shm_t *shm, uint32_t listing);

/*! 
  @param  shm      pointer to the shared memory structure
  @param  listing  the listing of the shared memory
  @return          pointer to the beginning of the packed listing
  */
inline uint8_t *
tmap_shm_get_buffer(tmap_shm_t *shm, uint32_t listing);

/*! 
  @param  shm  pointer to the shared memory structure
  */
inline void
tmap_shm_set_not_ready(tmap_shm_t *shm);
/*! 
  @param  shm  pointer to the shared memory structure
  */
inline void
tmap_shm_set_ready(tmap_shm_t *shm);

/*! 
  @param  shm  pointer to the shared memory structure
  */
inline void
tmap_shm_set_dead(tmap_shm_t *shm);

/*! 
  @param  key         the key of the shared memory 
  @param  size        the total size of the shared memory
  @param  create      1 if the process is to create the shared memory, 0 otherwise
  @return             a pointer to the initialized shared memory
  */
tmap_shm_t *
tmap_shm_init(key_t key, size_t size, int32_t create);

/*! 
  @param  shm    pointer to the shared memory structure
  @param  force  forces the shared memory to be destroyed
  */
void
tmap_shm_destroy(tmap_shm_t *shm, int32_t force);

#endif
