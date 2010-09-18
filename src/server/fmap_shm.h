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
  @param  key     the key of the shared memory 
  @param  size    the size of the shared memory
  @param  create  1 if the process is to create the shared memory, 0 otherwise
  @return         a pointer to the initialized shared memory
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
