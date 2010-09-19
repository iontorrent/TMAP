#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <stdint.h>
#include <unistd.h>

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "fmap_shm.h"

static int32_t
fmap_shmget(key_t key, size_t size, int32_t shmflg, int32_t create)
{
  int32_t shmid, i;

  if(0 == create) {
      // try a number of times before failing
      for(i=0,shmid=-1;shmid<0 && i<FMAP_SHMGET_RETRIES-1;i++) {
          if(0 <= (shmid = shmget(key, size, shmflg))) {
              return shmid;
          }
          // sleep and retry
          sleep(FMAP_SHMGET_SLEEP);
      }
  }
  if((shmid = shmget(key, size, shmflg)) < 0) {
      fmap_error(NULL, Exit, SharedMemoryGet);
  }

  return shmid;
}

static void *
fmap_shmat(int32_t shmid, const void *shmaddr, int32_t shmflg)
{
  void *shm = NULL;

  if((shm = shmat(shmid, shmaddr, shmflg)) == (char*)-1) {
      fmap_error(NULL, Exit, SharedMemoryAttach);
  }

  return shm;
}

static int32_t
fmap_shmctl(int32_t shmid, int32_t cmd, struct shmid_ds *buf)
{
  if(shmctl(shmid, cmd, buf) < 0) {
      fmap_error(NULL, Exit, SharedMemoryControl);
  }

  return 0;
}

static int32_t
fmap_shmdt(const void *shmaddr)
{
  if(shmdt(shmaddr) < 0) {
      fmap_error(NULL, Exit, SharedMemoryDetach);
  }
  return 0;
}

static inline void
fmap_shm_set_state(fmap_shm_t *shm, uint32_t state)
{
  ((volatile uint32_t*)shm->ptr)[0] = state;
}

inline uint32_t 
fmap_shm_get_state(fmap_shm_t *shm)
{
  return ((volatile uint32_t*)shm->ptr)[0];
}

inline void
fmap_shm_set_not_ready(fmap_shm_t *shm)
{
  fmap_shm_set_state(shm, FMAP_SHM_NOT_READY);
}

inline void
fmap_shm_set_ready(fmap_shm_t *shm)
{
  fmap_shm_set_state(shm, FMAP_SHM_READY);
}

inline void
fmap_shm_set_dead(fmap_shm_t *shm)
{
  fmap_shm_set_state(shm, FMAP_SHM_DEAD);
}

static inline uint32_t 
fmap_shm_get_listing(fmap_shm_t *shm)
{
  return ((volatile uint32_t*)shm->ptr)[1];
}

inline uint32_t
fmap_shm_listing_exists(fmap_shm_t *shm, uint32_t listing)
{
  if(listing == (listing & fmap_shm_get_listing(shm))) {
      return 1;
  }
  else {
      return 0;
  }
}

inline void
fmap_shm_add_listing(fmap_shm_t *shm, uint32_t listing, size_t size)
{
  ((volatile uint32_t*)shm->ptr)[1] |= listing;
  ((size_t*)(((volatile uint32_t*)shm->ptr) + 2))[fmap_log2(listing)] = size;
}

inline size_t
fmap_shm_get_listing_bytes(fmap_shm_t *shm, uint32_t listing)
{
  int32_t i;
  size_t s = 0;

  if(0 == fmap_shm_listing_exists(shm, listing)) return SIZE_MAX;

  for(i=1;i<listing;i<<=1) {
      if(1 == fmap_shm_listing_exists(shm, i)) {
          s += ((size_t*)(((volatile uint32_t*)shm->ptr) + 2))[fmap_log2(i)];
      }
  }
  return s;
}

inline uint8_t *
fmap_shm_get_buffer(fmap_shm_t *shm, uint32_t listing)
{
  size_t s = fmap_shm_get_listing_bytes(shm, listing);
  if(SIZE_MAX == s) return NULL;
  return (uint8_t*)(((volatile uint8_t*)(shm->buf)) + s);
}

fmap_shm_t *
fmap_shm_init(key_t key, size_t size, int32_t create)
{
  fmap_shm_t *shm = NULL;
  int32_t shmflg = 0;
  struct shmid_ds buf;

  shm = fmap_calloc(1, sizeof(fmap_shm_t), "shm");
  shm->key = key;
  shm->size = size;

  if(1 == create) {
      shm->size += sizeof(uint32_t); // add for synchronization
      shm->size += sizeof(uint32_t); // add for on/off bits for listing what is in memory
      shm->size += 32*sizeof(size_t); // add for the byte size of each listing
      shmflg = IPC_CREAT | IPC_EXCL | 0666;
  }
  else {
      shmflg = 0666;
  }

  // get the shared memory id
  shm->shmid = fmap_shmget(shm->key, shm->size, shmflg, create);

  // attach the shared memory
  shm->ptr = fmap_shmat(shm->shmid, NULL, 0);
  shm->buf = ((char*)shm->ptr);
  shm->buf += sizeof(uint32_t); // synchronization 
  shm->buf += sizeof(uint32_t) + 32*sizeof(size_t); // listings

  if(1 == create) {
      // check that the current process created the shared memory
      fmap_shmctl(shm->shmid, IPC_STAT, &buf);
      if(buf.shm_cpid != getpid() || FMAP_SHM_READY == fmap_shm_get_state(shm)) {
          fmap_error("shared memory was not created by the current process", Exit, OutOfRange);
      }
      fmap_shm_set_not_ready(shm);
      if(buf.shm_segsz != shm->size) {
          fmap_error("shared memory size does not match the expected size", Exit, OutOfRange);
      }
  }
  else {
      // set the size
      shm->size = buf.shm_segsz;
  }

  return shm;
}

void
fmap_shm_destroy(fmap_shm_t *shm, int32_t force)
{
  struct shmid_ds buf;

  // detach the shared memory
  fmap_shmdt(shm->ptr);

  fmap_shmctl(shm->shmid, IPC_STAT, &buf);
  if(1 == force || buf.shm_cpid == getpid()) { // delete the shared memory
      fmap_shmctl(shm->shmid, IPC_RMID, NULL);
  }
  free(shm);
}
