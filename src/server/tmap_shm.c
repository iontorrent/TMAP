/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <stdint.h>
#include <unistd.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_progress.h"
#include "tmap_shm.h"

static int32_t
tmap_shmget(key_t key, size_t size, int32_t shmflg, int32_t create)
{
  int32_t shmid, i;

  if(0 == create) {
      // try a number of times before failing
      for(i=0,shmid=-1;shmid<0 && i<TMAP_SHMGET_RETRIES-1;i++) {
          if(0 <= (shmid = shmget(key, size, shmflg))) {
              return shmid;
          }
          tmap_progress_print("could not get shared memory, %d more %s", 
                              TMAP_SHMGET_RETRIES-i-1,
                              (1 != TMAP_SHMGET_RETRIES-i-1) ? "retries" : "retry");
          tmap_progress_print("retrying in %d seconds", TMAP_SHMGET_SLEEP);
          // sleep and retry
          sleep(TMAP_SHMGET_SLEEP);
      }
  }
  if((shmid = shmget(key, size, shmflg)) < 0) {
      tmap_error(NULL, Exit, SharedMemoryGet);
  }

  return shmid;
}

static void *
tmap_shmat(int32_t shmid, const void *shmaddr, int32_t shmflg)
{
  void *shm = NULL;

  if((shm = shmat(shmid, shmaddr, shmflg)) == (char*)-1) {
      tmap_error(NULL, Exit, SharedMemoryAttach);
  }

  return shm;
}

static int32_t
tmap_shmctl(int32_t shmid, int32_t cmd, struct shmid_ds *buf)
{
  if(shmctl(shmid, cmd, buf) < 0) {
      tmap_error(NULL, Exit, SharedMemoryControl);
  }

  return 0;
}

static int32_t
tmap_shmdt(const void *shmaddr)
{
  if(shmdt(shmaddr) < 0) {
      tmap_error(NULL, Exit, SharedMemoryDetach);
  }
  return 0;
}

static inline void
tmap_shm_set_state(tmap_shm_t *shm, uint32_t state)
{
  ((volatile uint32_t*)shm->ptr)[0] = state;
}

inline uint32_t 
tmap_shm_get_state(tmap_shm_t *shm)
{
  return ((volatile uint32_t*)shm->ptr)[0];
}

inline void
tmap_shm_set_not_ready(tmap_shm_t *shm)
{
  tmap_shm_set_state(shm, TMAP_SHM_NOT_READY);
}

inline void
tmap_shm_set_ready(tmap_shm_t *shm)
{
  tmap_shm_set_state(shm, TMAP_SHM_READY);
}

inline void
tmap_shm_set_dead(tmap_shm_t *shm)
{
  tmap_shm_set_state(shm, TMAP_SHM_DEAD);
}

static inline uint32_t 
tmap_shm_get_listing(tmap_shm_t *shm)
{
  return ((volatile uint32_t*)shm->ptr)[1];
}

inline uint32_t
tmap_shm_listing_exists(tmap_shm_t *shm, uint32_t listing)
{
  if(listing == (listing & tmap_shm_get_listing(shm))) {
      return 1;
  }
  else {
      return 0;
  }
}

inline void
tmap_shm_add_listing(tmap_shm_t *shm, uint32_t listing, size_t size)
{
  ((volatile uint32_t*)shm->ptr)[1] |= listing;
  ((size_t*)(((volatile uint32_t*)shm->ptr) + 2))[tmap_log2(listing)] = size;
}

inline size_t
tmap_shm_get_listing_bytes(tmap_shm_t *shm, uint32_t listing)
{
  int32_t i;
  size_t s = 0;

  if(0 == tmap_shm_listing_exists(shm, listing)) return SIZE_MAX;

  for(i=1;i<listing;i<<=1) {
      if(1 == tmap_shm_listing_exists(shm, i)) {
          s += ((size_t*)(((volatile uint32_t*)shm->ptr) + 2))[tmap_log2(i)];
      }
  }
  return s;
}

inline uint8_t *
tmap_shm_get_buffer(tmap_shm_t *shm, uint32_t listing)
{
  size_t s = tmap_shm_get_listing_bytes(shm, listing);
  if(SIZE_MAX == s) return NULL;
  return (uint8_t*)(((volatile uint8_t*)(shm->buf)) + s);
}

tmap_shm_t *
tmap_shm_init(key_t key, size_t size, int32_t create)
{
  tmap_shm_t *shm = NULL;
  int32_t i, shmflg = 0;
  struct shmid_ds buf;

  shm = tmap_calloc(1, sizeof(tmap_shm_t), "shm");
  shm->key = key;
  shm->size = size;

  if(1 == create) {
      shm->size += sizeof(uint32_t); // add for synchronization
      shm->size += sizeof(uint32_t); // add for on/off bits for listing what is in memory
      shm->size += 32*sizeof(size_t); // add for the byte size of each listing
      shmflg = IPC_CREAT | IPC_EXCL | 0666;
      shm->creator = 1;
  }
  else {
      shmflg = 0666;
      shm->creator = 0;
  }

  // get the shared memory id
  shm->shmid = tmap_shmget(shm->key, shm->size, shmflg, create);

  // attach the shared memory
  shm->ptr = tmap_shmat(shm->shmid, NULL, 0);
  shm->buf = ((char*)shm->ptr);
  shm->buf += sizeof(uint32_t); // synchronization 
  shm->buf += sizeof(uint32_t) + 32*sizeof(size_t); // listings

  tmap_shmctl(shm->shmid, IPC_STAT, &buf);
  if(1 == create) {
      // check that the current process created the shared memory
      if(buf.shm_cpid != getpid() || TMAP_SHM_READY == tmap_shm_get_state(shm)) {
          tmap_error("shared memory was not created by the current process", Exit, OutOfRange);
      }
      tmap_shm_set_not_ready(shm);
      if(buf.shm_segsz != shm->size) {
          tmap_error("shared memory size does not match the expected size", Exit, OutOfRange);
      }
  }
  else {
      // try a number of times before failing
      for(i=0;i<TMAP_SHMGET_RETRIES;i++) {
          if(TMAP_SHM_READY == tmap_shm_get_state(shm)) {
              break;
          }
          tmap_progress_print("shared memory not ready, %d more retries", TMAP_SHMGET_RETRIES-i-1);
          tmap_progress_print("retrying in %d seconds", TMAP_SHMGET_SLEEP);
          // sleep and retry
          sleep(TMAP_SHMGET_SLEEP);
      }
      if(TMAP_SHMGET_RETRIES == i) {
          tmap_error("shared memory did not become available", Exit, SharedMemoryGet);
      }
      
      // set the size
      shm->size = buf.shm_segsz;
  }

  return shm;
}

void
tmap_shm_destroy(tmap_shm_t *shm, int32_t force)
{
  struct shmid_ds buf;

  // detach the shared memory
  tmap_shmdt(shm->ptr);

  tmap_shmctl(shm->shmid, IPC_STAT, &buf);
  if(1 == force || 1 == shm->creator) { // delete the shared memory
      tmap_shmctl(shm->shmid, IPC_RMID, NULL);
  }
  free(shm);
}
