/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/*-
 * Copyright 1997, 1998-2003 John-Mark Gurney.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *	$Id: fib.h,v 1.10 2003/01/14 10:11:30 jmg Exp $
 *
 */

/* Nils Homer - modified */

#ifndef TMAP_FIBHEAP_H_
#define TMAP_FIBHEAP_H_

/*! 
  */

/*! 
  a node in the fibonacci heap
 */
typedef struct __tmap_fibheap_element_t {
    int	tmap_fibheap_element_degree;
    int	tmap_fibheap_element_mark;
    struct __tmap_fibheap_element_t *tmap_fibheap_element_p;
    struct __tmap_fibheap_element_t *tmap_fibheap_element_child;
    struct __tmap_fibheap_element_t *tmap_fibheap_element_left;
    struct __tmap_fibheap_element_t *tmap_fibheap_element_right;
    int	tmap_fibheap_element_key;
    void	*tmap_fibheap_element_data;
} tmap_fibheap_element_t; 

/*! 
  a fibonacci heap
 */
typedef struct {
    int	(*tmap_fibheap_cmp_fnct)(void *, void *);
    int	tmap_fibheap_n;
    int	tmap_fibheap_Dl;
    tmap_fibheap_element_t **tmap_fibheap_cons;
    tmap_fibheap_element_t *tmap_fibheap_min;
    tmap_fibheap_element_t *tmap_fibheap_root;
    void *tmap_fibheap_neginf;
    int	tmap_fibheap_keys;
} tmap_fibheap_t;

/*! 
  the comparison prototype 
 */
typedef int (*tmap_fibheap_voidcmp)(void *, void *);

/*! 
  create a key heap
  @return    pointer to the initialized key heap
*/
tmap_fibheap_t *
tmap_fibheap_makekeyheap();
/*! 
  inserts an element with the given key into a key heap
  @param  h     pointer to the heap structure
  @param  key   the key of the element
  @param  data  the data to insert
  @return       pointer to the element inserted 
*/
tmap_fibheap_element_t *
tmap_fibheap_insertkey(tmap_fibheap_t *h, int key, void *data);
/*! 
  gets the minimum element's key from a key heap
  @param  h  pointer to the heap structure
  @return    the minimum element's key, or INT_MIN if the key is empty
*/
int 
tmap_fibheap_minkey(tmap_fibheap_t *h);
/*! 
  changes a given element's key in a key heap
  @param  h    pointer to the heap structure
  @param  x    the element whos key to update
  @param  key  the new key
  @return      the element's original key
*/
int 
tmap_fibheap_replacekey(tmap_fibheap_t *h, tmap_fibheap_element_t *x, int key);
/*! 
  changes a given element's key and data in a key heap
  @param  h    pointer to the heap structure
  @param  x    the element whos key and data to update
  @param  key  the new key
  @param  data the new data
  @return      the element's original data
*/
void *
tmap_fibheap_replacekeydata(tmap_fibheap_t *h, tmap_fibheap_element_t *x, int key, void *data);
/*! 
  changes a given element's key and data in a void heap
  @param  fnct  the comparison function of type tmap_fibheap_voidcmp
  @return       pointer to the initialized key heap
*/
tmap_fibheap_t *
tmap_fibheap_makeheap(tmap_fibheap_voidcmp fnct);
/*! 
  changes the void heap's comparison function
  @param  h     pointer to the heap structure
  @param  fnct  the comparison function of type tmap_fibheap_voidcmp
  @return       the previous comparison function
*/
tmap_fibheap_voidcmp 
tmap_fibheap_setcmp(tmap_fibheap_t *h, tmap_fibheap_voidcmp fnct);
/*! 
  changes the void heap's negative infinity data
  @param  h     pointer to the heap structure
  @param  data  the data to represent negative infinity 
  @return       the previous negative infinity data
*/
void *
tmap_fibheap_setneginf(tmap_fibheap_t *h, void *data);
/*! 
  insert the given data into a void heap
  @param  h     pointer to the heap structure
  @param  data  the data to insert
  @return       pointer to the element inserted 
*/
tmap_fibheap_element_t *
tmap_fibheap_insert(tmap_fibheap_t *h, void *data);
/*! 
  get the data on top of the heap
  @param  h  pointer to the heap structure
  @return    the minimum data
  @details   removes the minimum element from the heap
*/
void *
tmap_fibheap_extractmin(tmap_fibheap_t *h);
/*! 
  get the data on top of the heap
  @param  h    pointer to the heap structure
  @return      the minimum data
  @details     does not remove the minimum element from the heap
*/
void *
tmap_fibheap_min(tmap_fibheap_t *h);
/*! 
  replaces the given element's data in the heap
  @param  h    pointer to the heap structure
  @param  x    the element whos data to update
  @param  data the new data
  @return      the element's original data
*/
void *
tmap_fibheap_replacedata(tmap_fibheap_t *h, tmap_fibheap_element_t *x, void *data);
/*! 
  delete the given element from the heap
  @param  h    pointer to the heap structure
  @param  x    the element whos data to delete 
  @return      the deleted element's data
*/
void *
tmap_fibheap_delete(tmap_fibheap_t *h, tmap_fibheap_element_t *x);
/*! 
  delete the given the heap
  @param  h    pointer to the heap structure
  @details     does not destroy the data within
*/
void 
tmap_fibheap_deleteheap(tmap_fibheap_t *h);
/*! 
  merges two heaps
  @param  ha  pointer to the heap structure #1
  @param  hb  pointer to the heap structure #2
  @return     pointer to the merged heap
  @details    ha is used to store the merged heap
*/
tmap_fibheap_t *
tmap_fibheap_union(tmap_fibheap_t *ha, tmap_fibheap_t *hb);
#endif 
