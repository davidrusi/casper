/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2011 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: hash.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.12 $      $Date: 2010/12/16 04:08:55 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   A simple hash table implementation for strings, contributed by John Stone,
 *   derived from his ray tracer code.
 ***************************************************************************/
#ifndef HASH_H
#define HASH_H

typedef struct hash_node_t {
  int data;                           /* data in hash node */
  char * key;                   /* key for hash lookup */
  struct hash_node_t *next;           /* next node in hash chain */
} hash_node_t;

typedef struct hash_t {
  struct hash_node_t **bucket;        /* array of hash nodes */
  int size;                           /* size of the array */
  int entries;                        /* number of entries in table */
  int downshift;                      /* shift cound, used in hash function */
  int mask;                           /* used to select bits for hashing */
} hash_t;


#define HASH_FAIL -1

#if defined(VMDPLUGIN_STATIC)
#define VMDEXTERNSTATIC static
#include "hash.c"
#else

#define VMDEXTERNSTATIC 

#ifdef __cplusplus
extern "C" {
#endif

void hash_init(hash_t *tptr, int);

int hash_lookup (const hash_t *tptr, const char *);

int hash_insert (hash_t *tptr, const char *, int);

int hash_update(hash_t *tptr, const char *key, int data);
char* m_strdup(const char *o);
int hash_delete (hash_t *tptr, const char *);

void hash_destroy(hash_t *tptr);

char *hash_stats (hash_t *tptr);

//int hash(const hash_t *tptr, const char *key);
#ifdef __cplusplus
}
#endif

#endif

#endif
