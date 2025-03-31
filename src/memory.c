
/*
 *  memory.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 02.01.2010 00:11:46 CET
 *  
 */

#include "memory.h"
#include <stdlib.h>
#include <assert.h>

void*
bl_realloc(void *ptr, size_t sz) {
  ptr = realloc(ptr, sz);
  assert(ptr != NULL);
  return ptr;
}

void*
bl_calloc(void *ptr, size_t nelem, size_t sz) {
  ptr = calloc(nelem, sz);
  assert(ptr != NULL);
  return ptr;
}
