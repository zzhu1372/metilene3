#ifndef MEMORY_H
#define MEMORY_H
/*
#include "memman.h"
#include "memmac.h"
*/
#include <assert.h>
#include <stdlib.h>

#define ALLOCMEMORY(X,PTR,TYPE,SIZE) bl_realloc((PTR),sizeof(TYPE)*(SIZE))
#define CALLOCMEMORY(X,TYPE,SIZE) bl_realloc((PTR),(SIZE),sizeof(TYPE))

#define FREEMEMORY(X,PTR) free(PTR); PTR=NULL

void* bl_realloc(void *, size_t);
void* bl_calloc(void *, size_t, size_t);

#endif

