#ifndef FILEIO_H
#define FILEIO_H

/*
 * fileio.h
 * declarations for file io
 *
 * @author Steve Hoffmann
 * @date Sat 25 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: fileio.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/fileio.h $
 */

#ifndef ALLOCMEMORY
	#include "memory.h"
#endif
#include <math.h>
#include "stringutils.h"

typedef struct {
  char *filename;
  FILE *fp;
  char eof;
  char columndelim;
  char linedelim;
  size_t buffersize;
  size_t bufferfill;
  size_t ptr;
  char *buffer;
} fileiterator_t;

fileiterator_t* initFileIterator (void *space, char *filename);
void closeFileIterator (void *space, fileiterator_t *fi);
Uint readcsvlines(void *space, fileiterator_t *fi, char delim, Uint linecount, stringset_t ***out); 

#endif
