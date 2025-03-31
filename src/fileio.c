/*
 * fileio.c
 * functions to manipulate and read files
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: fileio.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/fileio.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "stringutils.h"
#include "basic-types.h"
#include "fileio.h"

#ifndef DIR_SEPARATOR
#define DIR_SEPARATOR '/'
#endif

#if defined (_WIN32) || defined (__MSDOS__) || defined (__DJGPP__) || \
  defined (__OS2__)
#define HAVE_DOS_BASED_FILE_SYSTEM
#ifndef DIR_SEPARATOR_2 
#define DIR_SEPARATOR_2 '\\'
#endif
#endif

/* Define IS_DIR_SEPARATOR.  */
#ifndef DIR_SEPARATOR_2
# define IS_DIR_SEPARATOR(ch) ((ch) == DIR_SEPARATOR)
#else /* DIR_SEPARATOR_2 */
# define IS_DIR_SEPARATOR(ch) \
  (((ch) == DIR_SEPARATOR) || ((ch) == DIR_SEPARATOR_2))
#endif /* DIR_SEPARATOR_2 */


int
bl_fgetsBuffered(void *space, fileiterator_t *fi, char **str) {
  char *line = NULL;
  Uint linebuffersize = 1000;
  Uint len = 0;


  line = ALLOCMEMORY(space, NULL, char, linebuffersize);
  
  //if the iterator is called for the first time with buffers
  if(!fi->buffer) { 
    fi->buffer = ALLOCMEMORY(space, NULL, char, fi->buffersize+1);
    fi->ptr = 0;
    fi->bufferfill =0;
//    fprintf(stderr, "initializing empty buffer\n");
  } else {
 //   fprintf(stderr, "bufferfill is at: %zd, ptr is at %zd\n", fi->bufferfill, fi->ptr);
  }

  //reload if pointer points to last character in buffer 
  if(fi->ptr == fi->bufferfill) {
    fi->bufferfill = fread(fi->buffer, 1, fi->buffersize, fi->fp);
 //   fprintf(stderr, "refilling buffer prior to loop.\n");
    fi->ptr =0;
  }
  
 // fprintf(stderr, "prior to loop current char is '%c'\n", fi->buffer[fi->ptr]);

  //as long as the buffer is not empty and there was no line break
  while(fi->ptr < fi->bufferfill && fi->buffer[fi->ptr] != '\n') {

    if(len+1 == linebuffersize) {
      linebuffersize += linebuffersize;
      line = ALLOCMEMORY(space, line, char, linebuffersize+1);
   //   fprintf(stderr, "reallocating linebuffer\n");
    }

    line[len] = (char) fi->buffer[fi->ptr];
    len++;
    fi->ptr++;

    //check whether it is necessary to reload the buffer
    if(fi->ptr == fi->bufferfill) {
      fi->bufferfill = fread(fi->buffer, 1, fi->buffersize, fi->fp);
      fi->ptr =0;
   //   fprintf(stderr, "refilling buffer inside loop.\n");
    } 
  }

  if(fi->ptr == fi->bufferfill) { 
    FREEMEMORY(NULL, line);
    return EOF;
  }
  
  assert(fi->buffer[fi->ptr]=='\n');
  fi->ptr++;

  line[len] = '\0';
  *str = line;

  return len;
}

/*----------------------------- initfileiterator -----------------------------
 *    
 * @brief initalize and open fileiterator
 * @author Steve Hoffmann 
 *   
 */

  fileiterator_t*
initFileIterator (void *space, char *filename)
{

  fileiterator_t *fi;   
  fi = ALLOCMEMORY(space, NULL, fileiterator_t, 1);
  fi->filename = filename;
  fi->fp = fopen(filename, "r");
  if (fi->fp == NULL){
    fprintf(stderr, "Opening of file %s failed. Exit forced.\n", filename);
    exit(EXIT_FAILURE);
  }

  fi->buffersize = 16384;
  fi->buffer = NULL;
  fi->ptr = 0;
  fi->bufferfill = 0;
  fi->eof = 0;
  return fi;
}


/*--------------------------- closeFileIterator ---------------------------
 *    
 * @brief destruct file iterator
 * @author Steve Hoffmann 
 *   
 */
 
void
closeFileIterator (void *space, fileiterator_t *fi)
{

  if(fi->fp) {
    fclose(fi->fp);
    fi->fp = NULL;
  }

  if(fi->buffer)
    FREEMEMORY(NULL, fi->buffer);

  return;
}


Uint
readcsvlines(void *space, 
    fileiterator_t *fi, char delim, 
    Uint linecount, stringset_t ***out) {

  Uint i=0, k=0, len;
  char *line; 
  stringset_t **csv=NULL;
    
  csv = ALLOCMEMORY(space, NULL, stringset_t*, linecount);
  for(k=0; k < linecount; k++) {
    csv[k] = NULL;
  }
  
  while(i < linecount && (len=bl_fgetsBuffered(space, fi, &line)) != EOF) {
    //fprintf(stderr, "len:%d -> '%s'\n", len, line);
    csv[i] = tokensToStringset(space, "\t", line, len);  
    FREEMEMORY(space, line);
    i++;
  }


  *out = csv;
  return i;
}
