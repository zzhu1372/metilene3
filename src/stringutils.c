/*
 * stringutils.c
 * functions to manipulate strings
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: stringutils.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/stringutils.c $
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include <assert.h>
 #include "stringutils.h"
 #include "basic-types.h"


char* strtok_bl(char *s, char *delim, char **saveptr) {

  char *ret;

  /* init */
  if (s == NULL){
    s = *saveptr;
    if (s == NULL){
      *saveptr = NULL;
      return NULL;
    }
  }
  /* skip delims at begin */
  while(*s && strchr(delim, *s)){
    s++;
  }

  if (*s == 0){
    *saveptr = NULL;
    return NULL;
  }

  /* locate next delim or end */
  ret = s;
  while(*s && !strchr(delim, *s)){
    s++;
  }

  if (*s == 0){
    *saveptr = NULL;
  }
  else {
    *s = 0;
    s++;
    *saveptr = s;
  }  

  return ret;
}


stringset_t* tokensToStringset(void *space, char* delim, char* toTokens, 
                                 Uint len){
	Uint toklen;
	char* token;
        char* saveptr;
	char* buffer;
	stringset_t *set;
	
	set = ALLOCMEMORY(space, NULL, stringset_t, 1);
	set->noofstrings = 0;
	set->strings = NULL;
									
	if (toTokens == NULL || len == 0)
	return set;
													
	buffer = ALLOCMEMORY(space, NULL, char, len+1);
	if (buffer == NULL) {
		fprintf(stderr, "copy tokenstring %s to buffer failed.\n", toTokens);
		exit(-1);	
	}
	buffer = memcpy(buffer, toTokens, len+1);
        buffer[len] = 0;

	
	token = strtok_bl(buffer, delim, &saveptr);
	
	while(token != NULL) {
		
			toklen = strlen(token);
			set->noofstrings++;
			set->strings = ALLOCMEMORY(space, set->strings, string_t, set->noofstrings);
			set->strings[set->noofstrings-1].str = ALLOCMEMORY(space, NULL, char, toklen+1);
			set->strings[set->noofstrings-1].str = memcpy(set->strings[set->noofstrings-1].str, token, toklen);
			set->strings[set->noofstrings-1].str[toklen]='\0'; 
			set->strings[set->noofstrings-1].len = toklen;
			token = strtok_bl(NULL, delim, &saveptr);
	}

	FREEMEMORY(space, buffer);
	return set;
}

void destructStringset(void *space, stringset_t *s) {

  if (s->strings) {
    for(Uint i=0; i < s->noofstrings; i++) {
      if(s->strings[i].str != NULL)	
        FREEMEMORY(space, s->strings[i].str);
    }

    FREEMEMORY(space, s->strings);
  }

  FREEMEMORY(space, s);
}

char *my_strdup(const char *str) {
  size_t len = strlen(str);
  char *x = (char *)malloc(len+1); /* 1 for the null terminator */
  if(!x) return NULL; /* malloc could not allocate memory */
  memcpy(x,str,len); /* copy the string into the new buffer */
  x[len]=0;
  return x;
}
