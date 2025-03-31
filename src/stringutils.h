#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include "basic-types.h"

#ifndef ALLOCMEMORY
	#include "memory.h"
#endif

#define TAB '\t' /*0x09*/
#define LF  '\n' /*0x0A*/
#define VT  '\v'
#define FF  '\f'
#define CR  '\r' /*0x0D*/
#define SP  ' '
#define DQT '\"'
#define SQT '\''


#define COPYSTR(M,D,S,L) 		D=ALLOCMEMORY(M,D,char,L+1); \
								strncpy(D,S,L);\
								D[L]='\0'

#define INDENT(s,x,c) 			{int p; for(p=0;p<x;p++) fprintf(s,"%c",c);}
#define ISWHITESPACE(C) 		(C == SP || C == TAB ||  C == LF \
								|| C == VT || C == CR || C == FF ) 
#define ISQUOTE(C) 				(C== DQT || C==SQT)
#define SETSTR(X,I) 			(X)->strings[(I)].str
#define SETSTRLEN(X,I)			(X)->strings[(I)].len

#define APPENDCHAR(S, X, L, Y) 	(X)=ALLOCMEMORY(S, X, char, L+1);\
										 X[L-1]=Y;\
										 X[L]=0 

typedef struct{

	char* str;
	Uint len;

} string_t;

typedef struct{

	string_t* strings;
	Uint noofstrings;

} stringset_t;

char * strtok_bl (char *, char *, char **);
stringset_t *tokensToStringset(void *, char *, char *, Uint);
void destructStringset(void *, stringset_t *);
char *my_strdup(const char *str);

#endif

