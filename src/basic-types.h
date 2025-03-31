#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#include <stdint.h>

#ifdef __CYGWIN__ 

#define CRLF '\r'

#else
#define CRLF ' '
#endif

#define MAXBUFFERSIZE 10000
#define BASEINC 10000
#define MAX_INT_LENGTH 50
typedef unsigned char Uchar;
typedef unsigned int Uint;
typedef signed long long Lint;
typedef signed long long int LLint;
typedef signed int Sint;
typedef unsigned char BOOL;

#define True 1
#define False 0

#ifndef TRUE
#define TRUE True
#endif

#ifndef FALSE
#define FALSE False
#endif

typedef struct {
  int  a, 
       b;
} PairSint; 

typedef struct{
  Lint a,
       b;
} PairLSint;

typedef struct {
  Uint  a, 
       b;
} PairUint; 


typedef struct {
  Uint a,
      b,
      c;
} TripleUint;


typedef struct {
  int a,
      b,
      c;
} TripleSint;

typedef struct {
  int a,
      b,
      c,
      d;
} QuadSint;




#endif

