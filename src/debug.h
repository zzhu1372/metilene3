 #ifndef DEBUG_H
 #define DEBUG_H

/*
 *
 *	debug.h
 *  debug messages
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 08/26/2007 07:17:44 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: debug.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/debug.h $
 */

 #include <stdarg.h>
 #include <stdio.h>
 #include <string.h>

#ifndef DBGLEVEL
#define DBGLEVEL 0
#endif

#ifndef DBGDEVICE
#define DBGDEVICE stderr
#endif

#define DBGL(L, X, ... ) debuglevel (__FILE__, __LINE__, L, X, __VA_ARGS__)
#define DBG(X, ...) debugmsg(__FILE__, __LINE__, X, __VA_ARGS__)
#define DBGEXIT(X, ...) { debugmsg (__FILE__, __LINE__, X, __VA_ARGS__); \
                    exit(-1); }


/*deprecated*/
#define DEBUG(X, ...) debugmsg(__FILE__, __LINE__, X, __VA_ARGS__)

int debugmsg(const char *, const int, const char *fmt, ...);
int debuglevel(const char *, const int, int, const char *fmt, ...);

#endif
