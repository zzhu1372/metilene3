#ifndef SORT_H
#define SORT_H

/*
 * sort.h
 * declarations for various sorting algorithms
 *
 * @author Steve Hoffmann
 * @date Mon 27 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: sort.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/sort.h $
 */
 #include "basic-types.h"


 Uint cmp_dbl(Uint, Uint, void *, void *);
 Uint *quickSort(void *, void *, Uint, Uint (*cmp)(Uint, Uint, void *, void *), void *);
#endif

