/*
	  This file is part of gdub.
	  (C) 2006 Steve Hoffmann 
 
	  gdub is free software; you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published
	  by the Free Software Foundation; either version 2, or (at your
	  option) any later version.
 
	  gdub is distributed in the hope that it will be useful, but
	  WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	  General Public License for more details.
 
	  You should have received a copy of the GNU General Public License
	  along with gdub; see the file COPYING.  If not, write to the
	  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
	  Boston, MA 02111-1307, USA.	
 
 */

/**
 * stack.h
 * implementation of a simple stack for int
 *
 * @author Steve Hoffmann
 * @email steve@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Tue Oct 28 10:42:34 CET 2008
 */

/*
 *  SVN
 *  Revision of last commit: $Rev: 73 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-10-29 10:03:28 +0100 (Wed, 29 Oct 2008) $
 *  Id: $Id: stack.h 73 2008-10-29 09:03:28Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/stack.h $
 */

#ifndef STACK_H
#define STACK_H

#include <stdio.h>
#include <stdlib.h>
#include "basic-types.h"
#include "metseg.h"

#define SEGMENTSTACK_NULL_TYPE NULL
#define SEGMENTSTACKINC 10000
#ifndef SEGMENTSTACKBASEINC
#define SEGMENTSTACKBASEINC SEGMENTSTACKINC
#endif

typedef segment_p_t Segmentstackelement;

typedef struct{
	Segmentstackelement* stackspace;
	Lint allocelem;
	Lint top;
} Segmentstack;

void bl_segmentstackInit(Segmentstack* stack, Lint allocelem);
void bl_segmentstackDestruct(Segmentstack *stack);
BOOL bl_segmentstackIsEmpty(Segmentstack *stack);
void bl_segmentstackPush(Segmentstack *stack, Segmentstackelement *elem);
Segmentstackelement* bl_segmentstackPop(Segmentstack *stack);
Segmentstackelement* bl_segmentstackTop(Segmentstack *stack);
Segmentstackelement* bl_segmentstackTopN(Segmentstack *stack, Lint n);
Lint bl_segmentstackSize(Segmentstack *stack);

#endif
