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
 * stack.c
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
 *  Id: $Id: stack.c 73 2008-10-29 09:03:28Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/stack.c $
 */

#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include "segmentstack.h"

/*----------------------------- bl_stackInit -----------------------------------
 *    
 * @brief 	init stack
 * @author 	Steve Hoffmann
 *   
 */
void bl_segmentstackInit(Segmentstack *stack, Lint allocelem) {
  if (allocelem <= 0){
    DBG("stack.c: Attempt to initialize a stack of size %d. Exit forced.\n",
	allocelem);
    exit(-1);
  }
  stack->stackspace = (Segmentstackelement *) malloc(sizeof(Segmentstackelement) * allocelem);
  if (stack->stackspace == NULL){
    DBG("stack.c: Memory allocation failed. Exit forced.\n", NULL);
    exit(-1);
  }
  stack->allocelem=allocelem;
  stack->top=-1;

}

/*--------------------------- bl_stackDestruct ---------------------------------
 *    
 * @brief 	destruct stack
 * @author 	Steve Hoffmann
 *   
 */
void bl_segmentstackDestruct(Segmentstack *stack) {
  free(stack->stackspace);
  stack->top = 0;
  stack->allocelem = 0;
}

/*---------------------------- bl_stackIsEmpty ---------------------------------
 *    
 * @brief 	returns if the stack is empty
 * @author 	Steve Hoffmann
 *   
 */
BOOL bl_segmentstackIsEmpty(Segmentstack *stack) {
  return (stack->top < 0);
}


/*----------------------------- bl_stackPush -----------------------------------
 *    
 * @brief 	pushs elements on the top of the stack
 * @author 	Steve Hoffmann
 *   
 */		
void bl_segmentstackPush(Segmentstack *stack, Segmentstackelement* elem) {
		
  if(stack->top >= stack->allocelem - 1) {
			
    stack->stackspace = (Segmentstackelement *) realloc(stack->stackspace, 
						 sizeof(Segmentstackelement) * 
						 (stack->allocelem + SEGMENTSTACKBASEINC)); 
    if (stack->stackspace == NULL || SEGMENTSTACKBASEINC <= 0){
      DBG("stack.c: Memory reallocation failed. Exit forced.\n", NULL);
      exit(-1);
    }
    stack->allocelem += SEGMENTSTACKBASEINC; 
  }

  stack->top++;
  memmove(&stack->stackspace[stack->top], elem, sizeof(Segmentstackelement));
}

/*------------------------------ bl_stackPop -----------------------------------
 *    
 * @brief 	pops the top of the stack
 * @author 	Steve Hoffmann
 *   
 */
Segmentstackelement* bl_segmentstackPop(Segmentstack *stack){
  if(!bl_segmentstackIsEmpty(stack)) {
    return &stack->stackspace[stack->top--];	
  }	
  return SEGMENTSTACK_NULL_TYPE;
}

/*------------------------------ bl_stackTop -----------------------------------
 *    
 * @brief 	returns top of the stack
 * @author 	Steve Hoffmann
 *   
 */
Segmentstackelement* bl_segmentstackTop(Segmentstack *stack){	
  if(!bl_segmentstackIsEmpty(stack)) {
    return &stack->stackspace[stack->top];
  }
  return SEGMENTSTACK_NULL_TYPE;
}
 
/*------------------------------ bl_stackTopN --------------------------------
 *    
 * @brief 	returns Nth highest object of the stack
 *              with N = 0,..,numofelems - 1
 * @author 	Steve Hoffmann
 *   
 */
Segmentstackelement* bl_segmentstackTopN(Segmentstack *stack, Lint n){	
  if(!bl_segmentstackIsEmpty(stack) && n >= 0 && n <= stack->top) {
    return &stack->stackspace[stack->top - n];
  }
  return SEGMENTSTACK_NULL_TYPE;
}

/*------------------------------ bl_stackSize ----------------------------------
 *    
 * @brief 	returns number of elements on the stack
 * @author 	Steve Hoffmann
 *   
 */
Lint bl_segmentstackSize(Segmentstack *stack) {
  return (stack->top + 1);    
}

