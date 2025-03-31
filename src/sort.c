
/*
 * sort.c
 * implementation of various sorting algorithms
 *
 * @author Steve Hoffmann
 * @date Mon 27 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 77 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-17 13:16:59 +0100 (Mon, 17 Nov 2008) $
 *
 *  Id: $Id: sort.c 77 2008-11-17 12:16:59Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sort.c $
 */

 #include <math.h>
 #include <string.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "vstack.h"
 #include "mathematics.h"
 #include "sort.h"
 #include "debug.h"

 Uint cmp_dbl(Uint a, Uint x, void *data, void *info) {
	double *d = (double*) data;

	/*if(floor(d[a]) > floor(d[x])) return 1;
	if(floor(d[a]) < floor(d[x])) return 2;
*/
	if (fabs((double) d[a] - d[x]) <= FLT_EPSILON * FLT_MIN) return 0;
	if ((double) d[a] - d[x] >  FLT_EPSILON * FLT_MIN) return 1;
	if ((double) d[x] - d[a] >  FLT_EPSILON * FLT_MIN) return 2;

	return 0;
 }

Uint *quickSort(void *space, void* toSort, Uint size, 
		Uint (*cmp)(Uint, Uint, void *, void*),
		void *info) {
  int left, left2, right, right2;
  PairSint ins, *lr;
  Uint i, resc, *sorted, x;
  VStack vstack;
	
  sorted = ALLOCMEMORY(space, NULL, Uint, size);
  for (i=0; i < size; i++) sorted[i]=i;
  ins.a = 0;
  ins.b = size-1;
  bl_vstackInit(&vstack, 10000, sizeof(PairSint));
  bl_vstackPush(&vstack, &ins);
   
  while (!bl_vstackIsEmpty(&vstack)){
    lr = (PairSint *) bl_vstackPop(&vstack, NULL);
    left = lr->a;
    right = lr->b;
    free(lr);
    while (left < right) {
      x=sorted[(left+right)/2];
      left2  = left;
      right2 = right;
	
      do {
	while(cmp(sorted[left2],  x, toSort, info)==2){	
	  left2++;
	}
	while(cmp(sorted[right2], x, toSort, info)==1){ 
	  right2--;
	}
			
	if(left2 <= right2) {
	  resc = sorted[right2];
	  sorted[right2]=sorted[left2];
	  sorted[left2]=resc;
	  left2++;
	  right2--;
	} 	
      } while (right2 >= left2);
			

      if ((left2-left) > (right-left2))  {		
	ins.a = left;
	ins.b = right2;
	bl_vstackPush(&vstack, &ins);
	left  = left2;
      } else {
	ins.a = left2;
	ins.b = right;
	bl_vstackPush(&vstack, &ins);
	right = right2;
      }
    }
  }
  bl_vstackDestruct(&vstack, NULL);
  return sorted;
}
