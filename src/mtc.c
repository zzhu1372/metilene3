/*
 *
 *  mtc.c
 *  
 * 
 *  @author Christian Otto
 *  @company ecSeq Bioinformatics
 *  @date 23/06/2018 12:58:01 PM CEST
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include "basic-types.h"
#include "info.h"
#include "stringutils.h"
#include "mathematics.h"
#include "sort.h"
#include "metseg.h" // to support struct list_out
#include "mtc.h"

/*----------------------- multiple_testing_correction --------------------------
 *    
 * @brief performs multiple testing correction to p-values
 * @author Christian Otto
 *   
 */
void multiple_testing_correction(list_out *list, int run_mode, int mtc_method){
  double *p, *q;
  int i;

  /* extract p-values */
  p = ALLOCMEMORY(NULL, NULL, double, list->i);
  for(i=0; i<list->i; i++){
    if (run_mode == 1 || run_mode == 2) {
      p[i] = list->segment_out[i].p;
    }
    else if (run_mode == 3) {
      p[i] = list->segment_out[i].mwu;
    }
  }

  /* perform correction */
  if (mtc_method == BONFERRONI){
    q = bonferroni_correction(p, list->i, list->numberTests);
  }
  else if (mtc_method == BENJAMINI_HOCHBERG){
    q = benjamini_hochberg_correction(p, list->i, list->numberTests);
  }
  else {
    fprintf(stderr, "Invalid method for multiple testing correction. Exit forced.\n");
    exit(-1);
  }
  
  /* store q-values */
  for(i=0; i<list->i; i++){
    list->segment_out[i].q = q[i];
  }

  /* clean-up */
  FREEMEMORY(NULL, p);
  FREEMEMORY(NULL, q);
}

/*-------------------------- bonferroni_correction -----------------------------
 *    
 * @brief perform Bonferroni correction, controlling the family-wise error rate (FWER)
 * @author Christian Otto
 *   
 */
double * bonferroni_correction(double *p, int n, int num_tests){

  double *q = ALLOCMEMORY(NULL, NULL, double, n);
  for (int i=0; i<n; i++){
    q[i] = p[i] * ((double) num_tests);
    if (q[i] > 1){
      q[i] = 1.0;
    }
  }

  return q;
}


/*----------------------- benjamini_hochberg_correction ------------------------
 *    
 * @brief perform Benjamini-Hochberg correction, controlling the false discovery rate (FDR)
 * @author Christian Otto
 *   
 */
double * benjamini_hochberg_correction(double *p, int n, int num_tests){
  int i;
  double *q;
  Uint *sorted;
  
  q = ALLOCMEMORY(NULL, NULL, double, n);
  sorted = quickSort(NULL, p, n, cmp_dbl, NULL);

  for (i=0; i<n; i++){
    q[sorted[i]] = DMIN(1.0, (double)num_tests/((double)i+1) * p[sorted[i]]);
  }
  for (i=1; i<n; i++){
    if (q[sorted[n-i-1]] > q[sorted[n-i]]){
      q[sorted[n-i-1]] = q[sorted[n-i]];
    }
  }
  /* for debugging
  for (i=0; i<n; i++){
    printf("num_tests=%i, i=%d, sorted[i]=%d, p[i]=%g, p[sorted[i]]=%g, factor=%g, q[sorted[i]]=%g\n", num_tests, i, sorted[i], p[i], p[sorted[i]], q[sorted[i]]);
  }
  for (i=0; i<n; i++){
    //printf("%g\n", q[i]);
    printf("i=%d\tp=%lf\tq=%lf\n", i, p[i], q[i]);
  }
  exit(-1);
  */
  FREEMEMORY(NULL, sorted);
  return q;
}
