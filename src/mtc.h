#ifndef MTC_H
#define MTC_H
/*
 *
 *  mtc.h
 *  
 * 
 *  @author Christian Otto
 *  @company ecSeq Bioinformatics
 *  @date 23/06/2018 12:58:01 PM CEST
 *
 */
#include "metseg.h"

#define BONFERRONI 1
#define BENJAMINI_HOCHBERG 2

void multiple_testing_correction(list_out *list, int runMode, int mtcMethod);
double* bonferroni_correction(double *p, int n, int numTests);
double* benjamini_hochberg_correction(double *p, int n, int numTests);

# endif
