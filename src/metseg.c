/*
 *  metseg.c
 *  
 *
 *  @author Frank Juehling and Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 09/09/2014 08:54:52 AM CEST
 *  
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include "basic-types.h"
#include "manopt.h"
#include "info.h"
#include "stringutils.h"
#include "fileio.h"
#include "metseg.h"
#include "mtc.h"
#include <string.h>
#include "mathematics.h"
#include "vstack.h"
#include "segmentstack.h"
#include <time.h>
#include <ctype.h>
#include <float.h>

#define MAXN 13
#define MAXM 13

char *version = "0.2-8";
unsigned char mute=0;
pthread_mutex_t updatemtx;
double get_ratio(double *a, int m, double *b, int n);
/*------------------------------- initSegment --------------------------------
 *    
 * @brief initialize new segment
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void
initSegment(segment_t *seg) {
  seg->n =0;
  seg->prob=-1;
  seg->test=-1;
  seg->length = 0;
//  seg->chr="";
  seg->init=0;
  seg->value = NULL;
  seg->pos = NULL;
  seg->next = NULL;
  seg->parent = NULL;
}

void
destructSegment(segment_t *seg) {
    for(int i=0;i<seg->n;i++) {
         FREEMEMORY(NULL, seg->chr);
         FREEMEMORY(NULL, seg->value[i]);
    }
  FREEMEMORY(NULL, seg->chr);
  FREEMEMORY(NULL, seg->pos);
  FREEMEMORY(NULL, seg->value);
  FREEMEMORY(NULL, seg->next);
  FREEMEMORY(NULL, seg->parent);
  FREEMEMORY(NULL, seg);
  
  
}

void
destructCpg(cpg_t *cpg) {
    
  FREEMEMORY(NULL, cpg->groupA);
  FREEMEMORY(NULL, cpg->groupB);
  FREEMEMORY(NULL, cpg->chr);
  FREEMEMORY(NULL, cpg);
}

void
destructSegmentSet(segmentset_t *set) {
    segment_t *tmp = set->head;
    segment_t *next;
    while (set->seg) {
        next = tmp->next;
        destructSegment(tmp);
        tmp = next;
    }
  FREEMEMORY(NULL, set->chr);
  FREEMEMORY(NULL, set->nextchr);
  FREEMEMORY(NULL, set);
}


void
initSegmentSet(segmentset_t *set) {
  set->n=0;
  set->firststop=-1;
  set->nextchr=NULL;
  set->nextstart=-1;
  set->seg=NULL;
  set->head=NULL;
  set->tail=NULL;
  set->chr=NULL;
}

segment_t 
*addNewSegmentToSet(segmentset_t *set) {
  segment_t *seg = NULL;
  seg = ALLOCMEMORY(NULL, NULL, segment_t, 1);
  initSegment(seg);
 
  seg->parent= set->tail;
  if(set->head != NULL) {
      set->tail->next = seg;
  }
  else {
      set->head=seg;
  }
  set->tail=seg;   
  set->n++;
  return seg;
}


segment_t 
*getSegment(segmentset_t *set, int n) {
    segment_t *s = NULL;
    if(n<0 || n>=set->n) return s;
    s = set->head; 
    for(int i=1;i<=n;i++) {
        s = s->next;
    }
    return s;
}


void
removeThisSegmentFromSet(segmentset_t *set,segment_t *seg) {
    if(seg->parent == NULL) {
        if(set->n==1) {
            set->head=NULL;
            set->tail=NULL;
            set->n--;

            seg->next=NULL;
            seg->parent=NULL;
            return;
        }
        else {
            set->head=set->head->next;
            set->head->parent=NULL;
            set->n--;

            seg->next=NULL;
            seg->parent=NULL;
            return;
        }        
    }
//remove tail        
    if(seg->next == NULL) {
        if(set->n==1) {
            set->head=NULL;
            set->tail=NULL;
            set->n--;

            seg->next=NULL;
            seg->parent=NULL;
            return;
        }
        else {
            set->tail = set->tail->parent;
            set->tail->next=NULL;
            set->n--;

            seg->next=NULL;
            seg->parent=NULL;
            return;
        }
    }
//remove an object within the list
    seg->parent->next=seg->next;
    seg->next->parent=seg->parent;
    set->n--;

    seg->next=NULL;
    seg->parent=NULL;
    return;
}
    

void
removeSegmentFromSet(segmentset_t *set, int n) {
    //no object to remove
    if(set->n == 0 || n>= set->n) return;
    
    //remove head
    if(n == 0) {
        if(set->n==1) {
            segment_t *seg = set->head;
            set->head=NULL;
            set->tail=NULL;
            set->n--;
            
            seg->next=NULL;
            seg->parent=NULL;
            return;
        }
        else {
            segment_t *seg = set->head;
            set->head=set->head->next;
            set->head->parent=NULL;
            set->n--;
            
            seg->next=NULL;
            seg->parent=NULL;
            return;
        }        
    }
    else 
//remove tail        
        if(n == set->n-1) {
            if(set->n==1) {
                segment_t *seg = set->tail;
                set->head=NULL;
                set->tail=NULL;
                set->n--;
            
                seg->next=NULL;
                seg->parent=NULL;
                return;
            }
            else {
                segment_t *seg = set->tail;
                set->tail = set->tail->parent;
                set->tail->next=NULL;
                set->n--;
             
                seg->next=NULL;
                seg->parent=NULL;
                return;
           }
        }
    //remove an object within the list
        else {
            segment_t *seg = set->head; 
            for(int i=1;i<=n;i++) {
                seg = seg->next;
            }
            seg->parent->next=seg->next;
            seg->next->parent=seg->parent;
            set->n--;
            
            seg->next=NULL;
            seg->parent=NULL;
            return;
        }
}





/*-------------------------------- setSegment --------------------------------
 *    
 * @brief set values in a segment
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void 
setSegment(segment_t *seg, char *chr, int start, int stop, 
    double prob, double meandiff, double test, double sigcp) {

  seg->chr = chr;
  seg->start = start;
  seg->stop = stop;
  seg->prob = prob;
  seg->test = test;
  seg->meandiff = meandiff;
  seg->sigcp = sigcp;

  // assert((stop-start+1>=10)||(prob>1));
  
}

/*--------------------------------- get_meandiff ---------------------------------
 *    
 * @brief compute ratio of two group means
 * @author Frank Juehling
 *   
 */

double 
get_meandiff(cpg_t *cpg,double *a, int m, double *b, int n) {
  int i;
  double mean1, mean2, meandiff;
  mean1=0;
  mean2=0;
  for(i=0; i<m; i++) {
    mean1+=a[i];
  }
  mean1/=(double) m;
  for(i=0; i<n; i++) {
    mean2+=b[i];
  }
  mean2/=(double) n;
  meandiff=mean1-mean2;
//  fprintf(stdout,"meandiff %f\n",meandiff);
  cpg->methA=mean1;
  cpg->methB=mean2;
  return meandiff;
}

/*--------------------------------- calcMax ----------------------------------
 *    
 * @brief get the maximum quadrant
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

double 
calcMax(double *l1, double *l2, double c1, double c2, int a, int b, 
    double m, double n) {

  double d = 0;

  d=MAX(fabs(((c1)/m)-((c2)/n)),fabs(((c1+l1[0])/m)-((c2+l2[0])/n)));
  d=MAX(fabs(((c1+l1[a])/m)-((c2+l2[a])/n)),d);
  d=MAX(fabs(((c1+l1[b])/m)-((c2+l2[b])/n)),d);
  d=MAX(fabs(((c1+l1[0]+l1[a])/m)-((c2+l2[0]+l2[a])/n)),d);
  d=MAX(fabs(((c1+l1[0]+l1[b])/m)-((c2+l2[0]+l2[b])/n)),d);
  d=MAX(fabs(((c1+l1[b]+l1[a])/m)-((c2+l2[b]+l2[a])/n)),d);
  d=MAX(fabs(((c1+l1[b]+l1[a]+l1[0])/m)-((c2+l2[b]+l2[a]+l2[0])/n)),d);

  return d; 
}

/*--------------------------------- counter ----------------------------------
 *    
 * @brief count the data points in the four quadrants with center (x,y)
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

double 
counter(double x, double y, double *x0, double *x1, int m, 
    double *y0, double *y1, int n) {

  double **c, **l;
  int i;

  c = ALLOCMEMORY(NULL, NULL, double*, 2);
  l = ALLOCMEMORY(NULL, NULL, double*, 2);

  c[0] = ALLOCMEMORY(NULL, NULL, double, 4);
  c[1] = ALLOCMEMORY(NULL, NULL, double, 4);
  l[0] = ALLOCMEMORY(NULL, NULL, double, 5);
  l[1] = ALLOCMEMORY(NULL, NULL, double, 5);
  memset(c[0], 0, sizeof(double)*4);
  memset(c[1], 0, sizeof(double)*4);
  memset(l[0], 0, sizeof(double)*5);
  memset(l[1], 0, sizeof(double)*5);


  for(i=0; i < m; i++) {
    //central point
    if(x0[i] == x && x1[i] == y) {l[0][0]++; continue;}  
    //   if(lx[i] == x || ly[i] == y) {continue;}  
    //points on borders of quadrants
    if(x0[i] == x && x1[i] < y) {l[0][1]++; continue;}  
    if(x0[i] == x && x1[i] > y) {l[0][2]++; continue;}  
    if(x0[i] < x && x1[i] == y) {l[0][3]++; continue;}  
    if(x0[i] > x && x1[i] == y) {l[0][4]++; continue;}  

    if(x0[i] > x) {
      if(x1[i] > y)
        c[0][0]++;
      else
        c[0][1]++;
    }
    else {
      if(x1[i] > y)
        c[0][2]++;
      else
        c[0][3]++;
    }
  }

  for(i=0; i < n; i++) {
    //central point
    if(y0[i] == x && y1[i] == y) {l[1][0]++; continue;}  
    //  if(lx[i] == x || ly[i] == y) {continue;}  
    //points on borders of quadrants
    if(y0[i] == x && y1[i] < y) {l[1][1]++; continue;}  
    if(y0[i] == x && y1[i] > y) {l[1][2]++; continue;}  
    if(y0[i] < x && y1[i] == y) {l[1][3]++; continue;}  
    if(y0[i] > x && y1[i] == y) {l[1][4]++; continue;}  

    if(y0[i] > x) {
      if(y1[i] > y)
        c[1][0]++;
      else
        c[1][1]++;
    }
    else {
      if(y1[i] > y)
        c[1][2]++;
      else
        c[1][3]++;
    }
  }

  double d[4];

  d[0]=calcMax(l[0],l[1],c[0][0],c[1][0],2,4,m,n);
  d[1]=calcMax(l[0],l[1],c[0][1],c[1][1],1,4,m,n);
  d[2]=calcMax(l[0],l[1],c[0][2],c[1][2],2,3,m,n);
  d[3]=calcMax(l[0],l[1],c[0][3],c[1][3],1,3,m,n);

  //fprintf(stdout, "d[0]:%f, d[1]:%f, d[2]:%f, d[3]:%f, d[4]:%f\n", d[0], d[1], d[2], d[3], d[4]);

  double D = MAX(MAX(MAX(d[0],d[1]),d[2]),d[3]);

  FREEMEMORY(NULL, c[0]);
  FREEMEMORY(NULL, c[1]);
  FREEMEMORY(NULL, l[0]);
  FREEMEMORY(NULL, l[1]);
  FREEMEMORY(NULL, c);
  FREEMEMORY(NULL, l);

  return D;  
}

/*--------------------------------- kstest2d ---------------------------------
 *    
 * @brief two-dimensional ks-test
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void 
kstest2d(double *x0, double *x1, int m, double *y0, double *y1, int n, 
    double *kstest) {

  double cor[] = {0.0,0.0};
  double ks[] = {0.0,0.0};
  double d[] = {0.0,0.0};
  double dl1, dl2, s;

  for(int j=0; j < m; j++) {
    d[0] = MAX(d[0],counter(x0[j],x1[j],x0,x1,m,y0,y1,n));
  }

  for(int j=0; j< n; j++) {
    d[1]=  MAX(d[1],counter(y0[j],y1[j],x0,x1,m,y0,y1,n));
  }

  ks[1]=(d[0]+d[1])*0.5;

  cor[0] = rho(NULL, x0, x1, m);
  cor[1] = rho(NULL, y0, y1, n);

  dl1 = m;
  dl2 = n;
  s = sqrt(dl1*dl2/(dl1+dl2));
  kstest[0]  = kscdf(ks[1]*s/(1.0+sqrt(1.0-0.5*(cor[0]*cor[0]+cor[1]*cor[1]))*(0.25-0.75/s)));

  //fprintf(stdout, "s:%f, cor[0]:%f, cor[1]:%f, ks[1]:%f, test:%f\n", s, 
  //cor[0], cor[1], ks[1], kstest[0]);

  return;
}

/*---------------------------------- kstest ----------------------------------
 *    
 * @brief calculated the ks test
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void kstest(segment_t *seg , int a, int b, char mindiff, char mincpgs, char test, 
    double *ks, int *grpA, int noA, int *grpB, int noB, metseg_t* nfo){

  int i,j;
  int l = b-a+1;
  double **la;
  double **lb;
  double *x;
  double *y;
  double mean1=0;
  double mean2=0;
  double meandiff = 0;
  double dl1=noA;
  double dl2=noB;
  double dl = l;
  double p = 2;
  double faskstest[2] = {0,0};

  //mincpgs set to false by default so this condition is never fullfilled

  // mindiff = 0;
  // mincpgs = 0;
  // if(l< nfo->mincpgs && mincpgs) {
  //   ks[0]=2;ks[1]=0;ks[2]=2;
  // }
  //debug
  //  fprintf(stderr, "kstest for: %d\n",l);
  la = ALLOCMEMORY(NULL, NULL, double*, 2);
  lb = ALLOCMEMORY(NULL, NULL, double*, 2);

  la[0] = ALLOCMEMORY(NULL, NULL, double, l*noA);
  la[1] = ALLOCMEMORY(NULL, NULL, double, l*noA);
  lb[0] = ALLOCMEMORY(NULL, NULL, double, l*noB);
  lb[1] = ALLOCMEMORY(NULL, NULL, double, l*noB);
  memset(la[0], 0, sizeof(double)*l*noA);
  memset(la[1], 0, sizeof(double)*l*noA);
  memset(lb[0], 0, sizeof(double)*l*noB);
  memset(lb[1], 0, sizeof(double)*l*noB);

  x = ALLOCMEMORY(NULL, NULL, double, l);
  y = ALLOCMEMORY(NULL, NULL, double, l);

  ks[0]=0;ks[1]=0;ks[2]=0;

  for(i=a; i<=b ; i++) {
    for(j=0; j < noA; j++) {
      la[0][((i-a)*noA)+j]=seg->value[i][grpA[j]];
      la[1][((i-a)*noA)+j]=seg->pos[i];
      mean1+=seg->value[i][grpA[j]];
      x[i-a] += seg->value[i][grpA[j]];
    }
    x[i-a]/=dl1;
  }
  for(i=a; i<=b; i++) {
    for(j=0;j<noB;j++) {
      lb[0][((i-a)*noB)+j]=seg->value[i][grpB[j]];
      lb[1][((i-a)*noB)+j]=seg->pos[i];
      mean2+=seg->value[i][grpB[j]];
      y[i-a] += seg->value[i][grpB[j]];
    }
    y[i-a]/=dl2;
  }
  int u = mannwhitney( la[0], l*noA , lb[0], l*noB);
  p = mannwhitneyPvalue(u, l*noA, l*noB, nfo->MWU, MAXM, MAXN);
 // fprintf(stdout,"pVALUE %f\n",p);
  

  mean1/=dl1*dl;
  mean2/=dl2*dl;
  meandiff = mean1-mean2;
  seg->methA=-1.0;
  seg->methB=-1.0;
  
  kstest2d(la[0], la[1], l*noA, lb[0],lb[1], l*noB, faskstest);

  ks[0]=faskstest[0];
  ks[1]=meandiff;
  ks[2]=p;

  FREEMEMORY(NULL, la[0]);
  FREEMEMORY(NULL, la[1]);
  FREEMEMORY(NULL, lb[0]);
  FREEMEMORY(NULL, lb[1]);
  FREEMEMORY(NULL, la);
  FREEMEMORY(NULL, lb);
  FREEMEMORY(NULL, x);
  FREEMEMORY(NULL, y);

}

/*---------------------------------- kstest ----------------------------------
 *    
 * @brief calculated the means
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void concatFloatsToString(char **x, float *y, size_t size, char sep) {
    size_t total_length = 0;
    for (size_t i = 0; i < size; ++i) {
        total_length += snprintf(NULL, 0, "%.3f", y[i]) + 1;
    }
    total_length -= 1;

    *x = (char *)malloc(total_length + 1);

    char *ptr = *x;
    for (size_t i = 0; i < size; ++i) {
        int written = snprintf(ptr, total_length + 1, "%.3f", y[i]);
        ptr += written;
        total_length -= written;
        if (total_length > 0) {
          *ptr = sep;
          ++ptr;
          --total_length;
        }
    }
}

void outputImputedValues(char **x, double *y, size_t size, char sep) {
    int n_useless = 2;
    size_t total_length = 0;
    for (size_t i = 0; i < (size-n_useless); ++i) {
        total_length += snprintf(NULL, 0, "%.3f", y[i+n_useless]) + 1;
    }
    total_length -= 1;

    *x = (char *)malloc(total_length + 1);

    char *ptr = *x;
    for (size_t i = 0; i < (size-n_useless); ++i) {
        int written = snprintf(ptr, total_length + 1, "%.3f", y[i+n_useless]);
        ptr += written;
        total_length -= written;
        if (total_length > 0) {
          *ptr = sep;
          ++ptr;
          --total_length;
        }
    }
}

void means(segment_t *seg , int a, int b, int ***groupID, int **groupSize, int groupNumber, int **subgroupID, int *subgroupSize, int subgroupNumber, char **meansA, char **meansB){

  int i,j;
  float meanvalues[subgroupNumber];
  float submeanvalues[subgroupNumber];

  for (int sgn = 0; sgn < subgroupNumber; sgn++)
  {
    int *grpA = subgroupID[sgn];
    int noA = subgroupSize[sgn];
    float dl = 0;
    float mean=0;
    for(i=a; i<=b ; i++) {
      for(j=0; j < noA; j++) {
        mean+=seg->value[i][grpA[j]];
        dl+=1;
      }
    }
  
    mean/=dl;
    submeanvalues[sgn]=mean;
  }

  // int ci=0; // combination (more that one group) index
  // for (int gn = 0; gn < subgroupNumber; gn++)
  // {
  //   int *grpA = groupID[1][gn];
  //   int noA = groupSize[1][gn];
  //   float dlA = 0;
  //   float meanA=0;

  //   for(i=a; i<=b ; i++) {
  //     for(j=0; j < noA; j++) {
  //       meanA+=seg->value[i][grpA[j]];
  //       dlA+=1;
  //     }
  //   }
  
  //   meanA/=dlA;
  //   meanvalues[ci]=meanA;
  //   ci++;
  // }

  char *subtmp =NULL;
  concatFloatsToString(&subtmp, submeanvalues, subgroupNumber, '|');
  *meansA = subtmp;

  // char *tmp =NULL;
  // concatFloatsToString(&tmp, meanvalues, subgroupNumber, '|');
  // *meansB = tmp;
  *meansB = "TBC";

}


/*---------------------------- calcSingleDiffSum -----------------------------
 *    
 * @brief calculate single difference sum
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

double ***
calcSingleDiffSum(segment_t *s, int ***groupID, int **groupSize, int groupNumber, double mindiff, double mindiff2){

  int i, j;
  double ***S, s1, s2;
  
  S = ALLOCMEMORY(NULL, NULL, double**, groupNumber);
  // S[groupCombination][pos][0]: absolute
  // S[groupCombination][pos][1]: original
  // S[groupCombination][pos][2]: sign
  // S[groupCombination][pos][3]: accumulated # of significantly different CpGs
  for (int groupCombination = 0; groupCombination < groupNumber; groupCombination++)
  {
    int *grpA = groupID[0][groupCombination];
    int noA = groupSize[0][groupCombination];
    int *grpB = groupID[1][groupCombination];
    int noB = groupSize[1][groupCombination];

    S[groupCombination] = ALLOCMEMORY(NULL, NULL, double*, s->n);

    for(i=0; i < s->n; i++) {
      S[groupCombination][i] = ALLOCMEMORY(NULL, NULL, double, 5);
      memset(S[groupCombination][i], 0, sizeof(double)*5);
    }

    for(i=0; i < s->n; i++) {
      s1 = 0;
      for(int j=0;j<noA;j++)
        s1+=s->value[i][grpA[j]];
      s1/=noA;

      s2 = 0;
      for(j=0; j < noB; j++)
        s2+=s->value[i][grpB[j]];
      s2/=noB;

      if(i==0) {

        S[groupCombination][i][1]=s1-s2;
        S[groupCombination][i][0]=fabs(S[groupCombination][i][1]);
        if(s1-s2==0) {
          S[groupCombination][i][2] = 0;

        }
        else{
          if(s1-s2>0)
            S[groupCombination][i][2] = 1;
          else
            S[groupCombination][i][2] = -1;
        }
        if (fabs(s1-s2)>=mindiff)
        {
          S[groupCombination][i][3] = 1;
        } else {
          S[groupCombination][i][3] = 0;
        }

        if (fabs(s1-s2)>=mindiff2)
        {
          S[groupCombination][i][4] = 1;
        } else {
          S[groupCombination][i][4] = 0;
        }
        

      } else {

        S[groupCombination][i][1]=S[groupCombination][i-1][1]+s1-s2;
        S[groupCombination][i][0]=fabs(S[groupCombination][i][1]);
        if(s1-s2==0)
          S[groupCombination][i][2] = 0;
        else{
          if(s1-s2>0)
            S[groupCombination][i][2] = S[groupCombination][i-1][2]+1;
          else
            S[groupCombination][i][2] = S[groupCombination][i-1][2]-1;
        }   

        if (fabs(s1-s2)>=mindiff)
        {
          S[groupCombination][i][3] = S[groupCombination][i-1][3]+1;
        } else {
          S[groupCombination][i][3] = S[groupCombination][i-1][3];
        } 

        if (fabs(s1-s2)>=mindiff2)
        {
          S[groupCombination][i][4] = S[groupCombination][i-1][4]+1;
        } else {
          S[groupCombination][i][4] = S[groupCombination][i-1][4];
        } 

      }
    }
  }
  

  return S;
}

/*----------------------------- calcSigCpGs ------------------------------
 *    
 * @brief calculate # of significantly different CpGs
 * @author zzhu
 *   
 */

double 
calcSigCpGs(double **S,int s, int t) {
  double sigCpGs; 
  if(s==0)
    sigCpGs = S[t][3];
  else
    sigCpGs = S[t][3]-S[s-1][3];
  return sigCpGs;
}

double 
calcSigCpGs2(double **S,int s, int t) {
  double sigCpGs; 
  if(s==0)
    sigCpGs = S[t][4]/(t-s+1);
  else
    sigCpGs = (S[t][4]-S[s-1][4])/(t-s+1);
  return sigCpGs;
}

/*----------------------------- calcSingleTrend ------------------------------
 *    
 * @brief calculate single segment trend
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

double 
calcSingleTrendAbs(double **S,int s, int t) {
  double ds = s;
  double dt=t;
  double trend; 
  if(s==0)
    trend = fabs(S[t][2])/dt; // should divided by (dt+1)???
  else
    trend = fabs(S[t][2]-S[s-1][2])/(dt-ds+1);

  return trend;
}

double 
calcSingleTrendAbs2(double **S,int s, int t) {
  double ds = s;
  double dt=t;
  double trend; 
  if(s==0)
    trend = fabs(S[t][1]);
  else
    trend = fabs(S[t][1]-S[s-1][1]);

  return trend;
}

/*--------------------------------- noValley ---------------------------------
 *    
 * @brief check for local valley
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

char 
noValley(double **S, int s, int t, metseg_t *nfo){

  if(t-s+1< nfo->mincpgs) return 1;

  int i;
  double ds = s;
  double dt=t;
  double Sst=S[t][1];
  double Si, mean, imean;
  int minlength=(nfo->mincpgs>10)?nfo->mincpgs:10;
  
  if(s>0) {
    Sst-=S[s-1][1];
  }

  mean=fabs(Sst/(dt-ds+1));
//check if mean of windowsize=mincpgs is < valleyfactor*mean
  for(i=s; i+minlength-1 <= t; i++) {
    Si=S[i+minlength-1][1];
    if(i>0) Si-=S[i-1][1];
    imean=fabs(Si/(minlength));

//    if(fabs(imean)<mean*0.7) {
    if(fabs(imean)<mean* nfo->valley) {
      return 0;
    }    
  }
  return 1;
}



/*--------------------------------- findMaxZ ---------------------------------
 *    
 * @brief find maximum Z score
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

double 
findMaxZ(double **Z, int s, int t, int *ab, metseg_t *nfo) {
  int a, b;
  double max=0;

  //fprintf(stderr, "FindMaxN\t%d\t%d\n",s,t); 
  for(a=s; a<=t; a++)
    for(b=a+nfo->mincpgs-1; b<=t; b++) {
      if(((a!=s)||(b!=t) ) && Z[a][b]>max && (IsFiniteNumber(Z[a][b]))){
        max=Z[a][b];
        ab[0]=a;
        ab[1]=b;
      }  
    }
  //fprintf(stderr, "FindMaxNnew\t%d\t%d\n",ab[0],ab[1]); 
  return max;
}

double 
findRandomZ(int s, int t, int *ab, metseg_t *nfo) {
  int a, b;
  double max=0;

  for(a=s; a<=t; a++)
    for(b=a+nfo->mincpgs-1; b<=t; b++) {
        ab[0]=a;
        ab[1]=b;
        return max;
    }
  return max;
  //fprintf(stderr, "FindMaxNnew\t%d\t%d\n",ab[0],ab[1]); 
}


/*---------------------------- calcSingleDiffZabs ----------------------------
 *    
 * @brief calculate single Z score
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void 
calcSingleDiffZabs(double **S, int Ssize, double **Z, int s, int t, 
    char zero, metseg_t *nfo) {

  int a, b;
  double ds = s;
  double dt=t;
  double da, db, Sst, Sab, u;

  for(a=s; a<=t; a++){ 
    for(b=a+nfo->mincpgs-1; b<=t; b++) {
      da =a;
      db=b;
      Sst=S[t][1];
      if(s>0) { 
        Sst-=S[s-1][1];
      }
      Sab=S[b][1];
      if(a>0) { 
        Sab-=S[a-1][1];
      }
      Sab=fabs(Sab);
      Sst=fabs(Sst);

      if(db-da-nfo->mincpgs == 0 || (a==s && b==t) || a==b) { 
        Z[a-s][b-s]=0;
      } else {
        if(zero) {
          u = Sab;
        } else {
          u = Sab-(((db-da+1)*(Sst))/(dt-ds+1));
        }
        u*=u;
        u/=(db-da+1)*((1-((db-da+1)/(dt-ds+1))));
        Z[a-s][b-s]=u;
      }
    }
  }

  return;
}    


/*------------------------------ pushSegment_p -------------------------------
 *    
 * @brief a helper function to push segments to a stack
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
 
void
pushSegment_p (Segmentstack *stack, int a, int b, int *ab, int child,
    double *ks1, double *ks2, double *ks3, double *KS, segment_t *max)
{

    segment_p_t segment;

    segment.a = a;
    segment.b = b;
    if(child) { 
    segment.ab1 = ab[0];
    segment.ab2 = ab[1];
    segment.child = child;

    segment.ks11 =ks1[0];
    segment.ks12 =ks1[1];
    segment.ks13 =ks1[2];
    segment.ks14 =ks1[3];
 
    segment.ks21 =ks2[0];
    segment.ks22 =ks2[1];
    segment.ks23 =ks2[2];
    segment.ks24 =ks2[3];

    segment.ks31 =ks3[0];
    segment.ks32 =ks3[1];
    segment.ks33 =ks3[2];
    segment.ks34 =ks3[3];

    segment.KS1 =KS[0];
    segment.KS2 =KS[1];
    segment.KS3 =KS[2];
    segment.KS4 =KS[3]; 
    }
    if(max) { 
      segment_t *copy = ALLOCMEMORY(NULL, NULL, segment_t, 1);
      memmove(copy, max, sizeof(segment_t));
      segment.max = copy;
    }
    bl_segmentstackPush(stack, &segment);

	return ;
}

/*------------------------------- popSegment_p -------------------------------
 *    
 * @brief a helper function to pop segments from a stack
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
 
void
popSegment_p (Segmentstack *stack, int *a, int *b, int *ab, int *child, 
    double *ks1, double *ks2, double *ks3, double *KS, segment_t **max)
{
  segment_p_t *segment;

  segment = bl_segmentstackPop(stack);
     
  *a = segment->a;
  *b = segment->b;

  if(child) { 
    *child = segment->child;

    ab[0] = segment->ab1;
    ab[1] = segment->ab2;

    ks1[0] = segment->ks11;
    ks1[1] = segment->ks12;
    ks1[2] = segment->ks13;
    ks1[3] = segment->ks14;

    ks2[0] = segment->ks21;
    ks2[1] = segment->ks22;
    ks2[2] = segment->ks23;
    ks2[3] = segment->ks24;

    ks3[0] = segment->ks31;
    ks3[1] = segment->ks32;
    ks3[2] = segment->ks33;
    ks3[3] = segment->ks34;

    KS[0] = segment->KS1;
    KS[1] = segment->KS2;
    KS[2] = segment->KS3;
    KS[3] = segment->KS4;
  }

  if(max) { 
    *max = segment->max;
  }

  
  return;
}


/*------------------------------ stackSegment_p ------------------------------
 *    
 * @brief a helper function to initalize a segment stack
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
 
  Segmentstack*
stackSegment_p(Segmentstack *stack)
{
    bl_segmentstackInit(stack, 300);

	return stack;
}

/*---------------------------- concatStrings -----------------------------
 *    
 * @brief concat two strings to one
 * @author zzhu
 *   
 */

void concatStrings(char **x, char *y, char sep) {
    size_t x_length = (*x != NULL) ? strlen(*x) : 0;
    size_t y_length = strlen(y);
    size_t new_length = x_length + 1 + y_length + 1;

    // Allocate memory for the concatenated string
    char *concatenated = (char *)malloc(new_length);
    if (*x != NULL) {
        strcpy(concatenated, *x);
        strcat(concatenated, &sep);
    }
    strcat(concatenated, y);
    if (*x != NULL) {
        free(*x);
    }
    *x = concatenated;
}

void concatIntsToString(char **x, int *y, size_t size, char sep) {
    size_t total_length = 0;
    for (size_t i = 0; i < size; ++i) {
        total_length += snprintf(NULL, 0, "%d", y[i]) + 1;
    }
    total_length -= 1;

    *x = (char *)malloc(total_length + 1);

    char *ptr = *x;
    for (size_t i = 0; i < size; ++i) {
        int written = snprintf(ptr, total_length + 1, "%d", y[i]);
        ptr += written;
        total_length -= written;
        if (total_length > 0) {
          *ptr = sep;
          ++ptr;
          --total_length;
        }
    }
}


/*------------------------------ clustering ------------------------------
 *    
 * @brief given a 1 vs 1 comparison id, return 2 clusters for comparison
 * @author zzhu
 *   
 */

void
clustering(int ***clusters, int *nclusters, int numberSubGroup, int **subgroupID, int *subgroupSize, 
  segment_t *seg, metseg_t *nfo, double ***S, double *ks, int s, int t)
{
  if (ks[3]!=-1)
  {

    // Increase the number of rows in X by 1
    (*nclusters)++;
    // Reallocate memory for X
    *clusters = realloc(*clusters, (*nclusters) * sizeof(int *));

    // Allocate memory for the new row and append [1, 2, 3]
    (*clusters)[(*nclusters) - 1] = malloc(numberSubGroup * sizeof(int));

    int gn = ks[3];
    int A = 0;
    int B = 0;
    for (int i = 0; i < numberSubGroup; i++)
    {
      if (gn < ((2*numberSubGroup-i-1)*i/2)) {continue;}
      A = i;
      B = gn - (((2*numberSubGroup-i-1)*i/2) -i -1);
      
    }
    // assert(0);
    for (int g = 0; g < numberSubGroup; g++)
    {
      if (g==A) 
      {
        // fprintf(stderr, "A.\n");
        (*clusters)[(*nclusters) - 1][g] = 0;
      } else if (g==B)
      {
        // fprintf(stderr, "B.\n");
        (*clusters)[(*nclusters) - 1][g] = 4;
      } else
      {
        (*clusters)[(*nclusters) - 1][g] = 2;
      }
    }


    for (int g = 0; g < numberSubGroup; g++)
    {
      int i, j;
      if (g==A) {continue;}  
      if (g==B) {continue;}    
      if (g>A)
      {
        i = A; j = g;
      } else 
      {
        i = g; j = A;
      }
      int newgn = ((2*numberSubGroup-i-1)*i/2) + (j-i-1);

      if (nfo->minDMR2 != 1)
      {
        if (calcSigCpGs2(S[newgn], s, t) < nfo->minDMR2)
        {
          // fprintf(stderr, "A+.\n");
          (*clusters)[(*nclusters) - 1][g] -= 1;
        }
      } else {
        if (calcSigCpGs(S[newgn], s, t) < nfo->minDMR)
        {
          // fprintf(stderr, "A+.\n");
          (*clusters)[(*nclusters) - 1][g] -= 1;
        }
      }
      
 
        
    }

    for (int g = 0; g < numberSubGroup; g++)
    {
      int i, j;
      if (g==A) {continue;}  
      if (g==B) {continue;}    
      if (g>B)
      {
        i = B; j = g;
      } else 
      {
        i = g; j = B;
      }
      int newgn = ((2*numberSubGroup-i-1)*i/2) + (j-i-1);

      if (nfo->minDMR2 != 1)
      {
        if (calcSigCpGs2(S[newgn], s, t) < nfo->minDMR2)
        {
          // fprintf(stderr, "B+.\n");
          (*clusters)[(*nclusters) - 1][g] += 1;
        }
      } else {
        if (calcSigCpGs(S[newgn], s, t) < nfo->minDMR)
        {
          // fprintf(stderr, "B+.\n");
          (*clusters)[(*nclusters) - 1][g] += 1;
        }
      }

 
    }

    int **cgroupID;
    int *cgroupSize;

    cgroupID = ALLOCMEMORY(NULL, NULL, int*, 2);
    cgroupSize = ALLOCMEMORY(NULL, NULL, int, 2);
    int grpindex[2];
    grpindex[0] = 0;
    grpindex[1] = 0;
    for (int i = 0; i < 2; i++)
    {
      cgroupID[i] = NULL;
      // cgroupID[i] = ALLOCMEMORY(NULL, cgroupID[i], int*, numberSubGroup);
      cgroupSize[i] = 0;
    }
    for (int j = 0; j < numberSubGroup; j++)
    {
      int i;
      if ((*clusters)[(*nclusters) - 1][j]<2)
      {
        i = 0;
      }else if ((*clusters)[(*nclusters) - 1][j]>2)
      {
        i = 1;
      } else { continue; }
      
      cgroupID[i] = ALLOCMEMORY(NULL, cgroupID[i], int, cgroupSize[i]+subgroupSize[j]);
      for (int k = 0; k < subgroupSize[j]; k++)
      {
        cgroupID[i][grpindex[i]] = subgroupID[j][k];
        grpindex[i]++;
      }
      cgroupSize[i] += subgroupSize[j];
    }
    if (cgroupSize[0]>0 && cgroupSize[1]>0)
    {kstest(seg, s, t, 0, 1, 1, ks, cgroupID[0], cgroupSize[0], cgroupID[1], cgroupSize[1], nfo);}
    
    ks[3] = (*nclusters) - 1;
    // fprintf(stderr,"ks test done%f\n",ks[3]);
    for (int i = 0; i < 2; i++)
    {
      FREEMEMORY(NULL, cgroupID[i]);
    }
    FREEMEMORY(NULL, cgroupID);
    FREEMEMORY(NULL, cgroupSize);
  }
  
  
}

/*-------------------------------- segment_p ---------------------------------
 *    
 * @brief calculate segment probability
 * @author Frank Juehling and Steve Hoffmann  
 *   
 */


segment_t*
segment_pSTKopt(segment_t *seg, segment_t *breaks, int *nbreaks, double ***XS, 
    int a, int b, double *KS, char first, 
    int ***groupID, int **groupSize, int groupNumber, int **subgroupID, int *subgroupSize, int ***clusters, int *nclusters, metseg_t *nfo) {
// zzhu$ a,b: the start pos and end pos of the seg. ab: {start of subseg, end of subseg}. KS: the test result for the seg.
  int i, n, m;
  int dimZ, child = 0; // zzhu$ child: the pointer to a child, start from the left child
  int ab[3] = {-1,0,-1}; // {a,b,group with min P}
  int ab_tmp[2] = {-1,0};
  double ks1[] = {2,0,2,-1}; // zzhu$ {ks p, mean diff, ranksums p,group with min P } for all combinations
  double ks2[] = {2,0,2,-1};
  double ks3[] = {2,0,2,-1};
  double ks1_tmp[] = {2,0,2}; // zzhu$ {ks p, mean diff, ranksums p} for one combination
  double ks2_tmp[] = {2,0,2};
  double ks3_tmp[] = {2,0,2};
  // double ks1[groupNumber][3];
  // for (int i = 0; i < groupNumber; i++)
  // {
  //   ks1[i][0] = 2;
  //   ks1[i][1] = 0;
  //   ks1[i][2] = 2;
  // }
  // double ks2[groupNumber][3];
  // for (int i = 0; i < groupNumber; i++)
  // {
  //   ks2[i][0] = 2;
  //   ks2[i][1] = 0;
  //   ks2[i][2] = 2;
  // }
  // double ks3[groupNumber][3];
  // for (int i = 0; i < groupNumber; i++)
  // {
  //   ks3[i][0] = 2;
  //   ks3[i][1] = 0;
  //   ks3[i][2] = 2;
  // }
  double init_kstmp[] = {2,0,2};
  double init_ks[] = {2,0,2,-1};

  double **XZ;
  double newp; // prob, meandiff;

  Segmentstack stack;
  stackSegment_p(&stack);

  while(!bl_segmentstackIsEmpty(&stack) || a != -1) { 

    if(a != -1 && child <= 2) { 
      if(ab[0] == -1) {  // zzhu$ if ab is not calculated yet
        double Z_max = -1;
        int ab_updated = 0; // if no KS is performed, ab will be updated with max Z
        int ab_Zmax[3] = {-1,0,-1}; // if no KS is performed, ab will be updated with max Z

        memmove(ks1, init_ks, sizeof(double)*4);
        memmove(ks2, init_ks, sizeof(double)*4);
        memmove(ks3, init_ks, sizeof(double)*4);

        int existSigGn = 0;
        for(int gn=0;gn<groupNumber;gn++){

          ab_tmp[0] = 0;
          ab_tmp[1] = 0;

          memmove(ks1_tmp, init_kstmp, sizeof(double)*3);
          memmove(ks2_tmp, init_kstmp, sizeof(double)*3);
          memmove(ks3_tmp, init_kstmp, sizeof(double)*3);


          // double ks1_tmp[] = {2,0,2};
          // double ks2_tmp[] = {2,0,2};
          // double ks3_tmp[] = {2,0,2};

          // fprintf(stderr,"init ks.%f,%f,%f\n", ks1_tmp[0],ks2_tmp[0],ks3_tmp[0]);


          // add a filter step here
          if (calcSigCpGs(XS[gn], a, b) < nfo->minDMR)
          {
            if ((gn==(groupNumber-1))&&(existSigGn==0))
            {
                // fprintf(stderr,"no,%d,%d\n",a,b);
                dimZ = b - a + 1;
                // XZ = ALLOCMEMORY(NULL, NULL, double*, dimZ);

                // for(i=0; i < dimZ; i++) {
                //   XZ[i] = ALLOCMEMORY(NULL, NULL, double, dimZ);
                //   memset(XZ[i], 0, sizeof(double)*dimZ);
                // }

                // calcSingleDiffZabs(XS[gn], seg->n, XZ , a, b, 0, nfo);
                double Z_max_tmp = findRandomZ(0, dimZ-1, ab_tmp, nfo);
                ab_tmp[0]+=a; 
                ab_tmp[1]+=a;
                ab_Zmax[0] = ab_tmp[0];
                ab_Zmax[1] = ab_tmp[1];
                ab_Zmax[2] = gn;
                // for(i=0; i < dimZ; i++) {
                // FREEMEMORY(NULL, XZ[i]);
                // }
                // // fprintf(stderr,"no,%d,%d\n",ab_tmp[0],ab_tmp[1]);
                // FREEMEMORY(NULL, XZ);
            }
            continue;
          }
          // fprintf(stderr,"yes\n");
          // end a filter step here

          dimZ = b - a + 1;
          XZ = ALLOCMEMORY(NULL, NULL, double*, dimZ);

          for(i=0; i < dimZ; i++) {
            XZ[i] = ALLOCMEMORY(NULL, NULL, double, dimZ);
            memset(XZ[i], 0, sizeof(double)*dimZ);
          }

          //calculated the Z scores and find the maximum interval
          calcSingleDiffZabs(XS[gn], seg->n, XZ , a, b, 0, nfo);
          double Z_max_tmp = findMaxZ(XZ, 0, dimZ-1, ab_tmp, nfo);
          for(i=0; i < dimZ; i++) {
            FREEMEMORY(NULL, XZ[i]);
          }
          FREEMEMORY(NULL, XZ);
          ab_tmp[0]+=a; 
          ab_tmp[1]+=a;

          if (Z_max_tmp>Z_max)
          {
            Z_max=Z_max_tmp;
            ab_Zmax[0] = ab_tmp[0];
            ab_Zmax[1] = ab_tmp[1];
            ab_Zmax[2] = gn;
          }

          if (nfo->clustering == 0)
          {  
            //check the left side of the maximum interval with ks
            n=a; m=ab_tmp[0]-1;
            if(ab_tmp[0] > 0 && m-n+1 >= nfo->mincpgs 
                && calcSingleTrendAbs(XS[gn],a,ab_tmp[0]-1) > nfo->trend 
                && noValley(XS[gn], a, ab_tmp[0]-1, nfo)) {

              kstest(seg, a, ab_tmp[0]-1, 0, 1, 1, ks1_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
            }

            //check the maximum interval interval with ks
            n=ab_tmp[0];m=ab_tmp[1];
            if(m-n+1 >= nfo->mincpgs 
                && calcSingleTrendAbs(XS[gn],ab_tmp[0],ab_tmp[1]) > nfo->trend 
                && noValley(XS[gn], ab_tmp[0], ab_tmp[1], nfo) ) {

              kstest(seg, ab_tmp[0], ab_tmp[1], 0, 1, 1, ks2_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
            }

            //check the right side of the maximum interval with ks
            n=ab_tmp[1]+1;m=b;
            if(m-n+1 >= nfo->mincpgs 
                && calcSingleTrendAbs(XS[gn],ab_tmp[1]+1,b)> nfo->trend 
                && noValley(XS[gn], ab_tmp[1]+1, b, nfo)) {

              kstest(seg, ab_tmp[1]+1, b, 0, 1, 1, ks3_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
            }

            if (MIN(ks1_tmp[0],MIN(ks2_tmp[0],ks3_tmp[0]))<MIN(ks1[0],MIN(ks2[0],ks3[0])))
            {
              // fprintf(stderr,"yes.");
              for ( i = 0; i < 3; i++)
              {
                ks1[i] = ks1_tmp[i];
                ks2[i] = ks2_tmp[i];
                ks3[i] = ks3_tmp[i];
              }
              ks1[3] = gn;
              ks2[3] = gn;
              ks3[3] = gn;
              ab[0] = ab_tmp[0];
              ab[1] = ab_tmp[1];
              ab[2] = gn;
              ab_updated = 1;
              // fprintf(stderr,"ab updated. a:%d,b:%d,ks1:%f,ks2:%f,ks3:%f\n", a,b, ks1[0],ks2[0],ks3[0]);
            }
          } // else
          // {  
          //   //check the left side of the maximum interval with ks
          //   // n=a; m=ab_tmp[0]-1;
          //   // if(ab_tmp[0] > 0 && m-n+1 >= nfo->mincpgs 
          //   //     && calcSingleTrendAbs(XS[gn],a,ab_tmp[0]-1) > nfo->trend 
          //   //     && noValley(XS[gn], a, ab_tmp[0]-1, nfo)) {

          //   //   kstest_fake(seg, a, ab_tmp[0]-1, 0, 1, 1, ks1_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
          //   // }

          //   // //check the maximum interval interval with ks
          //   // n=ab_tmp[0];m=ab_tmp[1];
          //   // if(m-n+1 >= nfo->mincpgs 
          //   //     && calcSingleTrendAbs(XS[gn],ab_tmp[0],ab_tmp[1]) > nfo->trend 
          //   //     && noValley(XS[gn], ab_tmp[0], ab_tmp[1], nfo) ) {

          //   //   kstest_fake(seg, ab_tmp[0], ab_tmp[1], 0, 1, 1, ks2_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
          //   // }

          //   // //check the right side of the maximum interval with ks
          //   // n=ab_tmp[1]+1;m=b;
          //   // if(m-n+1 >= nfo->mincpgs 
          //   //     && calcSingleTrendAbs(XS[gn],ab_tmp[1]+1,b)> nfo->trend 
          //   //     && noValley(XS[gn], ab_tmp[1]+1, b, nfo)) {

          //   //   kstest_fake(seg, ab_tmp[1]+1, b, 0, 1, 1, ks3_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
          //   // }

          //   if (Z_max_tmp>Z_max)
          //   {
          //     ks1[3] = gn;
          //     ks2[3] = gn;
          //     ks3[3] = gn;
          //   }
          // } 
          existSigGn++;
        }

        if(existSigGn==0){
          ab[0] = ab_Zmax[0];
          ab[1] = ab_Zmax[1];
          ab[2] = ab_Zmax[2];
          ks1[3] = ab_Zmax[2];
          ks2[3] = ab_Zmax[2];
          ks3[3] = ab_Zmax[2];
          // fprintf(stderr,"ab Z-updated. a:%d,b:%d,ks1:%f,ks2:%f,ks3:%f\n", ab[0],ab[1], ks1[0],ks2[0],ks3[0]);
        } else {
          if(ab_updated==0){
            ab[0] = ab_Zmax[0];
            ab[1] = ab_Zmax[1];
            ab[2] = ab_Zmax[2];
            ks1[3] = ab_Zmax[2];
            ks2[3] = ab_Zmax[2];
            ks3[3] = ab_Zmax[2];
            // fprintf(stderr,"ab Z-updated. a:%d,b:%d,ks1:%f,ks2:%f,ks3:%f\n", a,b, ks1[0],ks2[0],ks3[0]);
          }
        }
        
        // fprintf(stderr,"ab:%d,%d",ab[0] ,ab[1] );
                // re-calculate the KS for all subsegments
        if(existSigGn>0){
          if (nfo->clustering == 0)
          {
            for(int gn=0;gn<groupNumber;gn++){
              memmove(ks1_tmp, init_kstmp, sizeof(double)*3);
              memmove(ks2_tmp, init_kstmp, sizeof(double)*3);
              memmove(ks3_tmp, init_kstmp, sizeof(double)*3);
              
              // add a filter step here
              if (calcSigCpGs(XS[gn], a, b) < nfo->minDMR)
              {
                continue;
              }
              // end a filter step here
              
              //check the left side of the maximum interval with ks
              n=a; m=ab[0]-1;
              if(ab[0] > 0 && m-n+1 >= nfo->mincpgs 
                  && calcSingleTrendAbs(XS[gn],a,ab[0]-1) > nfo->trend 
                  && noValley(XS[gn], a, ab[0]-1, nfo)) {

                kstest(seg, a, ab[0]-1, 0, 1, 1, ks1_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
              }

              //check the maximum interval interval with ks
              n=ab[0];m=ab[1];
              if(m-n+1 >= nfo->mincpgs 
                  && calcSingleTrendAbs(XS[gn],ab[0],ab[1]) > nfo->trend 
                  && noValley(XS[gn], ab[0], ab[1], nfo) ) {

                kstest(seg, ab[0], ab[1], 0, 1, 1, ks2_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
              }

              //check the right side of the maximum interval with ks
              n=ab[1]+1;m=b;
              if(m-n+1 >= nfo->mincpgs 
                  && calcSingleTrendAbs(XS[gn],ab[1]+1,b)> nfo->trend 
                  && noValley(XS[gn], ab[1]+1, b, nfo)) {
                kstest(seg, ab[1]+1, b, 0, 1, 1, ks3_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
              }

              if (ks1_tmp[0]<ks1[0])
              {
                for ( i = 0; i < 3; i++)
                {
                  ks1[i] = ks1_tmp[i];
                }
                ks1[3] = gn;
              }

              if (ks2_tmp[0]<ks2[0])
              {
                for ( i = 0; i < 3; i++)
                {
                  ks2[i] = ks2_tmp[i];
                }
                ks2[3] = gn;
              }

              if (ks3_tmp[0]<ks3[0])
              {
                for ( i = 0; i < 3; i++)
                {
                  ks3[i] = ks3_tmp[i];
                }
                ks3[3] = gn;
              }
            }
          }
            
          if (nfo->clustering == 1)
          {
            // to be replaced by:
            // void clustering(int ***clusters, int *nclusters, int numberSubGroup, int **subgroupID, int *subgroupSize, 
            // segment_t *seg, metseg_t *nfo, double ***S, double *ks, int s, int t)

            n=a; m=ab[0]-1;
            if(ab[0] > 0 && m-n+1 >= nfo->mincpgs) {
              clustering(clusters, nclusters, nfo->groups, subgroupID, subgroupSize, seg, nfo, XS, ks1, a, ab[0]-1);
            }
            n=ab[0];m=ab[1];
            if(m-n+1 >= nfo->mincpgs) {
              clustering(clusters, nclusters, nfo->groups, subgroupID, subgroupSize, seg, nfo, XS, ks2, ab[0], ab[1]);
            }
            n=ab[1]+1;m=b;
            if(m-n+1 >= nfo->mincpgs) {
              clustering(clusters, nclusters, nfo->groups, subgroupID, subgroupSize, seg, nfo, XS, ks3, ab[1]+1, b);
            }

          }
        }
        // fprintf(stderr,"Final: a:%d,b:%d,ks1:%f,ks2:%f,ks3:%f\n", a,b, ks1[0],ks2[0],ks3[0]);
      }

      //if one of the children has a good ks we check all
      newp = MIN(ks1[0],MIN(ks2[0],ks3[0]));
      // assert(((newp<1)&&(ab[1]-ab[0]>=nfo->mincpgs)));
      //if the p-value of one of the three intervals is better than the original
      //recurse down to find the best subinterval
      if((newp<KS[0] || (newp>1 && KS[0]>1)) && (b-a >= nfo->mincpgs)) { // zzhu$ termination criteria (1)#CpG (2)p_child>p_parent. no need for newp>1???
        // assert((b-ab[1]+1>=10)||(ks3[0]>1));
        // fprintf(stderr,"here\n");
        pushSegment_p (&stack, a, b, ab, child+1, ks1, ks2, ks3, KS, NULL); // zzhu$ push the next child into the stack.

        //left interval child
        if(child == 0) {
          n = a;
          m = ab[0]-1;
          a = -1; // zzhu$ end this child if meet termination criteria
          if(ab[0] > 0 && n <= m) {
            if(m-n >= nfo->mincpgs) { 
              a = n;
              b = m;
              memmove(KS, ks1, sizeof(double)*4);
              child = 0;
              ab[0] = -1; // zzhu$ ab[0]==-1 means ab to be calculated.
            } else { 
              // for(int gn=0;gn<groupNumber;gn++){
              //     double ks4tmp[] = {2,0,2};
              //     kstest(seg, n, m, 0, 1, 1, ks4tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
              //     fprintf(stderr,"left: a%d,b%d,ks%f\tgroup%d,ks%f\n", a, b, ks1[0], gn, ks4tmp[0]);
              // }
              breaks = ALLOCMEMORY(NULL, breaks, segment_t, (*nbreaks)+1);
              setSegment(&breaks[(*nbreaks)], seg->chr, n, m, ks1[0], ks1[1], ks1[2], ks1[3]);
              (*nbreaks)+=1;
            }
          }
        }

        //middle interval child
        if(child == 1) {
          n = ab[0];
          m = ab[1];
          a = -1;
          if(n <= m) {
            if(m-n>= nfo->mincpgs) { 
              a = n;
              b = m;
              memmove(KS, ks2, sizeof(double)*4);
              child = 0;
              ab[0] =-1;
            } else {
              // for(int gn=0;gn<groupNumber;gn++){
              //     double ks4tmp[] = {2,0,2};
              //     kstest(seg, n, m, 0, 1, 1, ks4tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
              //     fprintf(stderr,"middle: a%d,b%d,ks%f\tgroup%d,ks%f\n", a, b, ks2[0], gn, ks4tmp[0]);
              // }
              breaks = ALLOCMEMORY(NULL, breaks, segment_t, (*nbreaks)+1);
              setSegment(&breaks[(*nbreaks)], seg->chr, n, m, ks2[0], ks2[1], ks2[2], ks2[3]);
              (*nbreaks)+=1;
            }  
          }
        }

        //right interval child
        if(child == 2) {
          n = ab[1]+1;
          m = b;
          a = -1;
          if(n<=m) {
            if (m-n>= nfo->mincpgs) { 
              a = n;
              b = m;
              memmove(KS, ks3, sizeof(double)*4);
              child = 0;
              ab[0] = -1;
            } else {
              // for(int gn=0;gn<groupNumber;gn++){
              //     double ks4tmp[] = {2,0,2};
              //     kstest(seg, n, m, 0, 1, 1, ks4tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
              //     fprintf(stderr,"right: a%d,b%d,ks%f\tgroup%d,ks%f\n", a, b, ks3[0], gn, ks4tmp[0]);
              // }
              breaks = ALLOCMEMORY(NULL, breaks, segment_t, (*nbreaks)+1);
              setSegment(&breaks[(*nbreaks)], seg->chr, n, m, ks3[0], ks3[1], ks3[2], ks3[3]);
              (*nbreaks)+=1;
            } 
          }
        }
      } else {
        if(child ==0) { 
          // for(int gn=0;gn<groupNumber;gn++){
          //   double ks4tmp[] = {2,0,2};
          //   kstest(seg, a, b, 0, 1, 1, ks4tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
          //   // calcSingleTrendAbs(XS[gn],a,ab[0]-1) > nfo->trend 
          //   //   && noValley(XS[gn], a, ab[0]-1, nfo)
          //   fprintf(stderr,"KS. a%d,b%d,ks%f\tgroup%d,ks%f,trend%f,novally%d\n", a, b, KS[0], gn, ks4tmp[0],calcSingleTrendAbs(XS[gn],a,b),noValley(XS[gn], a, b, nfo));
          // }
          // assert((a!=15)||(b!=26));
          
          breaks = ALLOCMEMORY(NULL, breaks, segment_t, (*nbreaks)+1);
          setSegment(&breaks[(*nbreaks)], seg->chr, a, b, KS[0], KS[1], KS[2], KS[3]);
          (*nbreaks)+=1;
        }

        a = -1;
        ab[0] = -1;
      }
      
    }  else {

      popSegment_p (&stack, &a, &b, ab, &child, ks1, ks2, ks3, KS, NULL);
      // fprintf(stderr,"a%d,b%d,ab%d,ab%d,ks1%f,ks2%f,ks3%f\n", a,b, ab[0], ab[1], ks1[0], ks2[0], ks3[0]);

      if(child == 3) a = -1;

    }
  }

  bl_segmentstackDestruct(&stack);

  return breaks;
}
/*-------------------------------- segmenter ---------------------------------
 *    
 * @brief first segmenter function
 * @author Frank Juehling and Steve Hoffmann 
 *    
 */
segment_t*
segmenterSTK(segment_t *seg, segment_t *globalbreaks, int *nglobal, double ***XS, 
    int a, int b, double *KS, int ***groupID, int **groupSize, int groupNumber, int **subgroupID, int *subgroupSize, int ***clusters, int *nclusters,
    metseg_t *nfo) {
  
  // zzhu$ a: the start pos of the region, b is the end pos.

  int nbreaks=0, i, n, m; //, *s, *t;
  double trend[groupNumber];
  segment_t *breaks=NULL, *max; // zzhu$ breaks: pre-segments. max: the break with most significant p value.

  Segmentstack stack; 
  stackSegment_p(&stack);
  //inorder traversal of tree of intervals
  //while(!bl_vstackIsEmpty(stack) || a != -1) {
  while(!bl_segmentstackIsEmpty(&stack) || a != -1) {

    if(a != -1) { // zzhu$ initial a will be 0, from the function segmentation 
      //push current interval node
      nbreaks = 0;
      breaks = NULL;
      breaks = segment_pSTKopt(seg, breaks, &nbreaks, XS, a, b, KS, 0, 
          groupID, groupSize, groupNumber, subgroupID, subgroupSize, clusters, nclusters, nfo);

      max = &breaks[0]; // zzhu$ max: the break with most significant p value.
      // zzhu$ the most significant break as DMR.
      for(i=0; i < nbreaks; i++) { 
        if(max->prob > breaks[i].prob) { 
          max = &breaks[i];
        }
      }

      pushSegment_p (&stack, a, b, NULL, 0, NULL, NULL, NULL, NULL, max); 

      //set next
      n = a;
      m = max->start-1;

      if(max->start > 0 && n <= m) {
        double ks[] = {2,0,2,-1}; // ks p, meandiff, ranksums p, group with min p
        double ks_tmp[] = {2,0,2};
        for (int gn = 0; gn < groupNumber; gn++)
        {
          // add a filter step here
          if (calcSigCpGs(XS[gn], n, m) < nfo->minDMR)
          {
            continue;
          }
          // end a filter step here
          trend[gn] = calcSingleTrendAbs(XS[gn], n, m);
          if(m-n+1 >= nfo->mincpgs && trend[gn] > nfo->trend && noValley(XS[gn], n, m, nfo)) { 
            kstest(seg, n, m, 0, 1, 1, ks_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
            if (ks_tmp[0]<ks[0])
            {
              ks[0]=ks_tmp[0];
              ks[1]=ks_tmp[1];
              ks[2]=ks_tmp[2];
              ks[3]=gn;
            }
          }
        }
        if(m-n+1 >= nfo->mincpgs){clustering(clusters, nclusters, nfo->groups, subgroupID, subgroupSize, seg, nfo, XS, ks, n, m);}
        a = n;
        b = m;
        KS = ks;
      } else {
        a = -1; // zzhu$ a==-1 means there is no left region of max, continue to right segment.
      }

      FREEMEMORY(NULL, breaks);

    } else {

       
      popSegment_p (&stack, &a, &b, NULL, 0, NULL, NULL, NULL, NULL, &max);

      //here the segment is registered!
      globalbreaks = ALLOCMEMORY(NULL, globalbreaks, segment_t, (*nglobal)+1);
      memmove(&globalbreaks[(*nglobal)], max, sizeof(segment_t)); // zzhu$ globalbreaks: the DMR set
      (*nglobal)+=1;

      n = max->stop+1;
      m = b;
      if(n<=m) {
        // trend = calcSingleTrendAbs(XS,n,m);
        // double ks[] = {2,0,2};
        // if(m-n+1 >= nfo->mincpgs && trend > nfo->trend && noValley(XS, n, m, nfo)) { 
        //   kstest(seg, n, m, 0, 1, 1, ks, groupID, groupSize, groupNumber, nfo);
        // }
        double ks[] = {2,0,2,-1}; // ks p, meandiff, ranksums p, group with min p
        double ks_tmp[] = {2,0,2};
        for (int gn = 0; gn < groupNumber; gn++)
        {
          // add a filter step here
          if (calcSigCpGs(XS[gn], n, m) < nfo->minDMR)
          {
            continue;
          }
          // end a filter step here

          trend[gn] = calcSingleTrendAbs(XS[gn], n, m);
          if(m-n+1 >= nfo->mincpgs && trend[gn] > nfo->trend && noValley(XS[gn], n, m, nfo)) { 
            kstest(seg, n, m, 0, 1, 1, ks_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
            if (ks_tmp[0]<ks[0])
            {
              ks[0]=ks_tmp[0];
              ks[1]=ks_tmp[1];
              ks[2]=ks_tmp[2];
              ks[3]=gn;
            }
          }
        }
        if(m-n+1 >= nfo->mincpgs){clustering(clusters, nclusters, nfo->groups, subgroupID, subgroupSize, seg, nfo, XS, ks, n, m);}
        a = n;
        b = m;
        KS = ks;
      } else {
        a = -1; 
      }

      FREEMEMORY(NULL, max);
    }
  }

  bl_segmentstackDestruct(&stack); 
  return globalbreaks;
}

/*---------------------------------- convert_sigcp2string ----------------------------------
 *    
 * @brief convert segment-specific sigcp to string
 * @author zzhu
 *   
 */
void convert_sigcp2string(int nclusters, int sigcp, int **clusters, int subgroupNumber, char **meansB){
  char *subtmp =NULL;
  if (sigcp>=0 && sigcp<nclusters)
  {
      concatIntsToString(&subtmp, clusters[sigcp], subgroupNumber, '|');
      *meansB = subtmp;
  }
}

/*---------------------------------- output ----------------------------------
 *    
 * @brief output routine
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void
output(segment_t *seg, segment_t *breaks, int nglobal, double ***XS, 
    int ***groupID, int **groupSize, int groupNumber, 
    int **subgroupID, int *subgroupSize, int subgroupNumber, int ***clusters, int *nclusters,
    metseg_t *nfo) {
  
  // fprintf(stderr,"start out2*****:\t%d\n",*nclusters);
  if(nfo->outputList->i >= nfo->outputList->n) {
        nfo->outputList->n+=1000000;
        nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
  }

  // for (int cl = 0; cl < nclusters; cl++)
  // {
  //         for (int i = 0; i < subgroupNumber; i++)
  //   {
  //     fprintf(stderr, "output%d,\t",(*clusters)[cl][i]);
  //   }
  // }
        
  int i;
  double trend;
  segment_t *b, *tmp=NULL;
  nfo->outputList->numberTests+=nglobal;
  for(i=0; i<nglobal; i++) {

    b = &breaks[i];

    if(b->prob > 1) {
      if(tmp == NULL) {
        // fprintf(stderr,"FirstNUL\n");
        tmp = ALLOCMEMORY(NULL, NULL, segment_t, 1);
        tmp->chr = b->chr;
        tmp->start=b->start;
        tmp->stop=b->stop;
        // tmp->prob=b->prob;
        // tmp->meandiff=b->meandiff;
        // tmp->test=b->test;
        tmp->sigcp=b->sigcp;
        
        char *me[] = {"-2","-2"};
        // fprintf(stderr, "Output 0: \t%s\n", me[0]);
        // fprintf(stderr, "Output 02: \t%s\n", me[1]);
        means(seg, tmp->start,tmp->stop, groupID, groupSize, groupNumber, subgroupID, subgroupSize, subgroupNumber, &me[0], &me[1]);
        // fprintf(stderr, "Output 1: \t%s\n", me[0]);
        // fprintf(stderr, "Output 11: \t%f,%d\n", b->sigcp, *nclusters);
      //         for (int i = 0; i < subgroupNumber; i++)
      // {
      //   fprintf(stderr, "sadfa%d,\t",(*clusters)[(int)tmp->sigcp][i]);
      // }

        convert_sigcp2string(*nclusters, (int)tmp->sigcp, *clusters, subgroupNumber, &me[1]);
        // fprintf(stderr, "Output 1: \t%s\n", me[0]);
        // fprintf(stderr, "Output 12: \t%s\n", me[1]);
        
        // tmp->methA=me[0];
        // tmp->methB=me[1]; //important!!!
        

      } else { // continue to find next break with p>1
        tmp->stop=b->stop;
      }

    } else { // if b.pval<1:

        
  // ks[0]=faskstest[0];
 // ks[1]=meandiff;
 // ks[2]=p;
     
        
        
        
      if(tmp != NULL) { // combine all segments until a segment with p<1;

        
        // double ks[] = {2,0,2};
        double ks[] = {2,0,2,-1}; // ks:best ks pval. {ks p, meandiff, ranksums p, group with min p}
        double ks_tmp[] = {2,0,2};
        int existSigGn = 0;
        for(int gn=0; gn<groupNumber; gn++){
          // add a filter step here
          if (calcSigCpGs(XS[gn], tmp->start,tmp->stop) < nfo->minDMR)
          {
            continue;
          }
          // end a filter step here

          trend = calcSingleTrendAbs(XS[gn],tmp->start,tmp->stop);
          if(tmp->stop-tmp->start + 1 >= nfo->mincpgs && trend>nfo->trend 
              && noValley(XS[gn], tmp->start, tmp->stop, nfo)) {

            kstest(seg, tmp->start,tmp->stop, 0, 1, 1, ks_tmp, 
                groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
            if (ks_tmp[0]<ks[0])
            {
              ks[0] = ks_tmp[0];
              ks[1] = ks_tmp[1];
              ks[2] = ks_tmp[2];
              ks[3] = gn;
            }
          }
          existSigGn++;
        }
        if (existSigGn>0 && tmp->stop-tmp->start + 1 >= nfo->mincpgs){
          clustering(clusters, nclusters, nfo->groups, subgroupID, subgroupSize, seg, nfo, XS, ks, tmp->start,tmp->stop);
        }

        if(ks[0]<2) {
            nfo->outputList->segment_out[nfo->outputList->i].chr = ALLOCMEMORY(NULL, NULL, char, strlen(seg->chr)+1);
            nfo->outputList->segment_out[nfo->outputList->i].chr = strcpy(nfo->outputList->segment_out[nfo->outputList->i].chr,seg->chr);
            nfo->outputList->segment_out[nfo->outputList->i].start = seg->pos[tmp->start]-1;
            // if ((tmp->stop-tmp->start+1+1)<10)
            // {
              // fprintf(stderr,"start1:%f,%d\n",(ks[0]),(tmp->stop-tmp->start+1));
            //   assert(0);
            // }
            
            // assert((ks[0]<1)&&((tmp->stop-tmp->start+1)<nfo->mincpgs));
            // if(nfo->outputList->i>=2){fprintf(stderr,"start1:%d,%d",seg->pos[tmp->start]-1,nfo->outputList->segment_out[2].start);}
            nfo->outputList->segment_out[nfo->outputList->i].stop = seg->pos[tmp->stop];
            nfo->outputList->segment_out[nfo->outputList->i].p = ks[0];
            nfo->outputList->segment_out[nfo->outputList->i].q = -1;
            nfo->outputList->segment_out[nfo->outputList->i].meandiff = ks[1];
            nfo->outputList->segment_out[nfo->outputList->i].mwu = ks[2];
            nfo->outputList->segment_out[nfo->outputList->i].length = (tmp->stop-tmp->start+1);
            nfo->outputList->segment_out[nfo->outputList->i].sigcp = ks[3];
            
            char *me[] = {"-2","-2"};
            means(seg,tmp->start,tmp->stop, groupID, groupSize, groupNumber, subgroupID, subgroupSize, subgroupNumber, &me[0], &me[1]);
            convert_sigcp2string(*nclusters, (int)tmp->sigcp, *clusters, subgroupNumber, &me[1]);
            nfo->outputList->segment_out[nfo->outputList->i].methA = ALLOCMEMORY(NULL, NULL, char, strlen(me[0])+1);
            nfo->outputList->segment_out[nfo->outputList->i].methB = ALLOCMEMORY(NULL, NULL, char, strlen(me[1])+1);
            nfo->outputList->segment_out[nfo->outputList->i].methA = strcpy(nfo->outputList->segment_out[nfo->outputList->i].methA,me[0]);
            nfo->outputList->segment_out[nfo->outputList->i].methB = strcpy(nfo->outputList->segment_out[nfo->outputList->i].methB,me[1]);
            
            nfo->outputList->i+=1;
            if(nfo->outputList->i >= nfo->outputList->n) {
                nfo->outputList->n+=1000000;
                nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
            }

        } 

        FREEMEMORY(NULL, tmp);
        tmp=NULL;
      }

      // if (( seg->pos[b->stop]-seg->pos[b->start]+1+1)<10)
      // {
        // fprintf(stderr,"start12:%f,%d\n",(b->prob),( seg->pos[b->stop]-seg->pos[b->start]+1+1));
      //   assert(0);
      // }
      
      nfo->outputList->segment_out[nfo->outputList->i].chr = ALLOCMEMORY(NULL, NULL, char, strlen(seg->chr)+1);
      nfo->outputList->segment_out[nfo->outputList->i].chr = strcpy(nfo->outputList->segment_out[nfo->outputList->i].chr,seg->chr);
      nfo->outputList->segment_out[nfo->outputList->i].start = seg->pos[b->start]-1;
      nfo->outputList->segment_out[nfo->outputList->i].stop = seg->pos[b->stop];
      nfo->outputList->segment_out[nfo->outputList->i].p = b->prob;
      nfo->outputList->segment_out[nfo->outputList->i].q = -1;
      nfo->outputList->segment_out[nfo->outputList->i].meandiff = b->meandiff;
      nfo->outputList->segment_out[nfo->outputList->i].mwu = b->test;
      nfo->outputList->segment_out[nfo->outputList->i].length = (b->stop-b->start+1);
      nfo->outputList->segment_out[nfo->outputList->i].sigcp = b->sigcp;
      // fprintf(stderr,"start12:%f,%d\n",(b->prob),( seg->pos[b->stop]-seg->pos[b->start]+1+1));
      char *me[] = {"-2","-2"};
      means(seg, b->start,b->stop, groupID, groupSize, groupNumber,subgroupID, subgroupSize, subgroupNumber, &me[0], &me[1]);
// fprintf(stderr,"start12:%f,%d\n",(b->prob),( seg->pos[b->stop]-seg->pos[b->start]+1+1));
      // for (int i = 0; i < subgroupNumber; i++)
      // {
      //   fprintf(stderr, "sadfa%d,\t",(*clusters)[(int)b->sigcp][i]);
      // }

      convert_sigcp2string(*nclusters, (int)b->sigcp, *clusters, subgroupNumber, &me[1]);
      // fprintf(stderr,"start12:%f,%d\n",(b->prob),( seg->pos[b->stop]-seg->pos[b->start]+1+1));
      nfo->outputList->segment_out[nfo->outputList->i].methA = ALLOCMEMORY(NULL, NULL, char, strlen(me[0])+1);
      nfo->outputList->segment_out[nfo->outputList->i].methB = ALLOCMEMORY(NULL, NULL, char, strlen(me[1])+1);
      nfo->outputList->segment_out[nfo->outputList->i].methA = strcpy(nfo->outputList->segment_out[nfo->outputList->i].methA,me[0]);
      nfo->outputList->segment_out[nfo->outputList->i].methB = strcpy(nfo->outputList->segment_out[nfo->outputList->i].methB,me[1]);
          
      // if(nfo->outputList->i>=2){fprintf(stderr,"start2:%d,%d\n",nfo->outputList->segment_out[2].start,nfo->outputList->segment_out[nfo->outputList->i].stop);}
      
      nfo->outputList->i+=1;
      if(nfo->outputList->i >= nfo->outputList->n) {
          nfo->outputList->n+=1000000;
          nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
      }
      // fprintf(stderr,"start13:%f,%d\n",(b->prob),( seg->pos[b->stop]-seg->pos[b->start]+1+1));
    }
  }


  if(tmp != NULL) {
    // fprintf(stderr,"NULL\n");
    // trend = calcSingleTrendAbs(XS,tmp->start,tmp->stop);
    // double ks[] = {2,0,2};

    // if(tmp->stop-tmp->start + 1 >= nfo->mincpgs && trend > nfo->trend 
    //     && noValley(XS, tmp->start, tmp->stop, nfo)) {
    //   kstest(seg,tmp->start,tmp->stop,0, 1, 1, ks, groupID, groupSize, groupNumber, nfo);
    // }

    double ks[] = {2,0,2,-1}; // ks:best ks pval. {ks p, meandiff, ranksums p, group with min p}
    double ks_tmp[] = {2,0,2};
    int existSigGn = 0;
    for(int gn=0; gn<groupNumber; gn++){
      // add a filter step here
      if (calcSigCpGs(XS[gn], tmp->start,tmp->stop) < nfo->minDMR)
      {
        continue;
      }
      // end a filter step here
      
      trend = calcSingleTrendAbs(XS[gn],tmp->start,tmp->stop);
      if(tmp->stop-tmp->start + 1 >= nfo->mincpgs && trend>nfo->trend 
          && noValley(XS[gn], tmp->start, tmp->stop, nfo)) {

        kstest(seg, tmp->start,tmp->stop, 0, 1, 1, ks_tmp, 
            groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
        if (ks_tmp[0]<ks[0])
        {
          ks[0] = ks_tmp[0];
          ks[1] = ks_tmp[1];
          ks[2] = ks_tmp[2];
          ks[3] = gn;
        }
      }
      existSigGn++;
    }
    if (existSigGn>0 && tmp->stop-tmp->start + 1 >= nfo->mincpgs){
      clustering(clusters, nclusters, nfo->groups, subgroupID, subgroupSize, seg, nfo, XS, ks, tmp->start,tmp->stop);
    }

    if(ks[0]<2) {
        // fprintf(stderr,"NULL:\n");
        // if ((seg->pos[tmp->stop]-seg->pos[tmp->start]+1+1)<10)
        // {
        //   // fprintf(stderr,"start1:%f,%d\n",(ks[0]),( seg->pos[tmp->stop]-seg->pos[tmp->start]+1+1));
        //   assert(0);
        // }
        nfo->outputList->segment_out[nfo->outputList->i].chr = ALLOCMEMORY(NULL, NULL, char, strlen(seg->chr)+1);
        nfo->outputList->segment_out[nfo->outputList->i].chr = strcpy(nfo->outputList->segment_out[nfo->outputList->i].chr,seg->chr);
        nfo->outputList->segment_out[nfo->outputList->i].start = seg->pos[tmp->start]-1;
        // fprintf(stderr,"start3:");
        // if(nfo->outputList->i>=2){fprintf(stderr,"start3:%d,%d",seg->pos[tmp->start]-1,nfo->outputList->segment_out[2].start);}
        nfo->outputList->segment_out[nfo->outputList->i].stop = seg->pos[tmp->stop];
        nfo->outputList->segment_out[nfo->outputList->i].p = ks[0];
        nfo->outputList->segment_out[nfo->outputList->i].q = -1;
        nfo->outputList->segment_out[nfo->outputList->i].meandiff = ks[1];
        nfo->outputList->segment_out[nfo->outputList->i].mwu = ks[2];
        nfo->outputList->segment_out[nfo->outputList->i].length = (tmp->stop-tmp->start+1);
        nfo->outputList->segment_out[nfo->outputList->i].sigcp = ks[3];
        
        char *me[] = {"-2","-2"};
        means(seg, tmp->start,tmp->stop, groupID, groupSize, groupNumber, subgroupID, subgroupSize, subgroupNumber, &me[0], &me[1]);
        convert_sigcp2string(*nclusters, (int)tmp->sigcp, *clusters, subgroupNumber, &me[1]);
        nfo->outputList->segment_out[nfo->outputList->i].methA = ALLOCMEMORY(NULL, NULL, char, strlen(me[0])+1);
        nfo->outputList->segment_out[nfo->outputList->i].methB = ALLOCMEMORY(NULL, NULL, char, strlen(me[1])+1);
        nfo->outputList->segment_out[nfo->outputList->i].methA = strcpy(nfo->outputList->segment_out[nfo->outputList->i].methA,me[0]);
        nfo->outputList->segment_out[nfo->outputList->i].methB = strcpy(nfo->outputList->segment_out[nfo->outputList->i].methB,me[1]);
        
        nfo->outputList->i+=1;
        if(nfo->outputList->i >= nfo->outputList->n) {
            nfo->outputList->n+=1000000;
            nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
            // if(nfo->outputList->i>=2){fprintf(stderr,"start4:%d",nfo->outputList->segment_out[2].start);}
        }
    }
    FREEMEMORY(NULL, tmp);
    tmp=NULL;
  }

  return;
}

//these variables need to be volatile to stop the compiler from optimizing
//as they are manipulated by the threads!
static volatile unsigned int idle;
volatile char *schedule = NULL;
static pthread_mutex_t out;
static pthread_mutex_t cnt;

/*------------------------------- segmentation -------------------------------
 *    
 * @brief do the segmentation
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */


int 
segmentation(char **chr, int *pos, double **value, int n, 
    int ***groupID, int **groupSize, int groupNumber, 
    int **subgroupID, int *subgroupSize,
    metseg_t *nfo) {

  double ***S;
  double ks[] = {2,0,2,-1}; // ks:best ks pval. {ks p, meandiff, ranksums p, group with min p}
  double ks_tmp[] = {2,0,2};
  // double ks[groupNumber][3];
  // for (int i = 0; i < groupNumber; i++)
  // {
  //   ks[i][0] = 2;
  //   ks[i][1] = 0;
  //   ks[i][2] = 2;
  // }
  double trend[groupNumber];
  char novalley[groupNumber];
  segment_t *seg, *global=NULL; // zzhu$ seg: containing all met. values in the two groups. global: the DMRs from segmenterSTK.
  int i, nglobal = 0;

  int **clusters = NULL;
  int nclusters = 0;

  seg = ALLOCMEMORY(NULL, NULL, segment_t, 1);

  initSegment(seg);
  seg->n = n; // zzhu$ number of CpGs in the region
  seg->chr = chr[0];
  seg->pos = pos;
  seg->value = value;

  S = calcSingleDiffSum(seg, groupID, groupSize, groupNumber, nfo->mindiff, nfo->mindiff2);
  int existSigGn = 0;
  for (int gn = 0; gn < groupNumber; gn++)
  {
    // add a filter step here
    if (calcSigCpGs(S[gn], 0, n-1) < nfo->minDMR)
    {
      // fprintf(stderr,"no sigcpgs gn:%d\n",gn);
      continue;
    }
    // end a filter step here
    trend[gn] = calcSingleTrendAbs(S[gn], 0, n-1);
    novalley[gn] = noValley(S[gn], 0, n-1, nfo);
    if(seg->n-1 >= nfo->mincpgs && trend[gn] > nfo->trend && novalley[gn]) { // should be seg->n >= mincpg???
      kstest(seg , 0, n-1, 0, 1, 1, ks_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
      if (ks_tmp[0]<ks[0])
      {
        ks[0]=ks_tmp[0];
        ks[1]=ks_tmp[1];
        ks[2]=ks_tmp[2];
        ks[3]=gn;
      }
    }
    existSigGn++;
  }
  if (existSigGn>0 && seg->n-1 >= nfo->mincpgs){
    // fprintf(stderr,"*****:\t%d\n",existSigGn);
    clustering(&clusters, &nclusters, nfo->groups, subgroupID, subgroupSize, seg, nfo, S, ks, 0, n-1);
    global = segmenterSTK(seg, global, &nglobal, S, 0, n-1, ks, 
      groupID, groupSize, groupNumber, subgroupID, subgroupSize, &clusters, &nclusters, nfo);
  }
  

  //set lock if necessary 
  if(nfo->threads >1) { 
    pthread_mutex_lock(&out);
  }
  //output here   
  // for (int cl = 0; cl<nclusters; cl++) {
  //   for (int i = 0; i < nfo->groups; i++)
  //   {
  //     fprintf(stderr, "%d,\t", clusters[cl][i]);
  //   }
  //   fprintf(stderr, "\n");
  // }
  // fprintf(stderr,"nfo->mincpgs:%d", nfo->mincpgs);
  output(seg, global, nglobal, S, groupID, groupSize, groupNumber, subgroupID, subgroupSize, nfo->groups, &clusters, &nclusters, nfo);
  
  //unlock if necessary
  if(nfo->threads > 1) {
    pthread_mutex_unlock(&out);
  }

  for(int gn=0;gn<groupNumber;gn++){
    for(i=0; i < seg->n; i++) {
      FREEMEMORY(NULL, S[gn][i]);
    }
    FREEMEMORY(NULL, S[gn]);
  }

//     for (int i = 0; i < nclusters; i++) {
//         for (int j = 0; j < 3; j++) {
// //            fprintf(stderr,"%d ", clusters[i][j]);
//         }
//         printf("\n");
//     }
  for (int i = 0; i < nclusters; i++) {
      // fprintf(stderr,"rm%d\n",i);
    free(clusters[i]);
  }
  free(clusters);
  


  FREEMEMORY(NULL, global);
  FREEMEMORY(NULL, seg);
  FREEMEMORY(NULL, S);
  // if(nfo->outputList->i>=2){fprintf(stderr,"start6:%d\n",nfo->outputList->segment_out[2].start);}
  return 0;
}



/*-------------------------------- segworker ---------------------------------
 *    
 * @brief for threaded segmentation
 * @author Frank Juehling and Steve Hoffma2nn 
 *   
 */
 
void*
segworker (void *args)
{
  int i;
  metseg_t *t;
  t = (metseg_t*) args;
   
  segmentation(t->chr, t->pos, t->value, t->n, t->groupID, t->groupSize, t->groupNumber, t->subgroupID, t->subgroupSize, t);
  
  //cleanup own data
  for(i=0; i < t->n; i++) {
    FREEMEMORY(NULL, t->chr[i]);
    FREEMEMORY(NULL, t->value[i]);
  }
  FREEMEMORY(NULL, t->chr);
  FREEMEMORY(NULL, t->pos);
  FREEMEMORY(NULL, t->value);


  pthread_mutex_lock(&cnt);
  schedule[t->threadno] = 0;
  //decrement of idle at creation time
  idle = idle+1;
  pthread_mutex_unlock(&cnt);

  pthread_exit(args);
}

/*------------------------------- regionTest -------------------------------
 *    
 * @brief test a region
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void 
regionTest(segment_t *seg, int ***groupID, int **groupSize, int groupNumber, 
int **subgroupID, int *subgroupSize, metseg_t *nfo) {
    double ks[] = {2,0,2,-1};
    double ks_tmp[] = {2,0,2};
    double ***S;
    int **clusters = NULL;
    int nclusters = 0;
    double maxZ = 0;
    if(seg->n>0) {
      if (nfo->clustering == 1)
      {
          S = calcSingleDiffSum(seg, groupID, groupSize, groupNumber, nfo->mindiff, nfo->mindiff2);
      }
      for (int gn = 0; gn < groupNumber; gn++)
      {
        if (nfo->clustering == 0)
        {
          kstest(seg , 0, seg->n-1, 0, 0, 1, ks_tmp, groupID[0][gn], groupSize[0][gn], groupID[1][gn], groupSize[1][gn], nfo);
          if (ks_tmp[0]<ks[0])
          {
            for (int i = 0; i < 3; i++)
            {
              ks[i] = ks_tmp[i];
            }
            ks[3] = gn;
          } 
        }
        if (nfo->clustering == 1)
        {
          double newZ = calcSingleTrendAbs2(S[gn], 0, seg->n-1);
          // fprintf(stderr,"Z%d,%f\n",gn,newZ);
          if (maxZ<newZ)
          {
            maxZ = newZ;
            ks[3] = gn;
          } 
        }
      }
      if (nfo->clustering == 1)
      {
        // fprintf(stderr,"ks test done%f,%f\n",ks[3],seg->sigcp);
        clustering(&clusters, &nclusters, nfo->groups, subgroupID, subgroupSize, seg, nfo, S, ks, 0, seg->n-1);
        // fprintf(stderr,"ks test done%f,%f\n",ks[3],seg->sigcp);
        for (int gn = 0; gn < groupNumber; gn++)
        {
          for (int i = 0; i < seg->n; i++)
          {
            FREEMEMORY(NULL, S[gn][i]);
          }
          FREEMEMORY(NULL, S[gn]);
        }
      }
    }
    char *me[] = {"-2","-2"};
    means(seg, 0, seg->n-1,groupID, groupSize, groupNumber, subgroupID, subgroupSize, nfo->groups, &me[0],&me[1]);
    if (nfo->clustering == 1)
    {
      convert_sigcp2string(nclusters, ks[3], clusters, nfo->groups, &me[1]);
    }
//    void kstest(segment_t *seg , int a, int b, char mindiff, char mincpgs, char test, 
//  (segment_t *seg , int a, int b, char mindiff, char mincpgs, char test, 
//    double *ks, int *grpA, int noA, int *grpB, int noB, metseg_t* nfo){
    
    
    if(nfo->threads >1) { 
      pthread_mutex_lock(&out);
    }
 
/*    
  //output here   
    fprintf(stdout, "#%s\t%d\t%d\t%.2g\t%f\t%d\t%.2g\n", 
	    seg->chr, seg->start-1, 
	    seg->stop, ks[0], ks[1],
	    seg->n, ks[2]);
    fflush(stdout); 
  */  
    
    if(nfo->outputList->i >= nfo->outputList->n-1) {
      nfo->outputList->n+=1000000;
      nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
    }
    
    nfo->outputList->segment_out[nfo->outputList->i].chr = ALLOCMEMORY(NULL, NULL, char, strlen(seg->chr)+1);
    nfo->outputList->segment_out[nfo->outputList->i].chr = strcpy(nfo->outputList->segment_out[nfo->outputList->i].chr,seg->chr);
    nfo->outputList->segment_out[nfo->outputList->i].start = seg->start-1;
    nfo->outputList->segment_out[nfo->outputList->i].stop = seg->stop;
    nfo->outputList->segment_out[nfo->outputList->i].p = ks[0];
    nfo->outputList->segment_out[nfo->outputList->i].meandiff = ks[1];
    nfo->outputList->segment_out[nfo->outputList->i].mwu = ks[2];
    nfo->outputList->segment_out[nfo->outputList->i].length = seg->n;
    nfo->outputList->segment_out[nfo->outputList->i].methA = me[0];
    nfo->outputList->segment_out[nfo->outputList->i].methB = me[1];
    nfo->outputList->segment_out[nfo->outputList->i].sigcp = ks[3];
    
    nfo->outputList->i+=1;
    nfo->outputList->numberTests+=1;
    
    
    
    
    
    //unlock if necessary
    if(nfo->threads > 1) {
    pthread_mutex_unlock(&out);
  }
    
    
    
    

    destructSegment(seg);
    return;
}
void 
cpgTest(char *chr, int start, int stop, double ratio, double p, metseg_t *nfo, double methA, double methB) {
//cpg->chr, cpg->start,cpg->stop,ratio,p);    
    
    if(nfo->threads >1) { 
    pthread_mutex_lock(&out);
  }
    if(nfo->outputList->i >= nfo->outputList->n-1) {
      nfo->outputList->n+=1000000;
      nfo->outputList->segment_out = ALLOCMEMORY(NULL, nfo->outputList->segment_out, segment_out, nfo->outputList->n);
    }
    nfo->outputList->segment_out[nfo->outputList->i].chr = ALLOCMEMORY(NULL, NULL, char, strlen(chr)+1);
    nfo->outputList->segment_out[nfo->outputList->i].chr = strcpy(nfo->outputList->segment_out[nfo->outputList->i].chr,chr);
    nfo->outputList->segment_out[nfo->outputList->i].start = start-1;
    nfo->outputList->segment_out[nfo->outputList->i].stop = stop;
    //nfo->outputList->segment_out[nfo->outputList->i].p = p;
    nfo->outputList->segment_out[nfo->outputList->i].meandiff = ratio;
    nfo->outputList->segment_out[nfo->outputList->i].mwu = p;
    nfo->outputList->segment_out[nfo->outputList->i].length = 1;
    // nfo->outputList->segment_out[nfo->outputList->i].methA = methA;
    // nfo->outputList->segment_out[nfo->outputList->i].methB = methB; // important!!!
    
    nfo->outputList->i+=1;
    nfo->outputList->numberTests+=1;
   //unlock if necessary
    if(nfo->threads > 1) {
    pthread_mutex_unlock(&out);
    }
    return;
}
/*-------------------------------- segworker_CpG -----------------------------
 *    
 * @brief for threaded segmentation in single CpG mode
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
void*
segworker_CpG (void *args)
{
    
  metseg_t *t;
  t = (metseg_t*) args;
  cpg_t *cpg = t->cpg;
  
  
  int ua = mannwhitney (cpg->groupA, cpg->noA, cpg->groupB, cpg->noB);
  double p= mannwhitneyPvalue(ua, cpg->noA, cpg->noB, t->MWU, MAXM, MAXN);
  //double p = mannwhitney (cpg->groupA, cpg->noA, cpg->groupB, cpg->noB);
  double ratio = get_meandiff(cpg, cpg->groupA, cpg->noA, cpg->groupB, cpg->noB);
  
  
  
  
  /*
  //output here   
  fprintf(stdout, "#%s\t%d\t%d\t.\t1\t%f\n", 
              cpg->chr, cpg->start-1, 
              cpg->stop,p);
  */
  cpgTest(cpg->chr, cpg->start,cpg->stop,ratio,p,t, cpg->methA,cpg->methB);
  
  
  destructCpg(t->cpg);
  
  pthread_mutex_lock(&cnt);
  schedule[t->threadno] = 0;
  //decrement of idle at creation time
  idle = idle+1;
  pthread_mutex_unlock(&cnt);

  pthread_exit(args);
}
/*-------------------------------- segworker_region --------------------------
 *    
 * @brief for threaded segmentation in region mode
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
void*
segworker_region (void *args)
{
  metseg_t *t;
  t = (metseg_t*) args;
  regionTest(t->seg, t->groupID, t->groupSize, t->groupNumber, t->subgroupID, t->subgroupSize, t);
          
  pthread_mutex_lock(&cnt);
  schedule[t->threadno] = 0;
  //decrement of idle at creation time
  idle = idle+1;
  pthread_mutex_unlock(&cnt);

  pthread_exit(args);
}



/*-------------------------------- checkSetNAN ---------------------------------
 *    
 * @brief for checking NaNs in input data
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */
 
int
checkSetNAN(stringset_t **csv, double *values){
//check line for NANs and set . to NAN
    int nan=0;
    values[0] = -1;
    values[1] = -1;
    for(int k=2; k < csv[0]->noofstrings; k++) { 
        char *s = csv[0]->strings[k].str;
        values[k] = atof(s);
        if(strcmp(".", s) == 0 || strcmp("-", s) == 0 || strcmp("", s) == 0 || strcmp(" ", s) == 0){
            values[k] = NAN;
        }
        int in=0;
        while (values[k] == values[k] && s[in]) {
          if (isalpha(s[in])) {
              values[k] = NAN;
              break;
          }
//          else printf ("character %c is not alphabetic\n",s[in]);
          in++;
        }  
        if(values[k] != values[k] ) {
            nan+=1;
        }
    }
    return nan;
}
        

   
// /*-------------------------------- fillNAN ---------------------------------
//  *    
//  * @brief for replacing NaNs with betaDist
//  * @author Frank Juehling and Steve Hoffmann 
//  *   
//  */
 
int 
fillNAN(double *values, int **subgroupID, int *subgroupSize, int numberSubGroup, metseg_t *nfo) {
  //    fprintf(stderr,"#Ueberhaupt schaetzen\n");
  for (int gn = 0; gn < numberSubGroup; gn++)
  {  
    int na=0;
    int noA = subgroupSize[gn];
    double *groupA = ALLOCMEMORY(NULL, NULL, double, noA);
    double varA;
    double meanA = 0.0;
    int j;

//    fprintf(stdout,"\ngroupA\t");
    j=0;
    for(int i=0; i<noA; i++) 
        if(values[subgroupID[gn][i]+2] == values[subgroupID[gn][i]+2]) {
          groupA[j] = values[subgroupID[gn][i]+2];
 //         fprintf(stdout,"%f\t",groupA[j]);
          na+=1;
          meanA+=groupA[j];
          j++;
    }
 //   fprintf(stdout,"\ngroupB\t");
   
    if(na<1 || na<nfo->minNoA) {
        FREEMEMORY(NULL, groupA);
	      //  fprintf(stderr,"#REMOVING POSITION CUTOFF\n");
        return 1;
    }
    //    fprintf(stderr,"#NOT REMOVING POSITION CUTOFF\n");
    meanA/=(double)na;
    if(na == 1) {
        varA=0.000001;
    }
    else {
        varA = var(groupA,na);
    }
  //  fprintf(stderr,"Group %d, #meanA: %f\tvarA: %f\n",gn,meanA,varA);
    
    
    if(meanA < 0.000001)
        meanA = 0.000001;
   
  //  fprintf(stderr,"#new As:\n");
    for(int i=0; i<noA; i++) {
      if(isnan(values[subgroupID[gn][i]+2])) {
          values[subgroupID[gn][i]+2] = rbeta_mv(meanA, varA);
          while(isnan(values[subgroupID[gn][i]+2])) {
//            values[subgroupID[gn][i]+2] = meanA;
            values[subgroupID[gn][i]+2] = rbeta_mv(meanA, varA);
	    //            fprintf(stderr,"#betaAnot %d\n",i);
          }
          //else {
            //  fprintf(stdout,"#betaA %d\n",i);
       //   }
        }
      //  fprintf(stderr,"%f\t",values[subgroupID[gn][i]+2]);
    }
    // fprintf(stderr,"\n");


//   fprintf(stdout,"#A:\t%f / %f\t%f / %f\n",meanA,varA,meanB,varB);
    
    FREEMEMORY(NULL, groupA);
  }
  return 0;
    
 //   double
//var (double *x, int n)
    
    
}








/*---------------------------- initProgramParams -----------------------------
 *    
 * @brief initialize program parameters
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */

void
initProgramParams (metseg_t *nfo)
{
  nfo->maxdist=300;
  nfo->maxseg=-1;
  nfo->mincpgs=10;
  nfo->threads=1;
  nfo->mode=1;
  nfo->mtc=1;
  // nfo->nameA = "g1";
  // nfo->nameB = "g2";
  nfo->groups = 2; // newcodes
  nfo->minDMR = 1; // newcodes
  nfo->mindiff = 0; // newcodes
  nfo->minDMR2 = 1; // newcodes
  nfo->mindiff2 = 0; // newcodes
  nfo->clustering = 0; // newcodes
  nfo->outputImputed = 0; // newcodes
  nfo->trend = 0.6;
  nfo->minNoA = -1;
  nfo->minNoB = -1;
  nfo->minFactor = 0.8;
  nfo->valley = 0.7;
  nfo->minMethDist = 0.1;
  nfo->MWU = generateMannWhitneyCDFMatrix(MAXM, MAXN);
  //  testMannWhitneyApprox (MAXM, MAXN, nfo->MWU);
  nfo->randomseed=26061981;
  return ;
}

// /*---------------------------- initProgramParams2 -----------------------------
//  *    
//  * @brief initialize program parameters that depend on group sizes
//  * @author Frank Juehling and Steve Hoffmann 
//  *   
//  */

// void
// initProgramParams2 (metseg_t *nfo, int noA, int noB)
// {
//   if(nfo->minNoA<0) {
//         nfo->minNoA = ceil((double)noA * nfo->minFactor);
//   }
  
//   if(nfo->minNoB<0) {
//         nfo->minNoB = ceil((double)noB * nfo->minFactor);
//   }
  
//   return ;
// }





/*---------------------------- calGroupNumber -----------------------------
 *    
 * @brief find all possible combinations(groups)
 * @author zzhu
 *   
 */

int calGroupNumber(int n, int ***grpA_subgroups, int ***grpB_subgroups, int clustering){
  if (clustering==0)
  {
    int Nc = (pow(3,n)-pow(2,n+1)+1)/2; // number of possible combinations
    // fprintf(stderr, "# Combination %d\n",Nc);
    int **A_subgroups;
    int **B_subgroups;
    A_subgroups = ALLOCMEMORY(NULL, NULL, int*, Nc);
    B_subgroups = ALLOCMEMORY(NULL, NULL, int*, Nc);
    for (int i = 0; i < Nc; i++)
    {
      A_subgroups[i] = NULL;
      A_subgroups[i] = ALLOCMEMORY(NULL, A_subgroups[i], int, n);
      B_subgroups[i] = NULL;
      B_subgroups[i] = ALLOCMEMORY(NULL, B_subgroups[i], int, n);
      for (int j = 0; j < n; j++)
      {
        A_subgroups[i][j] = 0;
        B_subgroups[i][j] = 0;
      }
    }
    int ii = 0; // index for effective combinations
    for (int i = 0; i < pow(3,n); i++)
    {
      int i_copy = i;
      int sumA = 0;
      int sumB = 0;
      int validAB=0; // to exclude duplicates
      for (int j = 0; j < n; j++)
      {
        if (i_copy%3 == 1) {
          A_subgroups[ii][j] = 1;
          sumA++;
          validAB = 0;
        }
        if (i_copy%3 == 2) {
          B_subgroups[ii][j] = 1;
          sumB++;
          validAB = 1;
        }
        i_copy /= 3;
      }
      
      if ((sumA>0)&&(sumB>0)&&(validAB))
      {
        fprintf(stdout, "# Combination %d:\tGroup A subgroups: ",ii);
        for (int j = 0; j < n; j++)
        {
          if (A_subgroups[ii][j] == 1) {fprintf(stdout, "%d,",j);}
        }
        fprintf(stdout, "\tGroup B subgroups: ");
        for (int j = 0; j < n; j++)
        {
          if (B_subgroups[ii][j] == 1) {fprintf(stdout, "%d,",j);}
        }
        fprintf(stdout, "\n");
        ii++;
        if (ii==Nc)
        {
          break;
        }
      } else {
        for (int j = 0; j < n; j++)
        {
          A_subgroups[ii][j] = 0;
          B_subgroups[ii][j] = 0;
        }
      }
    }
    *grpA_subgroups = A_subgroups;
    *grpB_subgroups = B_subgroups;
    return Nc;
  }
  else {
    int Nc = (n*(n-1))/2; // number of possible combinations
    // fprintf(stderr, "# Combination %d\n",Nc);
    int **A_subgroups;
    int **B_subgroups;
    A_subgroups = ALLOCMEMORY(NULL, NULL, int*, Nc);
    B_subgroups = ALLOCMEMORY(NULL, NULL, int*, Nc);
    for (int i = 0; i < Nc; i++)
    {
      A_subgroups[i] = NULL;
      A_subgroups[i] = ALLOCMEMORY(NULL, A_subgroups[i], int, n);
      B_subgroups[i] = NULL;
      B_subgroups[i] = ALLOCMEMORY(NULL, B_subgroups[i], int, n);
      for (int j = 0; j < n; j++)
      {
        A_subgroups[i][j] = 0;
        B_subgroups[i][j] = 0;
      }
    }
    int ii = 0; // index for effective combinations
    for (int i = 0; i < n; i++)
    {
      int i_copy = i;
      for (int j = i+1; j < n; j++)
      {
        A_subgroups[ii][i] = 1;
        B_subgroups[ii][j] = 1;
        // fprintf(stdout, "# Combination %d:\tGroup A subgroups: %d,\tGroup B subgroups: %d,\n",ii,i,j);
        ii++;
      }
    }
    *grpA_subgroups = A_subgroups;
    *grpB_subgroups = B_subgroups;
    return Nc;
  }
  
}



/*---------------------------- selectGroups -----------------------------
 *    
 * @brief select meaningful combinations(groups)
 * @author zzhu
 *   
 */

// int selectGroups(int Nc, int ***grpA_subgroups, int ***grpB_subgroups){
//   /* TBC */
//   return 0;
// }

/*----------------------------------- main -----------------------------------
 *    
 * @brief the main routine
 * @author Frank Juehling and Steve Hoffmann 
 *   
 */


int main(int argc, char** argv) {

  manopt_optionset optset;
  manopt_arg *args;
  manopt_intconstraint modeconstraint;
  manopt_intconstraint mtcconstraint;
  manopt_intconstraint clusteringconstraint;
  metseg_t nfo; // zzhu$ nfo(metseg_t): parameters for the whole process and input data.
  metseg_t *th_nfo;
  stringset_t **csv, **bedcsv, **headercsv; // zzhu$ input table
  fileiterator_t *fi, *bedfi, *headerfi;
  unsigned int i, j, k, ln, bedln, headerln;
  pthread_t *threads = NULL;
  pthread_attr_t tattr;
  char *bedfile = NULL;
  char *headerfile = NULL;


  char **chr = NULL;
  int *grpA = NULL, noA=0;
  int *grpB = NULL, noB=0;
  double **val = NULL;
  int *pos = NULL;

  modeconstraint.min = 1;
  modeconstraint.max = 3;
  mtcconstraint.min = 1;
  mtcconstraint.max = 2;
  clusteringconstraint.min = 0;
  clusteringconstraint.max = 1;
  
  initProgramParams(&nfo);
  int verbose = 0;

  int numnonnan = 0;

   //we want to detach the threads to automatically have their resources freed
  pthread_attr_init(&tattr); 
  pthread_attr_setdetachstate(&tattr,PTHREAD_CREATE_DETACHED);


  pthread_mutex_init(&out, NULL);
  pthread_mutex_init(&cnt, NULL);

// Options  
  manopt_initoptionset(&optset, argv[0], NULL, 
      "metilene - a tool for fast and sensitive detection of differential DNA methylation\n\nDataInputFile\t\tneeds to be SORTED for chromosomes and genomic positions",
      "Implemented by Frank Juehling and Steve Hoffmann\n  2015-2016 Bioinformatik Leipzig\n",
      version,
      "Please report bugs to [frank,steve]@bioinf.uni-leipzig.de");

  manopt(&optset, REQUINTOPT, 0, 'M', "maxdist", 
      "maximum distance", "<n>", NULL, &nfo.maxdist);
  manopt(&optset, REQUINTOPT, 0, 'G', "maxseg", 
      "maximum segment length in case of memory issues", "<n>", NULL, &nfo.maxseg);
  manopt(&optset, REQUINTOPT, 0, 'm', "mincpgs", 
      "minimum cpgs", "<n>", NULL, &nfo.mincpgs);
  manopt(&optset, REQDBLOPT, 0, 'd', "minMethDiff", 
      "minimum mean methylation difference", "<n>", NULL, &nfo.minMethDist);
  manopt(&optset, REQUINTOPT, 0, 't', "threads", 
      "number of threads", "<n>", NULL, &nfo.threads);
  manopt(&optset, REQUINTOPT, 0, 'f', "mode", 
      "number of method: 1: de-novo, 2: pre-defined regions, 3: DMCs", "<n>", &modeconstraint, &nfo.mode);
  manopt(&optset, REQUINTOPT, 0, 'c', "mtc", 
      "method of multiple testing correction: 1: Bonferroni, 2: Benjamini-Hochberg (FDR)", "<n>", &mtcconstraint, &nfo.mtc);
  // manopt(&optset, REQSTRINGOPT, 0, 'a', "groupA", 
  //     "name of group A", "<string>", NULL, &nfo.nameA);
  // manopt(&optset, REQSTRINGOPT, 0, 'b', "groupB", 
  //     "name of group B", "<string>", NULL, &nfo.nameB);
  manopt(&optset, REQSTRINGOPT, 0, 'B', "bed", 
      "bed-file for mode 2 containing pre-defined regions; needs to be SORTED equally to the DataInputFile", "<string>", NULL, &bedfile);
  manopt(&optset, REQUINTOPT, 0, 'X', "minNoA", 
      "minimal number of values in group A", "<n>", NULL, &nfo.minNoA);
  manopt(&optset, REQUINTOPT, 0, 'Y', "minNoB", 
      "minimal number of values in group B", "<n>", NULL, &nfo.minNoB);
  manopt(&optset, REQUINTOPT, 0, 's', "seed",
      "set seed for random generator", "<n>", NULL, &nfo.randomseed);
  manopt(&optset, REQDBLOPT, 0, 'v', "valley", 
      "valley filter (0.0 - 1.0)", "<n>", NULL, &nfo.valley);

  manopt(&optset, REQUINTOPT, 0, 'n', "groups", 
      "number of groups", "<n>", NULL, &nfo.groups);
  manopt(&optset, REQUINTOPT, 0, 'r', "minDMR", 
      "minimal DMR", "<n>", NULL, &nfo.minDMR);
  manopt(&optset, REQDBLOPT, 0, 'w', "mindiff", 
      "minimal difference", "<n>", NULL, &nfo.mindiff);

  manopt(&optset, REQDBLOPT, 0, 'e', "minDMR2", 
      "minimal DMR 2", "<n>", NULL, &nfo.minDMR2);
  manopt(&optset, REQDBLOPT, 0, 'q', "mindiff2", 
      "minimal difference 2", "<n>", NULL, &nfo.mindiff2);

  manopt(&optset, REQUINTOPT, 0, 'l', "clustering", 
      "clustering or not: 0: no, 1: yes", "<n>", &clusteringconstraint, &nfo.clustering);

  manopt(&optset, REQUINTOPT, 0, 'p', "verbose", 
      "print segmenting position: 0: no, 1: yes", "<n>", &clusteringconstraint, &verbose);

  manopt(&optset, REQSTRINGOPT, 0, 'H', "header", 
      "header", "<string>", NULL, &headerfile);

  manopt(&optset, REQUINTOPT, 0, 'O', "outputImputed", 
      "output the matrix with imputed values: 0: no, 1: yes", "<n>", &clusteringconstraint, &nfo.outputImputed);


  args = manopt_getopts(&optset, argc, argv);
  if(args->noofvalues == 1) {
    manopt_help(&optset, "no source file provided.\n");
  }


  srand ((unsigned) nfo.randomseed);


  fi = initFileIterator(NULL, args->values[1]);

  // newcodes
  // /* check if groups are prefixes from one another */
  // if (strncmp(nfo.nameA, nfo.nameB, strlen(nfo.nameA)) == 0 ||
  //     strncmp(nfo.nameB, nfo.nameA, strlen(nfo.nameB)) == 0) {
  //   fprintf(stderr, "Error: name of group A (%s) must not be a prefix of name of group B (%s) "
  //           "and vice versa. Exit forced.\n", nfo.nameA, nfo.nameB);
  //   exit(-1);
  // }

  /* read header line with ids and assign to groups (prefix matches only) */
  ln = readcsvlines(NULL, fi, '\t', 1, &csv);

  headerfi = initFileIterator(NULL, headerfile);
  headerln = readcsvlines(NULL, headerfi, '\t', 1, &headercsv);
  for(k=2; k < headercsv[0]->noofstrings; k++) {
    char idstr[]="toBeReplaced";
    strcpy(idstr, headercsv[0]->strings[k].str);
    char *token = strtok(idstr, "_");
    int i = atoi(token);
    if (nfo.groups<(i+1))
    {
      nfo.groups=i+1;
    }
  }

  // subgroups: user-defined groups
  int subgroupNames_int[nfo.groups];char *subgroupNames = NULL;
  int *subgroupID[nfo.groups];
  int subgroupSize[nfo.groups];
  for (i = 0; i < nfo.groups; i++)
  {
    subgroupID[i] = NULL;
    subgroupSize[i] = 0;
    subgroupNames_int[i] = i;
  }
  for(k=2; k < headercsv[0]->noofstrings; k++) {
    char idstr[]="toBeReplaced";
    strcpy(idstr, headercsv[0]->strings[k].str);
    char *token = strtok(idstr, "_");
    int i = atoi(token);
    subgroupID[i] = ALLOCMEMORY(NULL, subgroupID[i], int, subgroupSize[i]+1);
    subgroupID[i][subgroupSize[i]] = k-2;
    subgroupSize[i]++;
  }

  if(verbose){
    for (i = 0; i < nfo.groups; i++)
    {
      fprintf(stderr, "Group: %d; Size: %d. The following ids belong to this group:\n", i, subgroupSize[i]);
      for (j = 0; j < subgroupSize[i]; j++)
      {
        fprintf(stderr, "%u: column: %d, name:%s\n", j, subgroupID[i][j], 
          headercsv[0]->strings[subgroupID[i][j]+2].str);
      }
    }
  }
  
  concatIntsToString(&subgroupNames, subgroupNames_int, nfo.groups, '|');
  // fprintf(stderr, "Single groups %s:\n", subgroupNames);

  // all combinations of subgroups
  // char *combinationNames = NULL;
  if(verbose){fprintf(stderr, "start combination.\n");}
  // int groupNumber = nfo.groups*(nfo.groups-1)/2 + nfo.groups; // #one vs one + #one vs others
  int **grpA_subgroups=NULL;
  int **grpB_subgroups=NULL;
  int groupNumber = calGroupNumber(nfo.groups, &grpA_subgroups, &grpB_subgroups, nfo.clustering);

  int ***groupID;
  int **groupSize;
  groupID = ALLOCMEMORY(NULL, NULL, int**, 2);
  groupSize = ALLOCMEMORY(NULL, NULL, int*, 2);
  for (i = 0; i < 2; i++)
  {
    groupID[i] = NULL;
    groupID[i] = ALLOCMEMORY(NULL, groupID[i], int*, groupNumber);
    groupSize[i] = NULL;
    groupSize[i] = ALLOCMEMORY(NULL, groupSize[i], int, groupNumber);
  }
  for (i = 0; i < groupNumber; i++)
  {
    groupID[0][i] = NULL;
    groupID[1][i] = NULL;
    groupSize[0][i] = 0;
    groupSize[1][i] = 0;
  }

  for (i = 0; i < groupNumber; i++)
  {
    int grpAindex = 0;
    int grpBindex = 0;
    
    // int combinationNames_int[nfo.groups-1];
    int j2=0;
    for (j = 0; j < nfo.groups; j++)
    {
      if (grpA_subgroups[i][j]==1)
      {
        // combinationNames_int[j2++]=j;
        groupID[0][i] = ALLOCMEMORY(NULL, groupID[0][i], int, groupSize[0][i]+subgroupSize[j]);
        for (k = 0; k < subgroupSize[j]; k++)
        {
          groupID[0][i][grpAindex] = subgroupID[j][k];
          grpAindex++;
        }
        groupSize[0][i] += subgroupSize[j];
      }

      if (grpB_subgroups[i][j]==1)
      {
        // combinationNames_int[j2++]=j;
        groupID[1][i] = ALLOCMEMORY(NULL, groupID[1][i], int, groupSize[1][i]+subgroupSize[j]);
        for (k = 0; k < subgroupSize[j]; k++)
        {
          groupID[1][i][grpBindex] = subgroupID[j][k];
          grpBindex++;
        }
        groupSize[1][i] += subgroupSize[j];
      }

    }

    char *tmp=NULL;
    // concatIntsToString(&tmp, combinationNames_int, nfo.groups-1, ',');
    // concatStrings(&combinationNames, tmp, '|');
    // fprintf(stderr, "CombinedGroup: %s\n", combinationNames);
    if ((nfo.clustering == 0)&&(verbose))
    {
      fprintf(stderr, "CombinedGroup: %d; Size: %d. The following ids belong to the first combined group:\n", i, groupSize[0][i]);
      for (k = 0; k < groupSize[0][i]; k++)
      {
        fprintf(stderr, "%u: column: %d, name:%s\n", j, groupID[0][i][k], 
          headercsv[0]->strings[groupID[0][i][k]+2].str);
      }
      fprintf(stderr, "CombinedGroup: %d; Size: %d. The following ids belong to the second combined group:\n", i, groupSize[1][i]);
      for (k = 0; k < groupSize[1][i]; k++)
      {
        fprintf(stderr, "%u: column: %d, name:%s\n", j, groupID[1][i][k], 
          headercsv[0]->strings[groupID[1][i][k]+2].str);
      }
    }

  }
  
  if(verbose){fprintf(stderr, "end combination. # of combination:%d\n", groupNumber);}
  // assert(0);

  destructStringset(NULL, csv[0]);
  FREEMEMORY(NULL, csv);
  destructStringset(NULL, headercsv[0]);
  FREEMEMORY(NULL, headercsv);
  
  /* init schedules and idle counter */
  idle = nfo.threads;
  schedule = ALLOCMEMORY(NULL, NULL, volatile char, nfo.threads);
  for(i=0; i < nfo.threads; i++) {
    schedule[i]=0;
  }
  
  threads = ALLOCMEMORY(space, NULL, pthread_t, nfo.threads);
  th_nfo = ALLOCMEMORY(space, NULL, metseg_t, nfo.threads);
  
  for(i=0; i < nfo.threads; i++) {
    memmove(&th_nfo[i], &nfo, sizeof(metseg_t));
  }

 
  
//###################### SINGLE CpG mode #########################
  if(nfo.mode == 3) {
      
      nfo.outputList = ALLOCMEMORY(NULL, NULL, list_out, 1);
      nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, 1);
      nfo.outputList->n=1000000;
      nfo.outputList->i=0;
      nfo.outputList->numberTests=0;
      nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, nfo.outputList->n);
    
      
      
        ln = readcsvlines(NULL, fi, '\t', 1, &csv);
        j = 0;
        while(ln) {
//check missing numbers            
            double *values = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);
            int nan = checkSetNAN(csv, values);
            if(nan>0) {
          //      fprintf(stderr,"call fillNAN");
                nan = fillNAN(values, subgroupID, subgroupSize, nfo.groups, &nfo);
         //       fprintf(stderr,"...done\n");
                
            } 
	    else {
		numnonnan++;
	    }
            if(nan>0) {
                destructStringset(NULL, csv[0]);
                csv[0] = NULL;
                FREEMEMORY(NULL, csv);
                csv = NULL;
                ln = readcsvlines(NULL, fi, '\t', 1, &csv);
                FREEMEMORY(NULL, values); 
                continue;
            } else {
              if (nfo.outputImputed==1){
                char *subtmp =NULL;
                outputImputedValues(&subtmp, values, csv[0]->noofstrings, '\t');
                fprintf(stdout,"//Imputed:");
                fprintf(stdout, my_strdup(csv[0]->strings[0].str));
                fprintf(stdout,"\t");
                fprintf(stdout, csv[0]->strings[1].str);
                fprintf(stdout,"\t");
                fprintf(stdout, subtmp);
                fprintf(stdout,"\n");
              }
            }

            cpg_t *cpg = ALLOCMEMORY(NULL, NULL, cpg_t, 1);
            double *groupA = ALLOCMEMORY(NULL, NULL, double, noA);
            double *groupB = ALLOCMEMORY(NULL, NULL, double, noB);
            
            for(i=0; i<noA; i++) { 
//                  groupA[i] = atof(csv[0]->strings[grpA[i]+2].str);
                    groupA[i] = values[grpA[i]+2];
            }
            for(i=0; i<noB; i++) { 
//                  groupB[i] = atof(csv[0]->strings[grpB[i]+2].str);
                    groupB[i] = values[grpB[i]+2];
            }
            cpg->groupA=groupA;
            cpg->groupB=groupB;
            cpg->chr = ALLOCMEMORY(NULL, NULL, char, csv[0]->strings[0].len+1);
            strcpy(cpg->chr, csv[0]->strings[0].str);
            cpg->noA=noA;
            cpg->noB=noB;
            cpg->start=atoi(csv[0]->strings[1].str);
            cpg->stop=atoi(csv[0]->strings[1].str);
            
            if(nfo.threads > 1) { 
                //wait for a free thread
                while(idle == 0);
                //look for the free slot
                // for(i=0; i < nfo.threads; i++) {
                //   if(schedule[i] == 0) break;
                // }
                int nofreeslot = 1;
                while (nofreeslot)
                {
                  for(i=0; i < nfo.threads; i++) {
                    if(schedule[i] == 0) {
                      nofreeslot = 0;
                      break;
                    }
                  }
                }
                //this must always hold because of idle variable
                assert(i < nfo.threads);
     //           fprintf(stderr, "starting thread %d\n", i);
                //decrement the idle and set the schedule
                pthread_mutex_lock(&cnt);
                idle--;
                schedule[i]=1;
                pthread_mutex_unlock(&cnt);
                //assign task
                th_nfo[i].threadno = i;
                th_nfo[i].cpg=cpg;
                th_nfo[i].outputList = nfo.outputList;
                fprintf(stderr, "CpG testing %s-[%d]\n", cpg->chr, cpg->start);
                //create the thread (detached!)
                pthread_create(&threads[i], &tattr, segworker_CpG, &th_nfo[i]);
                //now we must make sure that each thread keeps his own chunk
                //of the input data, thus the three arrays are simply set to NULL
                //the thread is going to take care of the deallocation
                cpg = NULL;

              } else { 
                fprintf(stderr, "CpG testing %s-[%d]\n", cpg->chr, cpg->start);
                int ua = mannwhitney (cpg->groupA, noA, cpg->groupB, noB);
		double p= mannwhitneyPvalue(ua, noA, noB, nfo.MWU, MAXM, MAXN);
		double ratio = get_meandiff(cpg, cpg->groupA, noA, cpg->groupB, noB);
            //    fprintf(stdout, "#%s\t%d\t%d\t.\t%.2g\t%f\n",cpg->chr, cpg->start,cpg->stop,ratio,p);
                
                
                cpgTest(cpg->chr, cpg->start,cpg->stop,ratio,p,&nfo, cpg->methA, cpg->methB);
                destructCpg(cpg);
              }
            
            
            
            destructStringset(NULL, csv[0]);
            csv[0] = NULL;
            FREEMEMORY(NULL, csv);
            csv = NULL;
            FREEMEMORY(NULL, values); 
           ln =readcsvlines(NULL, fi, '\t', 1, &csv);
    }
        
        
        

        
    
  }


  


  
  
  
  



  
//###################### DEFINED REGIONS mode#####################
  if(nfo.mode == 2) {
      fprintf(stderr, "Mode 2 -- pre-defined regions\n");
      nfo.outputList = ALLOCMEMORY(NULL, NULL, list_out, 1);
      nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, 1);
      nfo.outputList->n=1000000;
      nfo.outputList->i=0;
      nfo.outputList->numberTests=0;
      nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, nfo.outputList->n);
    
      
      
      
    //  exit(EXIT_SUCCESS);
//init segmentset of active regions
      bedfi = initFileIterator(NULL, bedfile);
      segmentset_t *set = ALLOCMEMORY(NULL, NULL, segmentset_t, 1);
      initSegmentSet(set);
//first region      
      bedln = readcsvlines(NULL, bedfi, '\t', 1, &bedcsv);
      set->nextchr = ALLOCMEMORY(NULL, NULL, char, bedcsv[0]->strings[0].len+1);
      strcpy(set->nextchr, bedcsv[0]->strings[0].str);
      set->nextstart = atoi(bedcsv[0]->strings[1].str)+1;
      
      set->chr = ALLOCMEMORY(NULL, NULL, char, bedcsv[0]->strings[0].len+1);
      strcpy(set->chr, bedcsv[0]->strings[0].str);
      
      ln = readcsvlines(NULL, fi, '\t', 1, &csv);
      int l=-1;
      while(ln) {
          l++;
//add missing values          
          double *values = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);
            int nan = checkSetNAN(csv, values);
            if(nan>0) {
  //              fprintf(stdout,"call fillNAN");
                nan = fillNAN(values, subgroupID, subgroupSize, nfo.groups, &nfo);
            }
	    else {
		numnonnan++;
	    }
            if(nan>0) {
                destructStringset(NULL, csv[0]);
                csv[0] = NULL;
                FREEMEMORY(NULL, csv);
                csv = NULL;
                ln = readcsvlines(NULL, fi, '\t', 1, &csv);
                FREEMEMORY(NULL, values); 
                continue;
            } else {
              if (nfo.outputImputed==1){
                char *subtmp =NULL;
                outputImputedValues(&subtmp, values, csv[0]->noofstrings, '\t');
                fprintf(stdout,"//Imputed:");
                fprintf(stdout, my_strdup(csv[0]->strings[0].str));
                fprintf(stdout,"\t");
                fprintf(stdout, csv[0]->strings[1].str);
                fprintf(stdout,"\t");
                fprintf(stdout, subtmp);
                fprintf(stdout,"\n");
              }
            }
         
          
          int pos = atoi(csv[0]->strings[1].str);
          
          
      //   fprintf(stdout,"########## %s %d (currChrom %s ,firststop %d, nextChr %s, nextStart %d)\n",csv[0]->strings[0].str,pos,set->chr,set->firststop, set->nextchr,set->nextstart);
          
//remove filled segments if
//Size of set >0 AND (currentChrNotNULL  OR  curr.Positon > FirstStopInSet OR  ChromosomeChangeForCpGInput)
          if(set->n>0 && (set->chr == NULL || pos > set->firststop || (strcmp(set->chr, csv[0]->strings[0].str) != 0))) {

   //           if((strcmp(set->chr, csv[0]->strings[0].str) != 0))
   //               fprintf(stdout,"CHROMCHANGE\n");
              
              segment_t *seg = set->head;
              set->firststop=-1;
              while(seg) {
                  segment_t *tmp = seg->next;
                  if(set->n>0 && ( set->chr == NULL || (strcmp(set->chr, csv[0]->strings[0].str) != 0) || seg->stop<pos)) {
       //                 fprintf(stdout,"@@@@@@@@@@Removing seg %s:%d-%d next%d parent%d\n",seg->chr,seg->start,seg->stop,seg->next == NULL,seg->parent == NULL);
                        removeThisSegmentFromSet(set,seg);
                        if(nfo.threads > 1) { 
                            //wait for a free thread
                            while(idle == 0);
                            //look for the free slot
                            // for(i=0; i < nfo.threads; i++) {
                            //   if(schedule[i] == 0) break;
                            // }
                            int nofreeslot = 1;
                            while (nofreeslot)
                            {
                              for(i=0; i < nfo.threads; i++) {
                                if(schedule[i] == 0) {
                                  nofreeslot = 0;
                                  break;
                                }
                              }
                            }
                            //this must always hold because of idle variable
                            assert(i < nfo.threads);
//                            fprintf(stderr, "starting thread %d\n", i);
                            //decrement the idle and set the schedule
                            pthread_mutex_lock(&cnt);
                            idle--;
                            schedule[i]=1;
                            pthread_mutex_unlock(&cnt);
                            //assign task
                            th_nfo[i].seg = seg;
                            th_nfo[i].groupID = groupID;
                            th_nfo[i].groupSize = groupSize;
                            th_nfo[i].groupNumber = groupNumber;
                            th_nfo[i].subgroupID = subgroupID;
                            th_nfo[i].subgroupSize = subgroupSize;
                            th_nfo[i].threadno = i;
                            th_nfo[i].outputList = nfo.outputList;
          
                            fprintf(stderr, "region testing %s-[%d,%d]\n", seg->chr, seg->start, seg->stop);
                            //create the thread (detached!)
                            pthread_create(&threads[i], &tattr, segworker_region, &th_nfo[i]);
                            //now we must make sure that each thread keeps his own chunk
                            //of the input data, thus the three arrays are simply set to NULL
                            //the thread is going to take care of the deallocation
                            seg = NULL;

                          } else { 
                            fprintf(stderr, "region testing %s-[%d,%d]\n", seg->chr, seg->start, seg->stop);
                            regionTest(seg, groupID, groupSize, groupNumber, subgroupID, subgroupSize, &nfo);
                              
                          }
                  }
                  else {
                      if(set->firststop == -1) { set->firststop = seg->stop; }
                      else      { set->firststop = MIN(seg->stop, set->firststop); }
                  }
                  seg = tmp;
              }
              if(set->n < 1) {
                  FREEMEMORY(NULL, set->chr); 
                  set->chr=NULL;
              }
          }
          
            if((!set->chr) && set->nextchr && (strcmp(set->nextchr, csv[0]->strings[0].str) != 0)) {
                 destructStringset(NULL, csv[0]);
                 csv[0] = NULL;
                 FREEMEMORY(NULL, csv);
                 csv = NULL;
                 ln = readcsvlines(NULL, fi, '\t', 1, &csv);
                 FREEMEMORY(NULL, values); 
                 continue;
             }
          
          
          
//add new segments 
       //   if(set->nextchr)
       //           fprintf(stdout,"while %d\t%d\t%d\n",bedln, pos>=set->nextstart ,(strcmp(set->chr,set->nextchr) == 0));
          while(bedln && pos>=set->nextstart && (!set->nextchr  || !set->chr || (strcmp(set->chr,set->nextchr) == 0)) ) {
              segment_t *tmp = addNewSegmentToSet(set);
              initSegment(tmp);
              tmp->chr = ALLOCMEMORY(NULL, NULL, char, bedcsv[0]->strings[0].len+1);
              strcpy(tmp->chr, bedcsv[0]->strings[0].str);
              tmp->start = atoi(bedcsv[0]->strings[1].str)+1;
              tmp->stop = atoi(bedcsv[0]->strings[2].str);
              
              if(tmp->stop < set->firststop) {
                  set->firststop = tmp->stop;
              }
              
              FREEMEMORY(NULL, set->chr); 
              if(set->n > 0) {
                  set->chr = ALLOCMEMORY(NULL, NULL, char, bedcsv[0]->strings[0].len+1);
                  strcpy(set->chr, bedcsv[0]->strings[0].str);
                  
                  
              }
              
              
         //read in next region
              destructStringset(NULL, bedcsv[0]);
              bedcsv[0] = NULL;
              FREEMEMORY(NULL, bedcsv);
              bedcsv = NULL;
              
              bedln = readcsvlines(NULL, bedfi, '\t', 1, &bedcsv); 

       //       if(bedln)
       //               fprintf(stdout,"seg %s:%s-%s\n",bedcsv[0]->strings[0].str,bedcsv[0]->strings[1].str,bedcsv[0]->strings[2].str);

              FREEMEMORY(NULL, set->nextchr);
              if(bedln) {
                set->nextchr = ALLOCMEMORY(NULL, NULL, char, bedcsv[0]->strings[0].len+1);
                strcpy(set->nextchr, bedcsv[0]->strings[0].str);
                set->nextstart = atoi(bedcsv[0]->strings[1].str)+1;
              }
              else {
                  set->nextstart = -1;
              }
              
          }
      
          
//add CpGs to regions          
        //add CpGs to regions
         
          segment_t *seg = set->head;
        //  if(seg && set->firststop==-1)
        //        { set->firststop = seg->stop; }
        //  int Notbreaking=1;
/*
          if(seg && seg->chr && csv[0] && (strcmp(seg->chr,csv[0]->strings[0].str) || ( (!strcmp(seg->chr,csv[0]->strings[0].str)) && seg->start > atoi(csv[0]->strings[1].str)))) {
  //          if(seg && seg->chr && csv[0] && (strcmp(seg->chr,csv[0]->strings[0].str))) {
              
              
              int z=0;
              while(csv[0] && seg && seg->chr && (strcmp(seg->chr,csv[0]->strings[0].str) || ( (!strcmp(seg->chr,csv[0]->strings[0].str)) && seg->start > atoi(csv[0]->strings[1].str)))){
    //          while(csv[0] && (strcmp(seg->chr,csv[0]->strings[0].str ))){
                  z++;
                destructStringset(NULL, csv[0]);
                csv[0] = NULL;
                FREEMEMORY(NULL, csv);
                csv = NULL;
                FREEMEMORY(NULL, values); 
                ln = readcsvlines(NULL, fi, '\t', 1, &csv);
              }
              fprintf(stdout,"removed %d lines\n",z);
              fprintf(stdout,"continue: segChr=%s csvChr=%s segStart=%d csvStart=%d\n",seg->chr,csv[0]->strings[0].str,seg->start,atoi(csv[0]->strings[1].str)+1);
              Notbreaking = 0;
          }
          */
        //  pos = atoi(csv[0]->strings[1].str)+1;
       //   fprintf(stdout,"########## %s %d (currChrom %s ,firststop %d, nextChr %s, nextStart %d)\n",csv[0]->strings[0].str,pos,set->chr,set->firststop, set->nextchr,set->nextstart);
          
          
          
          while(seg)  {
//              fprintf(stdout,"@@@@@@@@@@@@@@@@@@@@@adding seqs now %s %s\n",seg->chr,csv[0]->strings[0].str);
              if(!strcmp(seg->chr,csv[0]->strings[0].str)) {
//                fprintf(stderr,"#Adding CpG %s:%s to region %s:%d-%d\n",csv[0]->strings[0].str,csv[0]->strings[1].str,seg->chr,seg->start,seg->stop);
                seg->pos = ALLOCMEMORY(NULL, seg->pos, int,    seg->n+1); //index
                seg->value = ALLOCMEMORY(NULL, seg->value, double*, seg->n+1); //cpgs  

                seg->pos[seg->n] = atoi(csv[0]->strings[1].str);
                seg->value[seg->n] = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);
                for(k=2; k < csv[0]->noofstrings; k++) { 
//                  seg->value[seg->n][k-2] = atof(csv[0]->strings[k].str);
                    seg->value[seg->n][k-2] = values[k];
                }
                seg->n++;
              }
          seg = seg->next;    
          }
       //   if(Notbreaking==1){
            destructStringset(NULL, csv[0]);
            csv[0] = NULL;
            FREEMEMORY(NULL, csv);
            csv = NULL;
            FREEMEMORY(NULL, values); 
            ln = readcsvlines(NULL, fi, '\t', 1, &csv);
    //      }      
      
       }
      
//remaining regions to test 
      segment_t *seg = set->head;
      while(seg) {
                  segment_t *tmp = seg->next;
        //          fprintf(stderr,"@@@@@@@@@@Removing seg %s:%d-%d\n",seg->chr,seg->start,seg->stop);
                  if(nfo.threads > 1) { 
                            //wait for a free thread
                            while(idle == 0);
                            //look for the free slot
                            // for(i=0; i < nfo.threads; i++) {
                            //   if(schedule[i] == 0) break;
                            // }
                            int nofreeslot = 1;
                            while (nofreeslot)
                            {
                              for(i=0; i < nfo.threads; i++) {
                                if(schedule[i] == 0) {
                                  nofreeslot = 0;
                                  break;
                                }
                              }
                            }
                            //this must always hold because of idle variable
                            assert(i < nfo.threads);
                            fprintf(stderr, "starting thread %u\n", i);
                            //decrement the idle and set the schedule
                            pthread_mutex_lock(&cnt);
                            idle--;
                            schedule[i]=1;
                            pthread_mutex_unlock(&cnt);
                            //assign task
                            th_nfo[i].seg = seg;
                            th_nfo[i].grpA = grpA;
                            th_nfo[i].grpB = grpB;
                            th_nfo[i].noA = noA;
                            th_nfo[i].noB = noB;
                            th_nfo[i].threadno = i;
                            
                            fprintf(stderr, "region testing %s-[%d,%d]\n", seg->chr, seg->start, seg->stop);
                            //create the thread (detached!)
                            pthread_create(&threads[i], &tattr, segworker_region, &th_nfo[i]);
                            //now we must make sure that each thread keeps his own chunk
                            //of the input data, thus the three arrays are simply set to NULL
                            //the thread is going to take care of the deallocation
                            seg = NULL;

                          } else { 
    //                        fprintf(stderr, "region testing %s-[%d,%d]\n", seg->chr, seg->start, seg->stop);
                            regionTest(seg, groupID, groupSize, groupNumber, subgroupID, subgroupSize, &nfo);
                              
                          }
                  seg = tmp;
              }
    //  if(seg)
              destructSegmentSet(set);
              
          
              
  }




//###################### SEGMENTER (main) mode ###########################
  if(nfo.mode == 1) {
   //   fprintf(stderr,"#MODE2\n");
    nfo.outputList = ALLOCMEMORY(NULL, NULL, list_out, 1);
    nfo.outputList->n=1000000;
    nfo.outputList->i=0;
    nfo.outputList->numberTests=0;
    nfo.outputList->segment_out = ALLOCMEMORY(NULL, NULL, segment_out, nfo.outputList->n);
    //fprintf(stderr,"output->n: %d\n",nfo.outputList->n);
      
      
    ln = readcsvlines(NULL, fi, '\t', 1, &csv); // zzhu$ reading the methyl table into variable 'csv'
    j = 0;
    while(ln) { 
        //fprintf(stderr,"#new LINE\n");
//check missing numbers            
        double *values = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);
        int nan = checkSetNAN(csv, values);
        if(nan>0) {
        // fprintf(stderr,"#call fillNAN for %d groups\n", nfo.groups);
            nan = fillNAN(values, subgroupID, subgroupSize, nfo.groups, &nfo);
    //      fprintf(stderr,"#...done\n");
        }
	else{
	    numnonnan++;
	}
     //   fprintf(stdout,"#LINES INPUT\n");
        if(nan>0) {
     //       fprintf(stdout,"#REMOVING LINE\n");
            destructStringset(NULL, csv[0]);
            csv[0] = NULL;
            FREEMEMORY(NULL, csv);
            csv = NULL;
            ln = readcsvlines(NULL, fi, '\t', 1, &csv);
            FREEMEMORY(NULL, values); 
            continue;
        }
        else {
    //            fprintf(stdout,"#LINE OKAY \n");
              if (nfo.outputImputed==1){
                char *subtmp =NULL;
                outputImputedValues(&subtmp, values, csv[0]->noofstrings, '\t');
                fprintf(stdout,"//Imputed:");
                fprintf(stdout, my_strdup(csv[0]->strings[0].str));
                fprintf(stdout,"\t");
                fprintf(stdout, csv[0]->strings[1].str);
                fprintf(stdout,"\t");
                fprintf(stdout, subtmp);
                fprintf(stdout,"\n");
              }
        }
      char *x = my_strdup(csv[0]->strings[0].str); //zzhu$ x: chromosome in current line
      int y = atoi(csv[0]->strings[1].str); //zzhu$ y: CpG position in current line

      if(j > 0 && (strcmp(x, chr[j-1]) || y > pos[j-1] + nfo.maxdist ||
                   (nfo.maxseg > 0 && j >= nfo.maxseg))) {

        if(nfo.threads > 1) { 
          //wait for a free thread
          while(idle == 0);
          //look for the free slot
          int nofreeslot = 1;
          while (nofreeslot)
          {
            for(i=0; i < nfo.threads; i++) {
              if(schedule[i] == 0) {
                nofreeslot = 0;
                break;
              }
            }
          }
          

          //this must always hold because of idle variable
          assert(i < nfo.threads);
  //        fprintf(stderr, "starting thread %d\n", i);
          //decrement the idle and set the schedule
          pthread_mutex_lock(&cnt);
          idle--;
          schedule[i]=1;
          pthread_mutex_unlock(&cnt);
          //assign task
  /*copied from segworker*/
  //         segmentation(t->chr, t->pos, t->value, t->n, t->groupID, t->groupSize, t->groupNumber, t->subgroupID, t->subgroupSize, t);
  
  // //cleanup own data
  // for(i=0; i < t->n; i++) {
  //   FREEMEMORY(NULL, t->chr[i]);
  //   FREEMEMORY(NULL, t->value[i]);
  // }
  // FREEMEMORY(NULL, t->chr);
  // FREEMEMORY(NULL, t->pos);
  // FREEMEMORY(NULL, t->value);
          th_nfo[i].chr = chr;
          th_nfo[i].pos = pos;
          th_nfo[i].value = val;
          th_nfo[i].n = j;
          th_nfo[i].groupID = groupID;
          th_nfo[i].groupSize = groupSize;
          th_nfo[i].groupNumber = groupNumber;
          th_nfo[i].subgroupID = subgroupID;
          th_nfo[i].subgroupSize = subgroupSize;
          th_nfo[i].threadno = i;
          th_nfo[i].outputList = nfo.outputList;
          
          
          if(verbose){fprintf(stderr, "Thread: %d segmenting %s-[%d,%d], %u CpGs\n", i, chr[0], pos[0], pos[j-1], j);}
          //create the thread (detached!)
          pthread_create(&threads[i], &tattr, segworker, &th_nfo[i]);
          //now we must make sure that each thread keeps his own chunk
          //of the input data, thus the three arrays are simply set to NULL
          //the thread is going to take care of the deallocation
          // fprintf(stderr, "segmented %s-[%d,%d], %u CpGs\n", chr[0], pos[0], pos[j-1], j);

          chr = NULL;
          pos = NULL;
          val = NULL;

        } else { 
          if(verbose){fprintf(stderr, "Segmenting %s-[%d,%d], %u CpGs\n", chr[0], pos[0], pos[j-1],j);}
          // segmentation(chr, pos, val, j, grpA, noA, grpB, noB, &nfo);
          segmentation(chr, pos, val, j, groupID, groupSize, groupNumber, subgroupID, subgroupSize, &nfo);
          for(i=0; i < j; i++) { 
            FREEMEMORY(NULL, chr[i]);
            FREEMEMORY(NULL, val[i]);
          }   
          // if(nfo.outputList->i>=2){fprintf(stderr,"start7:%d\n",nfo.outputList->segment_out[2].start);}
        }

        j = 0;
      } 

      chr = ALLOCMEMORY(NULL, chr, char*,   j+1); //chr
      pos = ALLOCMEMORY(NULL, pos, int,    j+1); //index
      val = ALLOCMEMORY(NULL, val, double*, j+1); //cpgs  

      chr[j] = x;
      pos[j] = y;
      val[j] = ALLOCMEMORY(NULL, NULL, double, csv[0]->noofstrings);

      for(k=2; k < csv[0]->noofstrings; k++) { 
//        val[j][k-2] = atof(csv[0]->strings[k].str);
        val[j][k-2] = values[k];
      }


      j+=1; // zzhu$ j: the number of CpGs in the segment
      destructStringset(NULL, csv[0]);
      csv[0] = NULL;
      FREEMEMORY(NULL, csv);
      csv = NULL;
      FREEMEMORY(NULL, values); 
      ln =readcsvlines(NULL, fi, '\t', 1, &csv);
    } 
  
    if(verbose){fprintf(stderr, "segmenting %s-[%d,%d], %u CpGs \n", chr[0], pos[0], pos[j-1],j);}
    segmentation(chr, pos, val, j, groupID, groupSize, groupNumber, subgroupID, subgroupSize, &nfo);
    for(i=0; i < j; i++) { 
        FREEMEMORY(NULL, chr[i]);
        FREEMEMORY(NULL, val[i]);
    }
    // if(nfo.outputList->i>=2){fprintf(stderr,"start8:%d\n",nfo.outputList->segment_out[2].start);}
  }
  
  
  
  
  
  
  
  
  
  
  //wait for all threads to terminate
  while(idle != nfo.threads);

  if(nfo.mode == 1 || nfo.mode == 2) {
    if(verbose){fprintf(stderr, "Number of Tests: %d\n", nfo.outputList->numberTests);}
    // if(nfo.outputList->i>=2){fprintf(stderr,"start91:%d\n",nfo.outputList->segment_out[2].start);}
    multiple_testing_correction(nfo.outputList, nfo.mode, nfo.mtc);
    // if(nfo.outputList->i>=2){fprintf(stderr,"start92:%d\n",nfo.outputList->segment_out[2].start);}
    if(verbose){fprintf(stderr, "Multiple testing correction done.\n");}
    fprintf(stderr, "Number of non-NA CpGs: %d\n", numnonnan);
    // fprintf(stdout, "chr\tstart\tstop\tq\tmeandiff\tlength\tmwu\tp\t%s\tsig.comparison\n",subgroupNames);
    fprintf(stdout, "chr\tstart\tstop\tq\tmeandiff\tlength\tmwu\tp\tmean\tsig.comparison\n",subgroupNames);
    for(int i=0;i<nfo.outputList->i;i++){
      // fprintf(stderr,"start10:%d\n",nfo.outputList->segment_out[2].start);
      // fprintf(stderr, "TEST %d: %d,%f.\n",i,nfo.outputList->segment_out[i].start,nfo.outputList->segment_out[i].meandiff);

      if (nfo.clustering==1)
      {
        if(nfo.outputList->segment_out[i].meandiff >= nfo.minMethDist || nfo.outputList->segment_out[i].meandiff <= -1* nfo.minMethDist) {
          // fprintf(stderr, "TEST %d: %d,%d.\n",i,nfo.outputList->segment_out[i].start,nfo.outputList->segment_out[i].stop);
          fprintf(stdout, "%s\t%d\t%d\t%.5g\t%f\t%d\t%.5g\t%.5g\t%s\t%s\n", 
                  nfo.outputList->segment_out[i].chr,
                  nfo.outputList->segment_out[i].start,
                  nfo.outputList->segment_out[i].stop,
                  nfo.outputList->segment_out[i].q,
                  nfo.outputList->segment_out[i].meandiff,
                  nfo.outputList->segment_out[i].length,                
                  nfo.outputList->segment_out[i].mwu,
                  nfo.outputList->segment_out[i].p,
                  nfo.outputList->segment_out[i].methA,
                  // nfo.outputList->segment_out[i].methB,
                  nfo.outputList->segment_out[i].methB);
        }
      } else {
        if(nfo.outputList->segment_out[i].meandiff >= nfo.minMethDist || nfo.outputList->segment_out[i].meandiff <= -1* nfo.minMethDist) {
          // fprintf(stderr, "TEST %d: %d,%d.\n",i,nfo.outputList->segment_out[i].start,nfo.outputList->segment_out[i].stop);
          fprintf(stdout, "%s\t%d\t%d\t%.5g\t%f\t%d\t%.5g\t%.5g\t%s\t%f\n", 
                  nfo.outputList->segment_out[i].chr,
                  nfo.outputList->segment_out[i].start,
                  nfo.outputList->segment_out[i].stop,
                  nfo.outputList->segment_out[i].q,
                  nfo.outputList->segment_out[i].meandiff,
                  nfo.outputList->segment_out[i].length,                
                  nfo.outputList->segment_out[i].mwu,
                  nfo.outputList->segment_out[i].p,
                  nfo.outputList->segment_out[i].methA,
                  // nfo.outputList->segment_out[i].methB,
                  nfo.outputList->segment_out[i].sigcp);
        }
      }
      
      
    }
  }
    
  if(nfo.mode == 3) {
    if(verbose){fprintf(stderr, "Number of Tests: %d\n", nfo.outputList->numberTests);}
    multiple_testing_correction(nfo.outputList, nfo.mode, nfo.mtc);
    if(verbose){fprintf(stderr, "Multiple testing correction done.\n");}
    for(int i=0;i<nfo.outputList->i;i++){
      if(nfo.outputList->segment_out[i].meandiff >= nfo.minMethDist || nfo.outputList->segment_out[i].meandiff <= -1* nfo.minMethDist) {
        fprintf(stdout, "%s\t%d\t%d\t%.5g\t%f\t%d\t%.5g\t%s\n", 
                nfo.outputList->segment_out[i].chr,
                nfo.outputList->segment_out[i].start,
                nfo.outputList->segment_out[i].stop,
                nfo.outputList->segment_out[i].q,
                nfo.outputList->segment_out[i].meandiff,
                nfo.outputList->segment_out[i].length,
                nfo.outputList->segment_out[i].mwu,
                nfo.outputList->segment_out[i].methA
                // nfo.outputList->segment_out[i].methB
                );
        }
    }
  }
 
  fflush(stdout); 
  
  
  FREEMEMORY(NULL, chr);
  FREEMEMORY(NULL, pos);
  FREEMEMORY(NULL, val);
  FREEMEMORY(NULL, grpA);
  FREEMEMORY(NULL, grpB);
  
  closeFileIterator(NULL, fi);
  FREEMEMORY(NULL, fi);

  if(csv && csv[0]) destructStringset(NULL, csv[0]);
  if(csv) FREEMEMORY(NULL, csv);
  for(int i=0;i<nfo.outputList->i;i++){
    FREEMEMORY(NULL, nfo.outputList->segment_out[i].chr);
  }
  FREEMEMORY(NULL, nfo.outputList->segment_out);
  FREEMEMORY(NULL, nfo.outputList);

  pthread_attr_destroy(&tattr);  
  FREEMEMORY(NULL, threads);
  FREEMEMORY(NULL, th_nfo);
  free((void*)schedule);

  destructMannWhitneyCDFMatrix(nfo.MWU, MAXM, MAXN);
  manopt_destructoptionset(&optset);
  manopt_destructarg(args);
  FREEMEMORY(NULL, args);

  exit(EXIT_SUCCESS);
}

