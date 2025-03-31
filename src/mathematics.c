
/*
 * mathematics.c
 * implemtation of various mathematical functions
 *
 * @author Steve Hoffmann
 * @date Wed 22 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 54 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-10 22:13:30 +0200 (Wed, 10 Sep 2008) $
 *
 *  Id: $Id: mathematics.c 54 2008-09-10 20:13:30Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/trunk/libs/mathematics.c $
 *  
 */


#include "mathematics.h"
#include "sort.h"
#include <float.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

 void *initArray(void *space, int size, size_t datatype) {
	void *ptr=NULL;

	/*dirty trick: sizeof(char) == 1*/
	ptr = ALLOCMEMORY(space, ptr, char, size*datatype);
	return ptr;
 }



/*---------------------------------- coldel ----------------------------------
 *    
 * @brief delete column for matrix
 * @author Steve Hoffmann 
 *   
 */
 
double*
coldel (void *space, double *a, Uint m, Uint n, Uint d) {
	
	double *t;
	Uint	i,
			j=-1,
			k=0,
			l=0;

  t = (double*) INITMATRIX2D(space, m, (n-1), sizeof(double));

  for(i=0; i < m*n; i++) {
	if(i % n == 0) { 
	  j++; k=0; l=0;
	} 	
	if(k++ != d) {
	  MATRIX2D(t, n-1, j, l++) = a[i];
	}
  }
	
  FREEMEMORY(space, a);
  return t;
}


/*---------------------------------- rowdel ----------------------------------
 *    
 * @brief delete row from matrix
 * @author Steve Hoffmann 
 *   
 */
 
double*
rowdel (void *space, double *a, Uint m, Uint n, Uint d) {
	
	double *t;
	Uint	i,
			j=-1,
			k=0,
			l=-1;

  t = (double*) INITMATRIX2D(space, (n-1), m, sizeof(double));

  for(i=0; i < m*n; i++) {
	if(i % n == 0) { 
	  j++; k=0;
	  l = (j != d) ? l+1 : l;
	} 	
	if(j != d) {
	  MATRIX2D(t, n, l, k++) = a[i];
	}
  }

  FREEMEMORY(space, a);
  return t;
}

/*----------------------------------- add ------------------------------------
 *    
 * @brief componentwise addition of a to a vector of length m
 * @author Steve Hoffmann 
 *   
 */

double*
add(double *x, Uint m, double a) {
  Uint i;

  for(i=0; i < m; i++) {
    x[i] += a;
    //fprintf(stdout, "add: %f -> %f\n", x[i]-a, x[i]);
  }
  return x;
}


/*----------------------------------- mean -----------------------------------
 *    
 * @brief calculate the arithmetic mean for a vector of length m
 * @author Steve Hoffmann 
 *   
 */

double
mean (double *x, Uint m) {
  Uint i;
  double sum=0;

  for (i=0; i < m; i++) {
    sum += x[i];
  }
  
  return sum/m;
}

/*---------------------------------- scalar ----------------------------------
 *    
 * @brief calculate the scalar product of two vectors of length m
 * @author Steve Hoffmann 
 *   
 */

double
scalar (double* x, double *y, Uint m) {
  double  p=0;
  Uint 	i;

  for (i=0; i < m; i++) {
    //fprintf(stdout, "scal: %f*%f + %f\n", x[i], y[i], p);
    p += x[i]*y[i];
  }
  
  return p;
}


/*----------------------------------- cov ------------------------------------
 *    
 * @brief get the covariance matrix (2x2) for two vectors of length m
 * @author Steve Hoffmann 
 *   
 */

double*
cov (void *space, double *x, double *y, Uint m) {
  double *c,
         xm,
         ym;

  c = (double*) INITMATRIX2D(space, 2, 2, sizeof(double));
  xm = mean(x, m);
  ym = mean(y, m);

  //fprintf(stdout, "meanx: %f\n", xm);
  //fprintf(stdout, "meany: %f\n", ym);

  /*center*/
  add(x, m, (-1)*xm);
  add(y, m, (-1)*ym);

  MATRIX2D(c, 2, 0, 0) = (double) scalar(x,x,m)/(m-1);
  MATRIX2D(c, 2, 0, 1) = MATRIX2D(c, 2, 1, 0) = (double) scalar(x,y,m)/(m-1);
  MATRIX2D(c, 2, 1, 1) = (double) scalar(y,y,m)/(m-1);

  return c;
}

/*----------------------------------- var ------------------------------------
 *    
 * @brief get the variance
 * @author Steve Hoffmann 
 *   
 */
 
double
var (double *x, Uint n)
{   
    int i;
    double m, r, sum=0;

    m=mean(x, n);
    for (i=0; i < n; i++) {
      r = x[i]-m;
      sum += (r*r);
    }

	return sum/n;
}


/*----------------------------------- rho ------------------------------------
 *    
 * @brief calculate correlation $\rho$ for two vectors of length m
 * @author Steve Hoffmann 
 *   
 */

double
rho (void *space, double *x, double *y, Uint m) {
  double *cv;
  double ret;

  cv = cov(space, x, y, m); 
  ret = (MATRIX2D(cv, 2, 0, 1)/sqrt(MATRIX2D(cv, 2, 0, 0)*MATRIX2D(cv, 2, 1, 1)));
  //fprintf(stdout, "%f %f %f -> %f\n",MATRIX2D(cv, 2, 0, 1), MATRIX2D(cv, 2, 0, 0), MATRIX2D(cv, 2, 1, 1), ret);
  FREEMEMORY(NULL, cv);

  return ret;
}


/*-------------------------------- lancozs gamma -------------------------------
 *    
 * @brief lancozs approximation for the gamma function G=5 N=6+1
 * @author Steve Hoffmann 
 *   
 */
 
double
gammaln (double x)
{
    Uint i;
    double y,tmp,sum;
    double base;
    double G = 5.0;
    int N = 7;
    
    double lancozscoef[7] = { 1.000000000190015, 
                             76.18009172947146, 
                            -86.50532032941677, 
                             24.01409824083091, 
                             -1.231739572450155,  
                              0.1208650973866179e-2, 
                             -0.5395239384953e-5}; 
    y = x;
    
    sum = lancozscoef[0];
    for(i=1; i < N; i++) {
      sum += lancozscoef[i]/++y;
    }

    base = x + G + 0.5;
    tmp = (base) - log(base)*(x+0.5);
    return -tmp+log(M_SQRT2PI*sum/x);
}


/*-------------------------------- kscdf ---------------------------------
 *    
 * @brief kscdf returns 1-x of cumulative distribution function of 
 * K distribution
 * @author Steve Hoffmann 
 *   
 */
 
double kscdf (double x)
{

  unsigned int k;
  double sum = 0.0, tmp=0.0, old=0.0, val;
  double base, coeff = 1.0;

  base = -2.0*(x*x);

  for(k=1; k <= 100; k++) {
    tmp = exp(base*k*k);
    sum += coeff*tmp;
    if (tmp <= EPSILON3*old || tmp <= EPSILON8*sum) { 
      val = (2.0*sum);   
      return val;
    }
    coeff *= -1;
    old = tmp;
  }

  return 1.0;
}


/*------------------------------ IsFiniteNumber ------------------------------
 *    
 * @brief check if x is a finite number
 * @author Steve Hoffmann 
 *   
 */
 

char 
IsFiniteNumber(double x) {
  return (x <= DBL_MAX && x >= -DBL_MAX); 
}    


double 
mannwhitneyGamma(double a, double x, double eps, int iter) {
  double an = 1.0 / a;
  double psum = an;
  double n = 0.0;
  while ( n < iter && fabs(an/psum) > eps) {
    n += 1;
    an *= x/(a+n);
    psum += an;
  }
  n = psum * exp((a * log(x)) - gammaln(a) - x );
  return n;
}

double
mannwhitneyPvalue(Uint u, Uint m, Uint n, double ***CDF, Uint maxm, Uint maxn) {

  double e, z, x, p;
  double u1;

  if(CDF != NULL && m <= maxm && n <= maxn) {
    if(m<n){
      p = CDF[m][n][u];
    }
    else {
      p = CDF[n][m][u];
    }
    p=1-p;
  } else { 
    u1 = u;
    z = (u1 - (m * n / 2.0)) / sqrt(m * n * (m + n + 1) / 12.0);
    x = z / sqrt(2.0);

    if (fabs(x) > 40) {
      if(x > 0) {e=1;} else {e=-1;}
    }
    else {
      e = mannwhitneyGamma(0.5, x * x, 1.0e-15, 10000);
      if(x < 0) {e=-e;}
    }
    if (fabs(z) > 40) {
      if(z<0) p=0.0; else p = 0.5;
    }
    else {
      p = 0.5 * (1 + e);
    }
//     p = 1-p;
  }

//TODO: please check one-sided and two-sided test with p*=2 there are p-values > 1!
    p*=2;

//to be safe ;D     
    p = MIN(p,1.0);
  return p;
}


/*------------------------------- mannwhitney --------------------------------
 *    
 * @brief mann whitney u test
 * @author Steve Hoffmann 
 *   
 */
 
Uint
mannwhitney (double *a, Uint m, double *b, Uint n)
{
    
    Uint i, *sorted;
    double *ranks, md, nd,sum1=0.0,u1, u2;
    double *c;
    c = ALLOCMEMORY(NULL, NULL, double, m+n);
    ranks = ALLOCMEMORY(NULL, NULL, double, m+n);
    memmove(c, a, m*sizeof(double));
    memmove(&c[m], b, n*sizeof(double));
    sorted = quickSort(NULL, c, m+n, cmp_dbl, NULL);
    
    Uint j=-1;
    double entry=NAN;
    for(int i=0;i<m+n;i++) {
        double rank = (double) i +1;
        if(c[sorted[i]] != entry) {
            if(j<i-1){
                double m=0;
                if(j>0) { 
                    m=j-1; 
                    m=(j*j+j)/2;
                }
                double oldrank = i;
                oldrank=(oldrank*oldrank+oldrank)/2;
                oldrank-=m;
                oldrank/=i-j;
                for(int z=j;z<i;z++)
                    ranks[sorted[z]] = oldrank;
            }
            j=i;
            ranks[sorted[i]] = rank;
        }
        entry = c[sorted[i]];
    }
    if(j<n+m-1){
        int i=n+m;
            double m=0;
            if(j>0) { 
                m=j-1; 
                m=(j*j+j)/2;
            }
            double oldrank = i;
            oldrank=(oldrank*oldrank+oldrank)/2;
            oldrank-=m;
            oldrank/=i-j;
            for(int z=j;z<i;z++)
                ranks[sorted[z]] = oldrank;
        }
/*      
    fprintf(stdout,"double X[%d] = {%f",m,(double) a[0]);
    for(int i=1;i<m;i++)
        fprintf(stdout,",%f",(double) a[i]);
    fprintf(stdout,"}\n");
    
    fprintf(stdout,"ranks X[%d] = {%f",m,ranks[0]);
    for(i=1; i <  m; i++)
        fprintf(stdout,",%f",ranks[i]);
    fprintf(stdout,"}\n");
    
    
    fprintf(stdout,"double Y[%d] = {%f",n,(double) b[0]);
    for(int i=1;i<n;i++)
        fprintf(stdout,", %f",(double) b[i]);
    fprintf(stdout,"}\n");
  
    
    
    fprintf(stdout,"original %f",(double) c[0]);
    for(int i=1;i<m+n;i++)
        fprintf(stdout,"\t%f",(double) c[i]);
    fprintf(stdout,"\n");
    
    fprintf(stdout,"ranks %f",(double) ranks[0]);
    for(int i=1;i<m+n;i++)
        fprintf(stdout,"\t%f",(double) ranks[i]);
    fprintf(stdout,"\n");
    
    fprintf(stdout,"sorted values %f",(double) c[sorted[0]]);
    for(int i=1;i<m+n;i++)
        fprintf(stdout,"\t%f",(double) c[sorted[i]]);
    fprintf(stdout,"\n");
    
    fprintf(stdout,"sorted ranks \n%f",(double) ranks[sorted[0]]);
    for(int i=1;i<m+n;i++)
        fprintf(stdout,"\t%f",(double) ranks[sorted[i]]);
    fprintf(stdout,"\n");
  
    */
    
    md = (double) m;
    nd = (double) n;
    //fprintf(stdout,"s\t%f\n",sum1);
    for(i=0; i <  m; i++) {
      sum1 +=   ranks[i];
  //    fprintf(stdout,"s\t%f\n",sum1);
    }
  //  fprintf(stdout,"sum\t%f\n",sum1);
    
    u1 = md*nd + (md*(md+1.0))/2.0 - sum1;
    u2 = md*nd - u1;
  //  fprintf(stdout,"Us: %d\t%d\n",u1,u2);
//    u2 = MIN(u1,u2);
    u1 = MAX(u1,u2);
    u2 = md*nd - u1;
    
    FREEMEMORY(NULL, c);
    FREEMEMORY(NULL, ranks);
    FREEMEMORY(NULL, sorted);
 
    return (int)u2;
}


/*-------------------------- generateMannWhitneyCDF --------------------------
 *    
 * @brief generate the U-Statistics CDF for n,m
 * @author Steve Hoffmann 
 *   
 */
 
double*
generateMannWhitneyCDF(Uint m, Uint n) {
  double *C;
  long long unsigned int *S;
//  long long unsigned int minU, maxU, i, j, k, iter=0, u=0;
  long long unsigned int maxU, i, j, k, iter=0, u=0;

  assert(m <=n );
  maxU = m*n + (m*(m+1))/2;

//replace Umax with Umin, from here maxU <- Umin
//  minU = m*n - maxU;
//  maxU = MIN(maxU,minU);
  
  //setup matrix
  C = ALLOCMEMORY(NULL, NULL, double, maxU);
  S = ALLOCMEMORY(NULL, NULL, long long unsigned int, m);
  memset(C, 0, sizeof(double)*maxU);
  C[0] = 1;
  
  for(i=0; i < m; i++) {
    S[i]=i+1;
  }

  i =0;
  j = m-1;

  while(1) {

    if(iter > 0) {
      Uint maxrank = m+n;
      j = m-1;

      while(j > 0 && S[j] == maxrank) { 
        maxrank--;
        j--;
      }
      //the highest pos has its max rank (=m+1) 
      //all U are calculated
      if(j==0 && S[j] == maxrank) 
        break;

      k = m-j-1;
      j = m-1;
      maxrank = m+n;

      while(j > 0 && S[j]==maxrank) {
        //this is the lowest element possible
        S[j] = S[j-k]+k+1;
        j--; 
        k--;
        maxrank--;
      }
      S[j]++; 
    }

    for(u=0, i=0; i < m; i++) {
      u += S[i];
    }

    u = maxU - u;
    C[u]++;
    iter++;
  }
    
  C[maxU-1] = C[maxU-1]/iter;
  
  for(i=maxU-2; i > 0; i--) {
    C[i] = C[i+1] + C[i]/iter;
  }
  
  C[0] = C[1];

  FREEMEMORY(NULL, S);

  return C;
}


/*----------------------- generateMannWhitneyCDFMatrix -----------------------
 *    
 * @brief generate the CDF matrix
 * @author Steve Hoffmann 
 *   
 */
 
double***
generateMannWhitneyCDFMatrix(Uint maxm, Uint maxn)
{
  double ***M;
  Uint i, j;

  M = ALLOCMEMORY(NULL, NULL, double**, maxm+1);
  
  for(i=0; i <= maxm; i++) {
    M[i] = ALLOCMEMORY(NULL, NULL, double*, maxn+1);
    for(j=i; j <= maxn; j++) {
      if(i >= 1) { 
	//      fprintf(stderr, "generating matrix for %d,%d\n", i, j);
      M[i][j] = generateMannWhitneyCDF(i,j);
      }
    }
  }

  return M;
}


/*----------------------- destructMannWhitneyCDFMatrix -----------------------
 *    
 * @brief destruct the MWU matrix
 * @author Steve Hoffmann 
 *   
 */
 
void
destructMannWhitneyCDFMatrix ( double ***CDF, Uint m, Uint n)
{
  Uint i,j;

  for(i=0; i <= m; i++){
    for(j=i; j <= n; j++){
      if(i >= 1)
      FREEMEMORY(NULL, CDF[i][j]);
    }
    FREEMEMORY(NULL, CDF[i]);
  }
	
  FREEMEMORY(NULL, CDF);
  return ;
}

double uniform_rand_range(double rangeLow, double rangeHigh) {

    double myRand = rand()/(1.0 + RAND_MAX);
    double range = rangeHigh - rangeLow ;
    double myRand_scaled = (myRand * range) + rangeLow;
    return myRand_scaled;
}

double uniform_rand() {
    return uniform_rand_range(0.0,1.0);
}

void randgauss(double *r1, double *r2)
{

  double u1, u2, r, f;

  do {
//    u1 = 2.0 * ((double)rand() / (double)RAND_MAX) -1.0;
//    u2 = 2.0 * ((double)rand() / (double)RAND_MAX) -1.0;
    u1 = 2.0 * uniform_rand()  -1.0;
    u2 = 2.0 * uniform_rand()  -1.0;
    
    r = u1*u1 + u2*u2;
  } while(r >= 1.0 || r ==0.0);

  f = sqrt(-2.0*log(r)/r);
  *r1 = f*u1;
  *r2 = f*u2;

  return ;
}

/*--------------------------------- randgam ----------------------------------
 *    
 * @brief generate random gamma
 * @author Steve Hoffmann
 * http://www.hongliangjie.com/2012/12/19/how-to-generate-gamma-random-variables/
 *   
 */

double randgam (double alpha, double beta)
{

  double x, z1=0, z2=0, u, v, d, c;
  char phase = 0;
 
  assert(alpha > 0);

  if(alpha >= 1.0) { // ? > or >=
    
    d = alpha-1.0/3.0;
    c = 1.0/sqrt(9.0*d);
    while(1) { 

      if(phase == 0) {  
        randgauss(&z1, &z2);
        phase = 1;
      } else {
        z1 = z2;
        phase = 0;
      }

      if (z1 > -1.0/c) {
        v = (1+c*z1);
        v = v*v*v;
//        u = (double)rand() / (double)RAND_MAX; // U(0,1)
        u = uniform_rand();
        if(log(u) < 0.5*(z1*z1)+d-d*v+d*log(v)) break;
      } 
    }
    x = d*v/beta;
  
  } else {
    
    x = randgam(alpha+1, beta);
//    u = (double)rand() / (double)RAND_MAX; // U(0,1)
    u = uniform_rand();
    
    x = x * pow(u, 1.0/beta);
  }

  return x;
}


/*--------------------------------- randbeta ---------------------------------
 *    
 * @brief random beta
 * @author Steve Hoffmann 
 *   
 */
 
double randbeta (double alpha, double beta)
{
  
  double x,y,u,v;

  if (alpha <= 1.0 && beta <= 1.0){
    while(1) {
//      u= (double)rand() / (double)RAND_MAX;
//      v= (double)rand() / (double)RAND_MAX;
      u= uniform_rand();
      v= uniform_rand();
      x = pow(u, 1.0/alpha);
      y = pow(v, 1.0/beta);
      if(x+y <= 1.0) {
        return x/(x+y);
      }
    }
  } else {
    x = randgam(alpha, 1);
    y = randgam(beta,  1);
    return x/(x+y);
  }
}

double rbeta_mv(double mu, double var) {
    if(mu <= 0.000001)
        mu=0.000001;
    if(var <= 0.000001)
        //var=0.0000000001;
    //if( var <= 0.0000000001)
        return mu;
//  fprintf(stdout,"no return\n");  
  double x = ((mu * (1-mu))/var)-1;
  double alpha = x * mu;
  double beta = (1-mu) * x;
//  fprintf(stderr,"#alpha: %f\tbeta: %f\n",alpha,beta);
    
  if(alpha <= 0.000001)
        alpha=0.000001;
    if(beta <= 0.000001)
        beta=0.000001;
 // fprintf(stderr,"#alpha: %f\tbeta: %f\n",alpha,beta);
  
  //double alpha = ((1 - mu) / var - 1 / mu) * (mu *mu);
  //double beta = alpha * (1 / mu - 1);
  //fprintf(stdout,"MuVar\t%f\t%f\n",mu,var);
 // fprintf(stdout,"AlphaBeta\t%f\t%f\n",alpha,beta);
//  double ret = rbeta(alpha, beta);
  double ret = randbeta(alpha, beta);
  
  return ret;
}






