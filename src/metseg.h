#ifndef METSEG_H
#define METSEG_H
/*
 *
 *	metseg.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 09/10/2014 01:46:32 PM CEST  
 *
 */


typedef struct segm{

  double prob;
  double test;
  double methA;
  double methB;
  double **value; //cpgs
  int *pos; //index
  int n;
  int start;
  int stop;
  double meandiff;
  int length;
  double sigcp;
  char *chr;
  char init;
  struct segm *next;
  struct segm *parent;
} segment_t;

typedef struct{
  char *chr;
  int start;
  int stop;
  double p;
  double q;
  double mwu;
  int n;
  double meandiff;
  int length;
  double sigcp;
  char *methA;
  char *methB;
} segment_out;

typedef struct{
  segment_out *segment_out;
  int n;
  int i;
  int numberTests;
} list_out;


typedef struct{
  segment_t *seg;
  int n;
  char *chr;
  
  int firststop;
  char *nextchr;
  int nextstart;
  segment_t *head;
  segment_t *tail;
} segmentset_t;

typedef struct{
  int start;
  int stop;
  char *chr;
  double *groupA;
  double *groupB;
  int noA;
  int noB;
  double methA;
  double methB;
} cpg_t;

typedef struct{
  int maxdist;
  int maxseg;
  int mincpgs;
  int minDMR; // newcodes: min length of one DMR
  double minDMR2; // newcodes: min length of one DMR
  int threads;
  int mode;
  int mtc;
  // char *nameA;
  // char *nameB;
  double trend;
  double mindiff; // newcodes: min difference in one CpG site
  double mindiff2; // newcodes: min difference in one CpG site
  double minFactor; // factor for setting min number of present values
  double minMethDist;
  int minNoA; //min number of present values to fill missing numbers in group A, below discard input line
  int minNoB;
  double valley;
  int groups;// newcodes
  int **subgroupID;// newcodes
  int *subgroupSize;// newcodes
  int ***groupID;// newcodes
  int **groupSize;// newcodes
  int groupNumber;// newcodes

  int clustering;// newcodes
  int outputImputed;// newcodes

  //only used for threaded segmentation
  char **chr;
  int *pos;
  double **value;
  int n;
  int *grpA;
  int noA;
  int *grpB;
  int noB;
  double ***MWU;
  segment_t *seg; 
  cpg_t *cpg;
  segment_t *List;
  int nList;
  list_out *outputList;
  
  int threadno;
  int randomseed;

} metseg_t;
// zzhu$ metseg_t: parameters for the whole process and input data.

void initSegment(segment_t *seg);


typedef struct{
  int a;
  int b;
  int ab1;
  int ab2;
  int child;
  segment_t *max;
  
  double ks11;
  double ks12;
  double ks13;
  double ks14;
  
  double ks21;
  double ks22;
  double ks23;
  double ks24;
  
  double ks31;
  double ks32;
  double ks33;
  double ks34;

  double KS1;
  double KS2;
  double KS3;
  double KS4;

} segment_p_t;


#endif
