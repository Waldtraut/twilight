
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


int compare11(const void *x, const void *y)
{
  double *a, *b;
  a=(double*)x;
  b=(double*)y;
  return ((*a>*b)-(*a<*b));
}


/* returns indeces of *a in "first" mode */
int* indexx(int n, double *a)
{
  int *indx;
  int i;

  indx=calloc(n,sizeof(int));

  for (i=0; i<n; i++){
    a[i]=fabs(a[i]);
    indx[i]=i;
  }
  
  rsort_with_index((double*)a,(int*)indx,(int)n);

  return indx;
  /* indx will be freed in kolmogoroff with free(indx) */
}


/* Returns maximum absolute difference of empirical and uniform distribution */
/* Kolmogoroff-Smirnoff-Test of uniformity */
/* Caution: Input *a are the absolute observed values, not p-values! */
double* kolmogoroff(double *a, int nrowa, int ncola)
{

  int i;
  int k;
  double *b;
  double *asort;
  double *pval;
  int n=nrowa*ncola;
  int *indx;
  double *eout;


  b=calloc(1,sizeof(double));
  asort=calloc(nrowa,sizeof(double));
  pval=calloc(n,sizeof(double));
  eout=calloc(ncola,sizeof(double));

  indx=indexx(n,a);

  for (i=0; i<n; i++){
    pval[indx[i]]=(double)(n-i)/((double)n);   /* returns pooled p-value */
  }

  for (k=0; k<ncola; k++){
    
    for (i=0; i<nrowa; i++){
      asort[i]=pval[i*ncola+k];
    }
    
    qsort((void*)asort,nrowa,sizeof(double),compare11);    

    eout[k]=0;
    
    for (i=1; i<nrowa; i++){
      
      if (asort[i]!=asort[i-1]){
	b[0]=fabs( asort[i] - ( (double) i/ ((double) nrowa) ) );
	if (b[0]>eout[k]){eout[k]=b[0];}
	b[0]=fabs( asort[i] - ( ((double) i+1)/ ((double) nrowa) ) );
	if (b[0]>eout[k]){eout[k]=b[0];}
      }
    }    
  }

  free(b); 
  free(asort);  
  free(pval);  
  free(indx);  

  return eout;

  /* eout will be freed in unpairedKSTEST with free(e) */
}



/* Computes test statistics, pooled p-values and KS statistics */
void unpairedKSTEST(int *id, int *nperm, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, double *s0, double *f)
{
  double *ex1, *ex0, *ex21, *ex20, *r, *s, *ssort, *stat, *statall, *e;
  int i, j, k, *test;
  int nrowa=*ngene, ncola=*nperm;

  ex1=calloc(*ngene,sizeof(double));
  ex0=calloc(*ngene,sizeof(double));
  ex21=calloc(*ngene,sizeof(double));
  ex20=calloc(*ngene,sizeof(double));
  r=calloc(*ngene,sizeof(double));
  s=calloc(*ngene,sizeof(double));
  ssort=calloc(*ngene,sizeof(double));
  stat=calloc(*ngene,sizeof(double));
  statall=calloc((*ngene)*(*nperm),sizeof(double));
  test=calloc(1,sizeof(int));


  for (k=0; k<*nperm; k++){

    for (j=0; j<*ngene; j++){
      ex1[j]=0;
      ex0[j]=0;
      ex21[j]=0;
      ex20[j]=0;
      r[j]=0;
      s[j]=0;
      ssort[j]=0;
      stat[j]=0;
    }

    /* compute first and second moment to calculate mean and variance */
    for (j=0; j<*ngene; j++){

      for (i=0; i<*nsample; i++){
	test[0] = id[k*(*nsample)+i];

	if (test[0]==0){ex0[j] += matrix[j*(*nsample)+i];}
	if (test[0]==1){ex1[j] += matrix[j*(*nsample)+i];}
      }
      for (i=0; i<*nsample; i++){
	test[0] = id[k*(*nsample)+i];

	if (test[0]==0){ex20[j] += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];}
	if (test[0]==1){ex21[j] += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];}
      }

    }

    for (j=0; j<*ngene; j++){

      ex0[j]=ex0[j]/(*n0);
      ex1[j]=ex1[j]/(*n1);
      ex20[j]=ex20[j]/(*n0);
      ex21[j]=ex21[j]/(*n1);
      
      /* difference in means */
      r[j]=ex1[j] - ex0[j];

      /* pooled variance */
      s[j]=(*n1)*(ex21[j]-ex1[j]*ex1[j]) + (*n0)*(ex20[j]-ex0[j]*ex0[j]);
      s[j]=sqrt( s[j]*((double)1/(*n1)+(double)1/(*n0))/((*nsample)-2) );


      /* t test statistic */
      if (*meth==1){
	stat[j]=r[j]/s[j];
      }

      /* fold change equivalent */
      if (*meth==3){
	stat[j]=r[j];
      }
      
      ssort[j]=s[j];
    }


    /* Z test statistic, needs calculation of median(s) */
    if (*meth==2){
            
      if (*s0==0){
	qsort((void*)ssort,*ngene,sizeof(double),compare11);

	if (fmod(*ngene,2)==1){
	  *s0=ssort[(*ngene-1)/2];
	}
	if (fmod(*ngene,2)==0){
	  *s0=(ssort[(*ngene)/2]+ssort[(*ngene)/2-1])/2;
	}
      }

      for (j=0; j<*ngene; j++){
	stat[j]=r[j]/(s[j] + *s0);
      }
    }


    /* puts the ngene*nperm matrix of scores into a vector to be passed over to kolmogoroff */
    for (j=0; j<*ngene; j++){
      statall[j*(*nperm)+k]=fabs(stat[j]);
    }
  }

  e=kolmogoroff(statall,nrowa,ncola);

  for (k=0; k<*nperm; k++){
    f[k]=e[k];
  }

  free(ex1);
  free(ex0);
  free(ex21);
  free(ex20);
  free(r);
  free(s);
  free(ssort);
  free(stat);
  free(statall);
  free(test);
  free(e);
}



/* Computes test statistics, pooled p-values and KS statistics */
void pairedKSTEST(int *id, int *nperm, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, int *which1, int *which0, double *s0, double *f)
{
  double *r, *s, *ssort, *ex2, *stat, *statall, *e;
  double *diff;
  int i, j, k;
  int nrowa=*ngene, ncola=*nperm;

  diff=calloc(*n1,sizeof(double));
  r=calloc(*ngene,sizeof(double));
  s=calloc(*ngene,sizeof(double));
  ssort=calloc(*ngene,sizeof(double));
  ex2=calloc(*ngene,sizeof(double));
  stat=calloc(*ngene,sizeof(double));
  statall=calloc((*ngene)*(*nperm),sizeof(double));


  for (k=0; k<*nperm; k++){

    for (j=0; j<*ngene; j++){
      r[j]=0;
      s[j]=0;
      ssort[j]=0;
      ex2[j]=0;
      stat[j]=0;
    }
    for (j=0; j<*n1; j++){
      diff[j]=0;
    }

    for (j=0; j<*ngene; j++){

      /* compute pairwise differences and second moment for variance calculation */
      for (i=0; i<*n0; i++){

	diff[i]=matrix[j*(*nsample)+which1[i]]-matrix[j*(*nsample)+which0[i]];

	/* especially important for permuted values: where are the original pairs and in which order do we have to substract them? */
	if (id[k*(*nsample)+which0[i]]==1){diff[i]= -diff[i];}

	r[j] += diff[i];
	ex2[j] += diff[i]*diff[i];
      }
    }

    for (j=0; j<*ngene; j++){

      r[j]=r[j]/(*n1);
      ex2[j]=ex2[j]/(*n1);

      s[j]=(*n1)*(ex2[j]-r[j]*r[j]);
      s[j]=sqrt( s[j]/((*n1)*((*n1)-1)) );

      /* t test statistic */
      if (*meth==1){
	stat[j]=r[j]/s[j];
      }

      /* fold change equivalent */
      if (*meth==3){
	stat[j]=r[j];
      }
      
      ssort[j]=s[j];
    }


    /* Z test statistic, needs calculation of median(s) */
    if (*meth==2){

      if (*s0==0){
	qsort((void*)ssort,*ngene,sizeof(double),compare11);
	
	if (fmod(*ngene,2)==1){
	  *s0=ssort[(*ngene-1)/2];
	}
	if (fmod(*ngene,2)==0){
	  *s0=(ssort[(*ngene)/2]+ssort[(*ngene)/2-1])/2;
	}
      }

      for (j=0; j<*ngene; j++){
	stat[j]=r[j]/(s[j] + *s0);
      }
    }

    /* puts the ngene*nperm matrix of scores into a vector to be passed over to kolmogoroff */
    for (j=0; j<*ngene; j++){
      statall[j*(*nperm)+k]=fabs(stat[j]);
    }
  }

  e=kolmogoroff(statall,nrowa,ncola);

  for (k=0; k<*nperm; k++){
    f[k]=e[k];
  }
  

  free(diff);
  free(r);
  free(s);
  free(ssort);
  free(ex2);
  free(stat);
  free(statall);
  free(e);
}






void correlationKSTEST(double *vecperm, int *nperm, double *matrix, int *ngene, int *nsample, double *f)
{
  double *ex0, *ex1, *ex20, *ex21, *exboth, *stat, *statall, *e;
  int i, j, k;
  int nrowa=*ngene, ncola=*nperm;

  ex0=calloc(1,sizeof(double));
  ex1=calloc(*ngene,sizeof(double));
  ex20=calloc(1,sizeof(double));
  ex21=calloc(*ngene,sizeof(double));
  exboth=calloc(*ngene,sizeof(double));
  stat=calloc(*ngene,sizeof(double));
  statall=calloc((*ngene)*(*nperm),sizeof(double));


  for (k=0; k<*nperm; k++){

    ex0[0] = 0;
    ex20[0] = 0;

    for (j=0; j<*ngene; j++){
      ex1[j]=0;
      ex21[j]=0;
      exboth[j]=0;
      stat[j]=0;
    }

    /* compute first and second moments */
    for (i=0; i<*nsample; i++){
      ex0[0]  += vecperm[k*(*nsample)+i];
      ex20[0] += vecperm[k*(*nsample)+i]*vecperm[k*(*nsample)+i];
    }
   
    /* compute mixed term */
    for (j=0; j<*ngene; j++){
      for (i=0; i<*nsample; i++){
	ex1[j]    += matrix[j*(*nsample)+i];
	ex21[j]   += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];
	exboth[j] += matrix[j*(*nsample)+i]*vecperm[k*(*nsample)+i];
      }

      /* correlation coefficient */
      stat[j]=(exboth[j]-ex0[0]*ex1[j]/(*nsample))/sqrt((ex20[0]-ex0[0]*ex0[0]/(*nsample))*(ex21[j]-ex1[j]*ex1[j]/(*nsample)));
    }
    
    /* puts the ngene*nperm matrix of scores into a vector to be passed over to kolmogoroff */
    for (j=0; j<*ngene; j++){
      statall[j*(*nperm)+k]=fabs(stat[j]);
    }
  }

  e=kolmogoroff(statall,nrowa,ncola);

  for (k=0; k<*nperm; k++){
    f[k]=e[k];
  }

    
  free(ex0);
  free(ex1);
  free(ex20);
  free(ex21);
  free(exboth);
  free(stat);
  free(statall);
  free(e);
}
