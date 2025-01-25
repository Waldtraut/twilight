
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>

int compare1(const void *x, const void *y)
{
  double *a, *b;
  a=(double*)x;
  b=(double*)y;
  return ((*a>*b)-(*a<*b));
}


void unpairedci(int *id, int *nperm, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, double *sobs, int *which1, int *which0, double *s0, double *e)
{
  double *ex1, *ex0, *ex21, *ex20, *r, *s, *stat;
  int i, j, k;
    
  ex1=calloc(*ngene,sizeof(double));
  ex0=calloc(*ngene,sizeof(double));
  ex21=calloc(*ngene,sizeof(double));
  ex20=calloc(*ngene,sizeof(double));
  r=calloc(*ngene,sizeof(double));
  s=calloc(*ngene,sizeof(double));
  stat=calloc(*ngene,sizeof(double));


  for (k=0; k<*nperm; k++){

    for (j=0; j<*ngene; j++){
      ex1[j]=0;
      ex0[j]=0;
      ex21[j]=0;
      ex20[j]=0;
      r[j]=0;
      s[j]=0;
      stat[j]=0;
    }

    /* compute first and second moment to calculate mean and variance */
    for (j=0; j<*ngene; j++){
      for (i=0; i<*nsample; i++){
	if (id[k*(*nsample)+i]==0){
	  ex0[j] += matrix[j*(*nsample)+i];
	  ex20[j] += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];
	}
	if (id[k*(*nsample)+i]==1){
	  ex1[j] += matrix[j*(*nsample)+i];
	  ex21[j] += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];
	}
      }
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
      
      /* Z test statistic */
      if (*meth==2){
	stat[j]=r[j]/(s[j] + *s0);
      }

      /* fold change equivalent */
      if (*meth==3){
	stat[j]=r[j];
      }
      
    }
    

    /* Maximum absolute difference between stat and sobs */
    qsort((void*)stat,*ngene,sizeof(double),compare1);
    qsort((void*)sobs,*ngene,sizeof(double),compare1);

    for (j=0; j<*ngene; j++){
      stat[j]-=sobs[j];
      stat[j]=fabs(stat[j]);
    }

    qsort((void*)stat,*ngene,sizeof(double),compare1);

    e[k]=stat[*ngene-1];
    
  }
  
  free(ex1);
  free(ex0);
  free(ex21);
  free(ex20);
  free(r);
  free(s);
  free(stat);
}




void pairedci(int *id, int *nperm, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, double *sobs, int *which1, int *which0, double *s0, double *e)
{
  double *r, *s, *ex2, *stat;
  double *g1, *g0, *diff;
  int i, j, k1, k0, k;

  g1=calloc(*n1,sizeof(double));
  g0=calloc(*n1,sizeof(double));
  diff=calloc(*n1,sizeof(double));
  r=calloc(*ngene,sizeof(double));
  s=calloc(*ngene,sizeof(double));
  ex2=calloc(*ngene,sizeof(double));
  stat=calloc(*ngene,sizeof(double));

  for (k=0; k<*nperm; k++){

    for (j=0; j<*ngene; j++){
      r[j]=0;
      s[j]=0;
      ex2[j]=0;
      stat[j]=0;
    }
    for (j=0; j<*n1; j++){
      g1[j]=0;
      g0[j]=0;
      diff[j]=0;
    }

    for (j=0; j<*ngene; j++){

      /* especially important for permuted values: where are the original pairs and in which order do we have to substract them? */

      k1=0;
      k0=0;
      
      for (i=0; i<*n0; i++){
	if (id[k*(*nsample)+which0[i]]==0){
	  g0[k0]=matrix[j*(*nsample)+which0[i]];
	  k0=k0+1;
	}
      }
      for (i=0; i<*n1; i++){
	if (id[k*(*nsample)+which1[i]]==0){
	  g0[k0]=matrix[j*(*nsample)+which1[i]];
	  k0=k0+1;
	}
      }
      
      for (i=0; i<*n1; i++){
	if (id[k*(*nsample)+which1[i]]==1){
	  g1[k1]=matrix[j*(*nsample)+which1[i]];
	  k1=k1+1;
	}
      }
      for (i=0; i<*n0; i++){
	if (id[k*(*nsample)+which0[i]]==1){
	  g1[k1]=matrix[j*(*nsample)+which0[i]];
	  k1=k1+1;
	}
      }
      
      
      /* compute pairwise differences and second moment for variance calculation */
      for (i=0; i<*n1; i++){
	diff[i] = g1[i]-g0[i];
	
	r[j] += diff[i];
	ex2[j] += diff[i]*diff[i];
      }
      
      r[j]=r[j]/(*n1);
      ex2[j]=ex2[j]/(*n1);
      
      s[j]=(*n1)*(ex2[j]-r[j]*r[j]);
      s[j]=sqrt( s[j]/((*n1)*((*n1)-1)) );
    
      /* t test statistic */
      if (*meth==1){
	stat[j]=r[j]/s[j];
      }

      /* Z test statistic */
      if (*meth==2){      
	stat[j]=r[j]/(s[j] + *s0);
      }

      /* fold change equivalent */
      if (*meth==3){
	stat[j]=r[j];
      }
      
    }
    

    /* Maximum absolute difference between stat and sobs */
    qsort((void*)stat,*ngene,sizeof(double),compare1);
    qsort((void*)sobs,*ngene,sizeof(double),compare1);

    for (j=0; j<*ngene; j++){
      stat[j]-=sobs[j];
      stat[j]=fabs(stat[j]);
    }

    qsort((void*)stat,*ngene,sizeof(double),compare1);

    e[k]=stat[*ngene-1];
    
  }

  free(g1);
  free(g0);
  free(diff);
  free(r);
  free(s);
  free(ex2);
  free(stat);
}

