
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>

int compare3(const void *x, const void *y)
{
  double *a, *b;
  a=(double*)x;
  b=(double*)y;
  return ((*a>*b)-(*a<*b));
}



void unpaired(int *id, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, int *which1, int *which0, double *s0, double *e, double *fudge)
{
  double *ex1, *ex0, *ex21, *ex20, *r, *s, *ssort;
  int i, j;

  ex1=calloc(*ngene,sizeof(double));
  ex0=calloc(*ngene,sizeof(double));
  ex21=calloc(*ngene,sizeof(double));
  ex20=calloc(*ngene,sizeof(double));
  r=calloc(*ngene,sizeof(double));
  s=calloc(*ngene,sizeof(double));
  ssort=calloc(*ngene,sizeof(double));
  

  /* compute first and second moment to calculate mean and variance */
  for (j=0; j<*ngene; j++){
    for (i=0; i<*nsample; i++){ 
      if (id[i]==0){ex0[j] += matrix[j*(*nsample)+i];}
      if (id[i]==1){ex1[j] += matrix[j*(*nsample)+i];}
    }
    for (i=0; i<*nsample; i++){ 
      if (id[i]==0){ex20[j] += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];}
      if (id[i]==1){ex21[j] += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];}
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
      e[j]=r[j]/s[j];
    }
    
    /* fold change equivalent */
    if (*meth==3){
      e[j]=r[j];
    }    
    
    ssort[j]=s[j];
  }
  
  /* Z test statistic, needs calculation of median(s) */
  if (*meth==2){
    
    if (*s0==0){
      
      qsort((void*)ssort,*ngene,sizeof(double),compare3);
      
      if (fmod(*ngene,2)==1){
	*s0=ssort[(*ngene-1)/2];
      }
      if (fmod(*ngene,2)==0){
	*s0=(ssort[(*ngene)/2]+ssort[(*ngene)/2-1])/2;
      }
    }
    
    for (j=0; j<*ngene; j++){
      e[j]=r[j]/(s[j] + *s0);
    }
  }
  
  fudge[0]=*s0;

  free(ex1);
  free(ex0);
  free(ex21);
  free(ex20);
  free(r);
  free(s);
  free(ssort);
}




void paired(int *id, int *n1, int *n0, double *matrix, int *ngene, int *nsample, int *meth, int *which1, int *which0, double *s0, double *e, double *fudge)
{
  double *r, *s, *ssort, *ex2;
  double *diff;
  int i, j;

  diff=calloc(*n1,sizeof(double));
  r=calloc(*ngene,sizeof(double));
  s=calloc(*ngene,sizeof(double));
  ssort=calloc(*ngene,sizeof(double));
  ex2=calloc(*ngene,sizeof(double));

  for (j=0; j<*ngene; j++){
    
    /* especially important for permuted values: where are the original pairs and in which order do we have to substract them? */

    /* compute pairwise differences and second moment for variance calculation */
    for (i=0; i<*n0; i++){
      if (id[which0[i]]==0){diff[i]=matrix[j*(*nsample)+which1[i]]-matrix[j*(*nsample)+which0[i]];}
      if (id[which0[i]]==1){diff[i]=matrix[j*(*nsample)+which0[i]]-matrix[j*(*nsample)+which1[i]];}

      r[j] += diff[i];
      ex2[j] += diff[i]*diff[i];
    }
        
    r[j]=r[j]/(*n1);
    ex2[j]=ex2[j]/(*n1);

    s[j]=(*n1)*(ex2[j]-r[j]*r[j]);
    s[j]=sqrt( s[j]/((*n1)*((*n1)-1)) );
    
    /* t test statistic */
    if (*meth==1){
      e[j]=r[j]/s[j];
    }

    /* fold change equivalent */
    if (*meth==3){
      e[j]=r[j];
    }    
    
    ssort[j]=s[j];
  }

  /* Z test statistic, needs calculation of median(s) */
  if (*meth==2){

    if (*s0==0){
      qsort((void*)ssort,*ngene,sizeof(double),compare3);
      
      if (fmod(*ngene,2)==1){
	*s0=ssort[(*ngene-1)/2];
      }
      if (fmod(*ngene,2)==0){
	*s0=(ssort[(*ngene)/2]+ssort[(*ngene)/2-1])/2;
      }
    }

    for (j=0; j<*ngene; j++){
      e[j]=r[j]/(s[j] + *s0);
    }
  }

  fudge[0]=*s0;

  free(diff);
  free(r);
  free(s);
  free(ssort);
  free(ex2);
}
