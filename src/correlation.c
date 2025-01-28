                      
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <stdlib.h>

int compare7(const void *x, const void *y)
{
  double *a, *b;
  a=(double*)x;
  b=(double*)y;
  return ((*a>*b)-(*a<*b));
}


void corsingle(double *vector, double *matrix, int *ngene, int *nsample, double *e)
{
  double *ex0, *ex1, *ex20, *ex21, *exboth;
  int i, j;

  ex0=calloc(1,sizeof(double));
  ex1=calloc(*ngene,sizeof(double));
  ex20=calloc(1,sizeof(double));
  ex21=calloc(*ngene,sizeof(double));
  exboth=calloc(*ngene,sizeof(double));

  /* compute first and second moments */
  for (i=0; i<*nsample; i++){
    ex0[0]  += vector[i];
    ex20[0] += vector[i]*vector[i];
  }

  for (j=0; j<*ngene; j++){
    for (i=0; i<*nsample; i++){
      ex1[j]    += matrix[j*(*nsample)+i];
      ex21[j]   += matrix[j*(*nsample)+i]*matrix[j*(*nsample)+i];
      exboth[j] += matrix[j*(*nsample)+i]*vector[i];
    }
    /* correlation coefficient */
    e[j]=(exboth[j]-ex0[0]*ex1[j]/(*nsample))/sqrt((ex20[0]-ex0[0]*ex0[0]/(*nsample))*(ex21[j]-ex1[j]*ex1[j]/(*nsample)));
  }
  
  free(ex0);
  free(ex1);
  free(ex20);
  free(ex21);
  free(exboth);
}



void corperm(double *vecperm, int *nperm, double *matrix, int *ngene, int *nsample, double *e, double *f)
{
  double *ex0, *ex1, *ex20, *ex21, *exboth, *stat, *dstat;
  int i, j, k, *indx;

  ex0=calloc(1,sizeof(double));
  ex1=calloc(*ngene,sizeof(double));
  ex20=calloc(1,sizeof(double));
  ex21=calloc(*ngene,sizeof(double));
  exboth=calloc(*ngene,sizeof(double));
  stat=calloc(*ngene,sizeof(double));
  indx=calloc((*nperm)*(*ngene),sizeof(int));
  dstat=calloc((*nperm)*(*ngene),sizeof(double));


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
    
    for (j=0; j<*ngene; j++){
      dstat[j+(*ngene)*k]=fabs(stat[j]);
    }

    qsort((void*)stat,*ngene,sizeof(double),compare7);
    for (j=0; j<*ngene; j++){
      e[j]+=stat[j];
    }
  }


  /* Compute p-values */    
  for (j=0; j<(*nperm)*(*ngene); j++){
    indx[j]=j;
  }

  rsort_with_index((double*)dstat,(int*)indx,(int)(*nperm)*(*ngene));
  
  /* First ngene values correspond to the original labeling */
  for (j=0; j<(*nperm)*(*ngene); j++){      		
    if (indx[j]<*ngene){
      f[indx[j]] = (*nperm)*(*ngene)-j;
    }
  }   

  for (j=0; j<*ngene; j++){ 
    e[j]=e[j]/(*nperm);
    f[j]=f[j]/((*nperm)*(*ngene));
  }
    
  free(ex0);
  free(ex1);
  free(ex20);
  free(ex21);
  free(exboth);
  free(stat);
  free(dstat);
  free(indx);
}



void corci(double *vecperm, int *nperm, double *matrix, int *ngene, int *nsample, double *sobs, double *e)
{
  double *ex0, *ex1, *ex20, *ex21, *exboth, *stat;
  int i, j, k;

  ex0=calloc(1,sizeof(double));
  ex1=calloc(*ngene,sizeof(double));
  ex20=calloc(1,sizeof(double));
  ex21=calloc(*ngene,sizeof(double));
  exboth=calloc(*ngene,sizeof(double));
  stat=calloc(*ngene,sizeof(double));

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
        
    /* Maximum absolute difference between stat and sobs */
    qsort((void*)stat,*ngene,sizeof(double),compare7);
    qsort((void*)sobs,*ngene,sizeof(double),compare7);

    for (j=0; j<*ngene; j++){
      stat[j]-=sobs[j];
      stat[j]=fabs(stat[j]);
    }

    qsort((void*)stat,*ngene,sizeof(double),compare7);

    e[k]=stat[*ngene-1];     
  }

    
  free(ex0);
  free(ex1);
  free(ex20);
  free(ex21);
  free(exboth);
  free(stat);
}


