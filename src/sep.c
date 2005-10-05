
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>


/* Returns maximum absolute difference of empirical and uniform distribution */
double empirical(double *a, int na)
{

  int i;
  double *b;
  double e;

  if ((b=malloc(sizeof(double)*1))==0) {printf("Error, could not allocate memory");}

  e=0;

  for (i=1; i<na; i++){

    if (a[i]!=a[i-1]){
      b[0]=fabs( a[i] - ( (double) i/ ((double) na) ) );
      if (b[0]>e){e=b[0];}
      b[0]=fabs( a[i] - ( ((double) i+1)/ ((double) na) ) );
      if (b[0]>e){e=b[0];}
    }

  }

  free(b);

  return e;
}


/* SEP */
void sep(double *xin, int *nxin, double *lambda, int *xout, double *funout)
{

  int *ix;
  int i;
  int j=0;
  int count=0;
  int nyin=0;
  int randnum;
  double *yin;
  double *objfunc;

  if ((ix=malloc(sizeof(int)*(*nxin)))==0) {printf("Error, could not allocate memory");}
  if ((objfunc=malloc(sizeof(double)*2))==0) {printf("Error, could not allocate memory");}

  for (i=0; i<*nxin; i++){ 
    ix[i]=1; 
  }

  objfunc[0]=empirical(xin,*nxin);
  
  if (objfunc[0]<=0.1)
    {
      funout[0]=objfunc[0];
    }

  if (objfunc[0]>0.1)
    {
      while (objfunc[0]>0.1 && count<2*(*nxin))
	{
	  randnum=(int)( ((double) *nxin)*rand()/(RAND_MAX+1.0) );
	  
	  ix[randnum]=fabs(1-ix[randnum]);
	  
	  for (i=0; i<*nxin; i++){ 
	    if (ix[i]==1){
	      nyin++;
	    } 
	  }
	  
	  if ((yin=malloc(sizeof(double)*nyin))==0) {printf("Error, could not allocate memory");}
	  
	  for (i=0; i<*nxin; i++){ 
	    if (ix[i]==1){
	      yin[j]=xin[i]; 
	      j++;
	    } 
	  }
	  
	  objfunc[1]=empirical(yin,nyin);
	  
	  if (objfunc[1]<objfunc[0])
	    {
	      objfunc[0]=objfunc[1];
	      count=0;
	      funout[0]=objfunc[0] + (*lambda)*(*nxin - (double)nyin)*log(*nxin - (double)nyin)/(*nxin);
	    }
	  else
	    {
	      ix[randnum]=fabs(1-ix[randnum]);
	      count++;
	    }
	  
	  nyin=0;
	  j=0;
	  free(yin);
	}
    }



  objfunc[0] = funout[0];
  count=0;

  while (count<2*(*nxin))
    {
      randnum=(int)( ((double) *nxin)*rand()/(RAND_MAX+1.0) );
      
      ix[randnum]=fabs(1-ix[randnum]);
      
      for (i=0; i<*nxin; i++){ 
	if (ix[i]==1){
	  nyin++;
	} 
      }
      
      if ((yin=malloc(sizeof(double)*nyin))==0) {printf("Error, could not allocate memory");}
      
      for (i=0; i<*nxin; i++){ 
	if (ix[i]==1){
	  yin[j]=xin[i]; 
	  j++;
	} 
      }
      
      objfunc[1]=empirical(yin,nyin) + (*lambda)*(*nxin - (double)nyin)*log(*nxin - (double)nyin)/(*nxin);
      
      if (objfunc[1]<objfunc[0])
	{
	  objfunc[0]=objfunc[1];
	  count=0;
	  funout[0]=empirical(yin,nyin);
	}
      else
	{
	  ix[randnum]=fabs(1-ix[randnum]);
	  count++;
	}
      
      nyin=0;
      j=0;
      free(yin);
    } 
  
  
  for (i=0; i<*nxin; i++)
    {
      xout[i]=ix[i];
    }
  
  
  free(ix);
  free(objfunc);
 
}
