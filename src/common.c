#include "common.h"


char *AllocCharArray(int n)
{
  char *v=NULL;
  v = (char *) calloc(n,sizeof(char));
  if (v == NULL)
    Error(MSG1,"AllocCharArray");
  return(v);
}

uchar *AllocUCharArray(int n)
{
  uchar *v=NULL;
  v = (uchar *) calloc(n,sizeof(uchar));
  if (v == NULL)
    Error(MSG1,"AllocUCharArray");
  return(v);
}

short *AllocShortArray(int n)
{
  short *v=NULL;
  v = (short *) calloc(n,sizeof(short));
  if (v == NULL)
    Error(MSG1,"AllocShortArray");
  return(v);
}

ushort *AllocUShortArray(int n)
{
  ushort *v=NULL;
  v = (ushort *) calloc(n,sizeof(ushort));
  if (v == NULL)
    Error(MSG1,"AllocUShortArray");
  return(v);
}

int *AllocIntArray(int n)
{
  int *v=NULL;
  v = (int *) calloc(n,sizeof(int));
  if (v == NULL)
    Error(MSG1,"AllocIntArray");
  return(v);
}

uint *AllocUIntArray(int n)
{
  uint *v=NULL;
  v = (uint *) calloc(n,sizeof(uint));
  if (v == NULL)
    Error(MSG1,"AllocUIntArray");
  return(v);
}

ullong *AllocULLongArray(int n)
{
  ullong *v=NULL;
  v = (ullong *) calloc(n,sizeof(ullong));
  if (v == NULL)
    Error(MSG1,"AllocULLongArray");
  return(v);
}

float *AllocFloatArray(int n)
{
  float *v=NULL;
  v = (float *) calloc(n,sizeof(float));
  if (v == NULL)
    Error(MSG1,"AllocFloatArray");
  return(v);
}

double *AllocDoubleArray(int n)
{
  double *v=NULL;
  v = (double *) calloc(n,sizeof(double));
  if (v == NULL)
    Error(MSG1,"AllocDoubleArray");
  return(v);
}

Complex *AllocComplexArray(int n)
{
  Complex *v=NULL;
  v = (Complex *) calloc(n,sizeof(Complex));
  if (v == NULL)
    Error(MSG1,"AllocComplexArray");
  return(v);
}

long double *AllocLongDoubleArray(int n)
{
  long double *v=NULL;
  v = (long double *) calloc(n,sizeof(long double));
  if (v == NULL)
    Error(MSG1,"AllocLongDoubleArray");
  return(v);
}

void Error(char *msg,char *func){ 
  fprintf(stderr,"Error in %s: \n%s. \n",func,msg);
  exit(-1);
}

void Warning(char *msg,char *func){ 
  fprintf(stdout,"Warning in %s: \n%s. \n",func,msg);
}

timer *Tic(){ /* It marks the initial time */
  timer *tic=NULL;
  tic = (timer *)malloc(sizeof(timer));
  gettimeofday(tic,NULL); 
  return(tic);
}

timer *Toc(){ /* It marks the final time */
  timer *toc=NULL;
  toc = (timer *)malloc(sizeof(timer));
  gettimeofday(toc,NULL);
  return(toc);
}

float CompTime(timer *tic, timer *toc) /* It computes the time difference */
{ 
  float t=0.0;
  if ((tic!=NULL)&&(toc!=NULL)){
    t = (toc->tv_sec-tic->tv_sec)*1000.0 + 
      (toc->tv_usec-tic->tv_usec)*0.001;
    free(tic);free(toc);
    tic=NULL; toc=NULL;
  }
  return(t);
}
