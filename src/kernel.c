#include "kernel.h"

Kernel *CreateKernel(AdjRel *A)
{
  Kernel *K=NULL;
  int     i;

  /* Create Kernel */

  K = (Kernel *)calloc(1,sizeof(Kernel ));
  K->A = CreateAdjRel(A->n);
  K->w = AllocFloatArray(A->n);

  /* Copy adjacency relation */

  for (i=0; i < A->n; i++) 
    K->A->adj[i] = A->adj[i];

  return(K);
} 

void DestroyKernel(Kernel *K)
{

 if (K!=NULL){
   DestroyAdjRel(K->A);
   free(K->w);
   free(K);
 }

}

Kernel *GaussianKernel(float r)
{
  Kernel *K=NULL;
  AdjRel *A=Spheric(r);
  int     i;
  float   dx,dy,dz,sigma=r/3.0,S=2*sigma*sigma;
  float   sum=0.0;

  K = CreateKernel(A);
  DestroyAdjRel(A);

  /* compute coefficients */

  for (i=0; i < K->A->n; i++) {
    dx      = K->A->adj[i].dx;
    dy      = K->A->adj[i].dy;
    dz      = K->A->adj[i].dz;
    K->w[i] = expf(-(dx*dx + dy*dy + dz*dz)/S);
    sum    += K->w[i];
  }

  /* normalize coefficients */

  for (i=0; i < K->A->n; i++) {
    K->w[i] /= sum;
  }

  return(K);
}

