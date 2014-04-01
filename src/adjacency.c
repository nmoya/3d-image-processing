#include "adjacency.h"

Voxel GetAdjacentVoxel(AdjRel *A, Voxel u, int i)
{
  Voxel v;

  v.x = u.x + A->adj[i].dx;
  v.y = u.y + A->adj[i].dy;
  v.z = u.z + A->adj[i].dz;

  return(v);
}

AdjRel *CreateAdjRel(int n)
{
   AdjRel *A=NULL;

   A      = (AdjRel *)   calloc(1,sizeof(AdjRel ));
   A->adj = (AdjVoxel *) calloc(n,sizeof(AdjVoxel));
   A->n   = n;

  return(A);
} 

void DestroyAdjRel(AdjRel *A)
{
  if (A != NULL) {
    free(A->adj);
    free(A);
  }
}

AdjRel  *Spheric(float r) 
{
  AdjRel *A=NULL;
  int i,j,k,n,r0,d,dx,dy,dz,i0=0;
  float *dr,aux,r2;

  n=0;
  r0  = (int)r;
  r2  = r*r;
  for(dz=-r0;dz<=r0;dz++)
    for(dy=-r0;dy<=r0;dy++)
      for(dx=-r0;dx<=r0;dx++)
      if(((dx*dx)+(dy*dy)+(dz*dz)) <= r2)
  n++;
  
  A = CreateAdjRel(n);
  i=0;
  for(dz=-r0;dz<=r0;dz++)
    for(dy=-r0;dy<=r0;dy++)
      for(dx=-r0;dx<=r0;dx++)
  if(((dx*dx)+(dy*dy)+(dz*dz)) <= r2){
    A->adj[i].dx=dx;
    A->adj[i].dy=dy;
    A->adj[i].dz=dz;
    if ((dx==0)&&(dy==0)&&(dz==0))
      i0 = i;
    i++;    
  }

  for (i=i0; i > 0; i--) {
    dx = A->adj[i].dx;
    dy = A->adj[i].dy;
    dz = A->adj[i].dz;
    A->adj[i].dx = A->adj[i-1].dx;
    A->adj[i].dy = A->adj[i-1].dy;
    A->adj[i].dz = A->adj[i-1].dz;
    A->adj[i-1].dx = dx;
    A->adj[i-1].dy = dy;
    A->adj[i-1].dz = dz;
  }
  
  dr = AllocFloatArray(A->n);  
  for (i=0; i < A->n; i++) {
    dr[i] = A->adj[i].dx*A->adj[i].dx + A->adj[i].dy*A->adj[i].dy + A->adj[i].dz*A->adj[i].dz;  
  }

  for (i=1; i < A->n-1; i++){
    k = i;
    for (j=i+1; j < A->n; j++)
      if (dr[j] < dr[k]){
  k = j;
      }
    aux   = dr[i];
    dr[i] = dr[k];
    dr[k] = aux;
    d            = A->adj[i].dx;
    A->adj[i].dx = A->adj[k].dx;
    A->adj[k].dx = d;
    d            = A->adj[i].dy;
    A->adj[i].dy = A->adj[k].dy;
    A->adj[k].dy = d;
    d            = A->adj[i].dz;
    A->adj[i].dz = A->adj[k].dz;
    A->adj[k].dz = d;
  }

  free(dr);
  return(A);
}

AdjRel *Spheric2(float r)
{
  AdjRel *A=NULL;
  int     dx=0,dy=0,dz,n,i,R=(int)r,i0=0;
  float   R2=r*r;

  if (r < 1.0)
    Error("Radius must be >= 1.0","Spheric");
  
  /* Compute the number of adjacent voxels */

  n=0;
  for (dz=-R; dz <= R; dz++)
    for (dy=-R; dy <= R; dy++)
      for (dx=-R; dx <= R; dx++)
	if ((dx*dx + dy*dy + dz*dz)<=R2)
	  n++;


  /* Allocate memory to store adjacent voxels */

  A = CreateAdjRel(n);

  /* Store adjacent voxels */

  i=0;
  for (dz=-R; dz <= R; dz++)
    for (dy=-R; dy <= R; dy++)
      for (dx=-R; dx <= R; dx++)
	if ((dx==0)&&(dy==0)&&(dz==0)) i0 = i;
	if ((dx*dx + dy*dy + dz*dz)<=R2){
	  A->adj[i].dx = dx;
	  A->adj[i].dy = dy;
	  A->adj[i].dz = dz;
	  i++;
	}

  /* sh to right and place central voxel at first */

  for (i=i0; i > 0; i--) {
    dx = A->adj[i].dx;
    dy = A->adj[i].dy;
    dz = A->adj[i].dz;
    A->adj[i].dx = A->adj[i-1].dx;
    A->adj[i].dy = A->adj[i-1].dy;
    A->adj[i].dz = A->adj[i-1].dz;
    A->adj[i-1].dx = dx;
    A->adj[i-1].dy = dy;
    A->adj[i-1].dz = dz;
  }

  return(A);
}


