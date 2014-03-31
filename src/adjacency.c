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


AdjRel *Spheric(float r)
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

  /* shift to right and place central voxel at first */

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


