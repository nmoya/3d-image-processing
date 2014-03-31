#ifndef _ADJACENCY_H_
#define _ADJACENCY_H_

/* Translation-invariant adjacency relations */

#include "common.h"

/* Displacements to an adjacent voxel */

typedef struct _adjvoxel {
  int dx,dy,dz;    
} AdjVoxel;

/* Adjacency relation */

typedef struct _adjrel { 
  AdjVoxel *adj; /* list of displacements to the adjacent voxels */
  int        n;  /* number of adjacent voxels */
} AdjRel;

Voxel   GetAdjacentVoxel(AdjRel *A, Voxel u, int adj); /* Access adjacent voxel */

AdjRel *CreateAdjRel(int n);      /* Allocate memory to store
				     adjacency relation */ 
void    DestroyAdjRel(AdjRel *A); /* Free allocated memory */

/* ------------------ Translation-Invariant Adjacency Relations --------------*/

AdjRel *Spheric(float r); /* Create a ball adjacency with radius r >=
			     1.0 voxel */

#endif
