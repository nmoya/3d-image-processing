#ifndef _KERNEL_H_
#define _KERNEL_H_

#include "common.h"
#include "adjacency.h"

/* Kernel definitions */

typedef struct _kernel { 
  AdjRel *A; /* adjacency relation  */
  float  *w; /* a list of weights to the adjacent voxels */
} Kernel;


Kernel   *CreateKernel(AdjRel *A); /* Allocate memory to store a kernel  */
 
void      DestroyKernel(Kernel *K); /* Free memory */

Kernel   *GaussianKernel(float r); /* Create a Gaussian kernel with radius r and ball shape */


#endif
