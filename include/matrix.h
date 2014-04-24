#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "common.h"


typedef struct _matrix {
  float *val;
  int   ncols, nrows;
  int   n;
} Matrix;

typedef struct _cubefaces {
	Matrix *orthogonal;
	Matrix *center;
} CubeFaces;

/* Linear to 2d convertion */
#define GetMatrixCol(m,i) ((i) % (m)->ncols)
#define GetMatrixRow(m,i) ((i) / (m)->ncols)
#define GetMatrixIndex(m,c,r) ((m)->ncols*r + c)

/* Allocate memory to store image */ 
Matrix        *CreateMatrix(int nrows, int ncols);
void           DestroyMatrix(Matrix *m); /* Free memory */
Matrix 		  *VoxelToMatrix(Voxel v);
Matrix 		  *MatrixMultiply(Matrix *A, Matrix *B);
Matrix 		  *TranslationMatrix(float tx, float ty, float tz);
Matrix 		  *RotationMatrix(char axis, float angle);
Matrix 		  *ScaleMatrix(float sx, float sy, float sz);
float		   MatrixDot(Matrix *A, Matrix *B);
void		   DestroyCubeFaces(CubeFaces *cf);




#endif
