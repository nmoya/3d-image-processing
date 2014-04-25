#include "matrix.h"

Matrix        *CreateMatrix(int nrows, int ncols)
{
  Matrix *M = NULL;
  M = (Matrix*) malloc(sizeof(Matrix));
  M->n = nrows * ncols;
  M->nrows = nrows;
  M->ncols = ncols;
  M->val = AllocDoubleArray(nrows*ncols);
  return M;
}
void          DestroyMatrix(Matrix *M)
{
  if (M != NULL)
  {
    free(M->val);
    free(M);
    M = NULL;
  }
}
Matrix      *MatrixMultiply(Matrix *A, Matrix *B)
{
  /*NOT WORKING*/
  if (A->ncols != B->nrows)
    Error("Incompatible matrices for multiplication", "MatrixMultiply");

  Matrix *M = CreateMatrix(A->nrows, B->ncols);
  int i, j, k;
  float sum=0, val_A, val_B;

  for ( i = 0 ; i < M->nrows ; i++ )
  {
    for ( j = 0 ; j < M->ncols ; j++ )
    {
      M->val[GetMatrixIndex(M, i, j)] = 0.0;
      sum = 0.0;
      for ( k = 0 ; k < A->ncols ; k++ )
      {
        val_A = A->val[GetMatrixIndex(A, i, k)];
        //printf("Pegando o valor %f de A[%d, %d]\n", val_A, i, k);
        val_B = B->val[GetMatrixIndex(B, k, j)];
        //printf("Pegando o valor %f de B[%d, %d]\n", val_B, k, j);
        sum = sum + (val_A * val_B);
      }
      //printf("Colocando %f em M[%d, %d]\n", sum, i, j);
      M->val[GetMatrixIndex(M, i, j)] = sum;
    }
  }
  return M;
}
// Matrix *MultMatrices(Matrix *A, Matrix *B)
// {
//   Matrix *M = NULL; /* M = alpha A*B + beta M */
//   double  alpha=1.0, beta=0.0; 
 
//   if(A->ncols!=B->nrows)
//     Error("Cannot multiply matrices","iftMultMatrices");

//    //Compute multiplication between matrices 

//   M = CreateMatrix(A->nrows, B->ncols);

//   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->nrows, B->ncols, \
//         A->ncols, alpha, A->val, A->ncols, B->val, B->ncols, beta, \
//         M->val, B->ncols);

//   return(M);
// }
Matrix      *TranslationMatrix(float tx, float ty, float tz)
{
  Matrix *A;

  A = CreateMatrix(4,4);

  int i, j;
  for (i=0; i<=2; i++)
  {
    for(j=0; j<=2; j++)
    {
      if (i == j)
        A->val[GetMatrixIndex(A,i,j)] = 1.0;
      else
        A->val[GetMatrixIndex(A,i,j)] = 0.0;
    }
  }
  A->val[GetMatrixIndex(A,0,3)] = tx;
  A->val[GetMatrixIndex(A,1,3)] = ty;
  A->val[GetMatrixIndex(A,2,3)] = tz;
  A->val[GetMatrixIndex(A,3,3)] = 1.0;

  return(A);
}
Matrix      *RotationMatrix(char axis, float angle)
{
  Matrix *A;
  float cosine,sine, radians;

  A         = CreateMatrix(4,4);
  radians   = angle * PI/180.0;
  cosine    = cosf(radians);
  sine      = sinf(radians);

  if (axis == 'X')
  {
    A->val[GetMatrixIndex(A,0,0)] = 1.0;
    A->val[GetMatrixIndex(A,1,0)] = 0.0;
    A->val[GetMatrixIndex(A,2,0)] = 0.0;
    A->val[GetMatrixIndex(A,3,0)] = 0.0;

    A->val[GetMatrixIndex(A,0,1)] = 0.0;
    A->val[GetMatrixIndex(A,1,1)] = cosine;
    A->val[GetMatrixIndex(A,2,1)] = sine;
    A->val[GetMatrixIndex(A,3,1)] = 0.0;

    A->val[GetMatrixIndex(A,0,2)] = 0.0;
    A->val[GetMatrixIndex(A,1,2)] = -sine;
    A->val[GetMatrixIndex(A,2,2)] = cosine;
    A->val[GetMatrixIndex(A,3,2)] = 0.0;

    A->val[GetMatrixIndex(A,0,3)] = 0.0;
    A->val[GetMatrixIndex(A,1,3)] = 0.0;
    A->val[GetMatrixIndex(A,2,3)] = 0.0;
    A->val[GetMatrixIndex(A,3,3)] = 1.0;
  }
  else if (axis == 'Y')
  {
    A->val[GetMatrixIndex(A,0,0)] = cosine;
    A->val[GetMatrixIndex(A,1,0)] = 0.0;
    A->val[GetMatrixIndex(A,2,0)] = sine;
    A->val[GetMatrixIndex(A,3,0)] = 0.0;

    A->val[GetMatrixIndex(A,0,1)] = 0.0;
    A->val[GetMatrixIndex(A,1,1)] = 1.0;
    A->val[GetMatrixIndex(A,2,1)] = 0.0;
    A->val[GetMatrixIndex(A,3,1)] = 0.0;

    A->val[GetMatrixIndex(A,0,2)] = -sine;
    A->val[GetMatrixIndex(A,1,2)] = 0.0;
    A->val[GetMatrixIndex(A,2,2)] = cosine;
    A->val[GetMatrixIndex(A,3,2)] = 0.0;

    A->val[GetMatrixIndex(A,0,3)] = 0.0;
    A->val[GetMatrixIndex(A,1,3)] = 0.0;
    A->val[GetMatrixIndex(A,2,3)] = 0.0;
    A->val[GetMatrixIndex(A,3,3)] = 1.0;
  }
  else if (axis == 'Z')
  {
    A->val[GetMatrixIndex(A,0,0)] = cosine;
    A->val[GetMatrixIndex(A,1,0)] = -sine;
    A->val[GetMatrixIndex(A,2,0)] = 0.0;
    A->val[GetMatrixIndex(A,3,0)] = 0.0;

    A->val[GetMatrixIndex(A,0,1)] = sine;
    A->val[GetMatrixIndex(A,1,1)] = cosine;
    A->val[GetMatrixIndex(A,2,1)] = 0.0 ;
    A->val[GetMatrixIndex(A,3,1)] = 0.0;

    A->val[GetMatrixIndex(A,0,2)] = 0.0;
    A->val[GetMatrixIndex(A,1,2)] = 0.0;
    A->val[GetMatrixIndex(A,2,2)] = 1.0;
    A->val[GetMatrixIndex(A,3,2)] = 0.0;

    A->val[GetMatrixIndex(A,0,3)] = 0.0;
    A->val[GetMatrixIndex(A,1,3)] = 0.0;
    A->val[GetMatrixIndex(A,2,3)] = 0.0;
    A->val[GetMatrixIndex(A,3,3)] = 1.0;
  }
  else
    Error("Wrong axis. Please provide X, Y or Z.", "RotationMatrix");

  return(A);
}
float      MatrixDot(Matrix *A, Matrix *B)
{
  float result = 0;

  if (A->nrows != 1 && B->nrows != 1)
    Error("Cannot dot matrices", "Dot");
  if (A->ncols != B->ncols)
    Error("Cannot dot matrices", "Dot");

  for (int i = 0; i < A->ncols; i++)
    result += (A->val[i] * B->val[i]);

  return result;
}


Matrix      *ScaleMatrix(float sx, float sy, float sz)
{
  Matrix *A;

  A = CreateMatrix(4,4);
  float scales[3];
  scales[0] = sx;
  scales[1] = sy;
  scales[2] = sz;

  int i, j, current_axis;
  current_axis = 0;
  for(i=0; i<=2; i++)
  {
    for(j=0; j<=2; j++)
    {
      if (i == j)
        A->val[GetMatrixIndex(A,i,j)] = scales[current_axis++];
      else
        A->val[GetMatrixIndex(A,i,j)] = 0.0;
    }
  }  
  A->val[GetMatrixIndex(A,0,3)] = 0.0;
  A->val[GetMatrixIndex(A,1,3)] = 0.0;
  A->val[GetMatrixIndex(A,2,3)] = 0.0;
  A->val[GetMatrixIndex(A,3,3)] = 1.0;

  return(A);
}
void PrintMatrix(Matrix *M)
{
  int i,c,r;

  i=0; 
  fprintf(stdout,"\n");
  for (r=0; r < M->nrows; r++)  {
    for (c=0; c < M->ncols; c++) {
      fprintf(stdout,"%6.5lf ", M->val[i]);
      i++;
    }      
    fprintf(stdout,"\n");
  }
}

Matrix      *VoxelToMatrix(Voxel v)
{
  Matrix *M = CreateMatrix(4, 1);
  M->val[AXIS_X] = v.x;
  M->val[AXIS_Y] = v.y;
  M->val[AXIS_Z] = v.z;
  M->val[AXIS_H] = 1.0;
  return M;
}
void       DestroyCubeFaces(CubeFaces *cf)
{
  int i = 0;
  int number_of_faces_in_a_cube = 6;
  if (cf != NULL)
  {
    for (i=0; i< number_of_faces_in_a_cube; i++)
    {
      DestroyMatrix(cf[i].orthogonal);
      DestroyMatrix(cf[i].center);
    }
    free(cf);
  }
}