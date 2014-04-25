#include "image.h"


char ValidVoxel(Image *img, Voxel u)
{
  if ((u.x >= 0)&&(u.x < img->xsize)&&
      (u.y >= 0)&&(u.y < img->ysize)&&
      (u.z >= 0)&&(u.z < img->zsize) )
    return(1);
  else
    return(0);
}

Voxel GetVoxelCoord(Image *img, int p)
{
  Voxel u;

  u.x = GetXCoord(img,p);
  u.y = GetYCoord(img,p);
  u.z = GetZCoord(img,p);

  return(u);
}
int VoxelValue(Image *img, Voxel v)
{
  return img->val[GetVoxelIndex(img, v)];
}

Image *CreateImage(int xsize, int ysize, int zsize)
{
  Image *img=NULL;
  int    i;

  img = (Image *) calloc(1,sizeof(Image));
  img->n     = xsize*ysize*zsize;
  img->xsize = xsize; 
  img->ysize = ysize; 
  img->zsize = zsize; 
  img->dx    = img->dy = img->dz = 1.0;
  img->val   = AllocIntArray(img->n);
  img->tby   = AllocIntArray(img->ysize);
  img->tbz   = AllocIntArray(img->zsize);

  img->tby[0] = 0;
  for (i=1; i < img->ysize; i++)
    img->tby[i]=img->tby[i-1]+img->xsize;

  img->tbz[0] = 0;
  for (i=1; i < img->zsize; i++)
    img->tbz[i]=img->tbz[i-1]+img->xsize*img->ysize;

 
 return(img);
}

void DestroyImage(Image *img)
{

  if(img != NULL){
    free(img->val); 
    free(img->tby);
    free(img->tbz);
    free(img);
  }

}

Image *ReadImage(char *filename)
{
  FILE           *fp=NULL;
  unsigned short *data16=NULL;
  unsigned char  *data8=NULL;
  char            type[10];
  int             p,xsize,ysize,zsize,depth;
  Image          *img=NULL;

  fp = fopen(filename,"rb");
  if (fp == NULL){
    Error(MSG2,"ReadImage");
  }

  fscanf(fp,"%s\n",type);
  if((strcmp(type,"SCN")==0)){
    fscanf(fp,"%d %d %d\n",&xsize,&ysize,&zsize);
    img = CreateImage(xsize,ysize,zsize);
    fscanf(fp,"%f %f %f\n",&img->dx,&img->dy,&img->dz);
    fscanf(fp,"%d\n",&depth);    


    switch(depth) {
    case 8:
      data8 = AllocUCharArray(img->n); 
      fread(data8, sizeof(unsigned char), img->n, fp);
      for (p=0; p < img->n; p++) 
	img->val[p] = (int) data8[p];
      free(data8);
      break;
    case 16:
      data16 = AllocUShortArray(img->n); 
      fread(data16, sizeof(unsigned short), img->n, fp);
      for (p=0; p < img->n; p++) 
	img->val[p] = (int) data16[p];
      free(data16);
      break;
    default:
      Error("Invalid image depth","ReadImage");
    }
  }else{
    Error("Invalid image format","ReadImage");
  }


  fclose(fp);

  return(img);
}

void WriteImage(Image *img,char *filename)
{
  FILE           *fp;
  int             p, Imax, Imin;
  unsigned char  *data8; 
  unsigned short *data16; 

  fp = fopen(filename,"wb");
  if (fp == NULL){
    Error(MSG2,"WriteImage");
  }

  /* shift intensities if necessary */
  
  Imin = MinimumValue(img); 

  if (Imin < 0) {
    char msg[200];
    sprintf(msg,"Min. intensity is %d. Shifting intensities...\n",Imin);
    Warning(msg,"WriteImage");
    for (p=0; p < img->n; p++) 
      img->val[p] = img->val[p] - Imin; 
  }
    
  /* computing depth */ 

  Imax = MaximumValue(img); 
  
  if (Imax < 65536) { /* 16-bit data */

    fprintf(fp,"SCN\n");
    fprintf(fp,"%d %d %d\n",img->xsize,img->ysize,img->zsize);
    fprintf(fp,"%f %f %f\n",img->dx,img->dy,img->dz);
    fprintf(fp,"16\n");
    
    data16 = AllocUShortArray(img->n); 
    for (p=0; p < img->n; p++) 
      data16[p] = (unsigned short) img->val[p]; 

    fwrite(data16, sizeof(unsigned short), img->n, fp); 
    free(data16); 

  }else{

    if (Imax < 256) { /* 8-bit data */
      fprintf(fp,"SCN\n");
      fprintf(fp,"%d %d %d\n",img->xsize,img->ysize,img->zsize);
      fprintf(fp,"%f %f %f\n",img->dx,img->dy,img->dz);
      fprintf(fp,"8\n");
      
      data8 = AllocUCharArray(img->n); 
      for (p=0; p < img->n; p++) 
	data8[p] = (unsigned char) img->val[p]; 
      
      fwrite(data8, sizeof(unsigned char), img->n, fp); 
      free(data8); 
    } else {
      Error("Image depth is different from 8 and 16","WriteImage");
    }
  }
  fclose(fp);
}

void WriteImageP2(Image *img, char filename[])
{
    FILE *fp=NULL;
    int   p;

    fp = fopen(filename,"w"); 
    if (!fp) 
        Error("Failed to open P2 image","WriteImageP2");

    fprintf(fp,"P2\n");
    fprintf(fp,"%d %d\n",img->xsize,img->ysize);

    img->maxval = MaximumValue(img);

    fprintf(fp,"%d\n",img->maxval);
    for (p=0; p < img->n; p++) 
        fprintf(fp,"%d ",img->val[p]);
    fclose(fp);
}

/* ----------------- Voxel-based algorithms ---------------------*/

int MinimumValue(Image *img)
{
  int p, Imin = INFINITY_INT;

  for (p=0; p < img->n; p++) 
    if (img->val[p] < Imin) 
      Imin = img->val[p];

  img->minval = Imin;
  return Imin;
}

int MaximumValue(Image *img)
{
  int p, Imax = -INFINITY_INT;

  for (p=0; p < img->n; p++) 
    if (img->val[p] > Imax) 
      Imax = img->val[p];
  
  img->maxval = Imax; 
  return Imax;
}

Image  *WindowAndLevel(Image *img, int window, int level, int H)
{
  int    p, Imin, Imax; 
  Image *img_corr = CreateImage(img->xsize,img->ysize,img->zsize);  
  float  K; 

  if (window < 0) 
    Error("window must be >= 0","WindowAndLevel"); 

  if (H <= 0) 
    Error("H must be > 0","WindowAndLevel"); 

  Imin = level - window/2;
  Imax = level + window/2;

  if (window > 0)  { 
    K    = (float)H /(float)(Imax - Imin); 
  
    for (p=0; p < img->n; p++) {
      if (img->val[p] <= Imin) {
	img_corr->val[p] = 0; 
      } else {
	if (img->val[p] >= Imax) {
	  img_corr->val[p] = H; 
	} else { /* compute linear stretching */
	  img_corr->val[p] = K * (img->val[p] - Imin);
	}
      }  
    }

  } else {
    for (p=0; p < img->n; p++) {
      if (img->val[p] >= level) 
	  img_corr->val[p] = H; 
    }
  }

  return(img_corr);
}

/* ----------------- Adjacency-based algorithms ---------------------*/


Image *LinearFilter(Image *img, Kernel *K)
{
  Image *img_filt=CreateImage(img->xsize,img->ysize,img->zsize);
  int    p, q, i;
  Voxel  u, v;
  float  filt;

  for (p=0; p < img->n; p++) {
    u = GetVoxelCoord(img,p);
    filt = 0.0;
    for (i=0; i < K->A->n; i++) {
      v = GetAdjacentVoxel(K->A,u,i);
      if (ValidVoxel(img,v)){
	q = GetVoxelIndex(img,v);
	filt = filt + img->val[q] * K->w[i];
      }
    }
    img_filt->val[p] = ROUND(filt);
  }

  return(img_filt); 
}

/* ----------------- Connectivity-based algorithms ---------------------*/

Image        *LabelBinaryImage(Image *img, AdjRel *A)
{
  //QUEUE
  Voxel *QUEUE = (Voxel*) calloc(img->n, sizeof(Voxel));
  if (!QUEUE)
    Error("Error allocing queue", "LabelBinaryImage");
  int NEXT = 0; int LAST = 0;

  //Algorithm
  Image *labeled_image = CreateImage(img->xsize, img->ysize, img->zsize);
  int curr_label = 0;
  int p, q, i;
  Voxel u, v;

  for(p=0;  p< img->n; p++)
  {
    if (img->val[p] && !labeled_image->val[p]) //If it is UP and not labeled
    {
      u = GetVoxelCoord(img, p);
      curr_label++;
      
      //Insert in the queue
      QUEUE[LAST] = u;
      LAST++;

      while (NEXT != LAST) // While the queue is not empty.
      {
        u = QUEUE[NEXT];  //Remove from the queue
        NEXT++;

        q = GetVoxelIndex(img, u);
        labeled_image->val[q] = curr_label; //Set the label

        for(i=1; i<A->n; i++)
        {
          v = GetAdjacentVoxel(A,u,i);
          if (ValidVoxel(img, v))
          {
            q = GetVoxelIndex(img, v);  
            if (!labeled_image->val[q]) //Insert the uarked neighbours in the queue
            {
              QUEUE[LAST] = v;
              LAST++;
            }
          }
        }
      }
    }
  }
  free(QUEUE);
  printf("Labels: %d\n", MaximumValue(labeled_image));
  return (labeled_image);
}


/* ----------------- Task 1 ---------------------*/
Image * GetAxialSlice(Image *img, int slice)
{
    int y, x, index;
    Voxel v;
    Image *tempImg = CreateImage(img->xsize, img->ysize, 1);
    v.z = slice;
    int counter = 0;
    for (y=0; y<img->ysize; y++)
    {
        for (x=0; x<img->xsize; x++)
        {
            v.x = x; v.y = y;
            index = GetVoxelIndex(img, v);
            tempImg->val[counter] = img->val[index];
            counter++;
        }
    }
    return tempImg;
}
Image * GetCoronalSlice(Image *img, int slice)
{
    int z, x, index;
    Voxel v;
    Image *tempImg = CreateImage(img->zsize, img->xsize, 1);
    v.y = slice;
    int counter = 0;
    for (x=0; x<img->xsize; x++)
    {
        for (z=0; z<img->zsize; z++)
        {
            v.x = x; v.z = z;
            index = GetVoxelIndex(img, v);
            tempImg->val[counter] = img->val[index];
            counter++;
        }
    }
    return tempImg;
}
Image * GetSagitalSlice(Image *img, int slice)
{
    int y, z, index;
    Voxel v;
    Image *tempImg = CreateImage(img->ysize, img->zsize, 1);
    v.x = slice;
    int counter = 0;
    for (z=0; z<img->zsize; z++)
    {
        for (y=0; y<img->ysize; y++)
        {
            v.y = y; v.z = z;
            index = GetVoxelIndex(img, v);
            tempImg->val[counter] = img->val[index];
            counter++;
        }
    }
    return tempImg;
}

/* ----------------- Task 2 ---------------------*/
int CompareVoxels(Voxel v1, Voxel v2)
{
  return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}
FVoxelList * DDAAlgorithm(Voxel p1, Voxel pn)
{
  int n, k;
  float dx=0, dy=0, dz=0, DX, DY, DZ;
  FVoxelList *output = NULL;

  if (CompareVoxels(p1, pn))
    n = 1;
  else
  {
    //This computes the projection along each axis
    DX = pn.x - p1.x; DY = pn.y - p1.y; DZ = pn.z - p1.z;
    //Check for the longest projection
    if (DX >= DY && DX >= DZ)
    {
      n = DX+1;
      dx = SIGN(DX);
      dy = (dx * DY) / DX;
      dz = (dx * DZ) / DX;
    }
    else if (DY >= DX && DY >= DZ)
    {
      n = DY+1;
      dy = SIGN(DY);
      dx = (dy * DX) / DY;
      dz = (dy * DZ) / DY;
    }
    else
    {
      n = DZ + 1;
      dz = SIGN(DZ);
      dx = (dz * DX) / DZ;
      dy = (dz * DY) / DZ;
    }
  }
  //The the first position receives the first point.
  
  output = CreateFVoxelList(n);
  output->val[0].x = p1.x;
  output->val[0].y = p1.y;
  output->val[0].z = p1.z;
  
  for (k=1; k<n; k++)
  {
    output->val[k].x = output->val[k-1].x + dx;
    output->val[k].y = output->val[k-1].y + dy;
    output->val[k].z = output->val[k-1].z + dz;
  }
  return output;
}
int LinearInterpolationValue(Image *img, FVoxel v)
{
  Voxel u[8];
  int p[8], i, value;
  float dx=1.0;
  float dy=1.0;
  float dz=1.0;
  float value_array[6];
  
  if ((int) (v.x + 1.0) == img->xsize) 
    dx = 0.0; 
  if ((int) (v.y + 1.0) == img->ysize) 
    dy = 0.0;
  if ((int) (v.z + 1.0) == img->zsize) 
    dz = 0.0;
  
  //Computes the closest neighbour in each direction
  u[0].x = (int)v.x;      u[0].y = (int)v.y;       u[0].z = (int)v.z;
  u[1].x = (int)(v.x+dx); u[1].y = (int)v.y;       u[1].z = (int)v.z;
  u[2].x = (int)v.x;      u[2].y = (int)(v.y+dy);  u[2].z = (int)v.z;
  u[3].x = (int)(v.x+dx); u[3].y = (int)(v.y+dy);  u[3].z = (int)v.z;
  u[4].x = (int)v.x;      u[4].y = (int)v.y;       u[4].z = (int)(v.z+dz);
  u[5].x = (int)(v.x+dx); u[5].y = (int)v.y;       u[5].z = (int)(v.z+dz);
  u[6].x = (int)v.x;      u[6].y = (int)(v.y+dy);  u[6].z = (int)(v.z+dz);
  u[7].x = (int)(v.x+dx); u[7].y = (int)(v.y+dy);  u[7].z = (int)(v.z+dz);
  
  for (i=0; i < 8; i++) 
  {
    if (ValidVoxel(img, u[i])) 
      p[i] = GetVoxelIndex(img, u[i]);
    else  //It any neighbour is outside the image domain
    {
      u[0].x = ROUND(v.x);
      u[0].y = ROUND(v.y);
      u[0].z = ROUND(v.z);
      p[0]   = GetVoxelIndex(img, u[0]);
      //Return the closest neighbour
      return img->val[p[0]];
    }    
  }

  //Ref: http://cl.ly/image/0m070b0k1F14
  value_array[0] = (float)img->val[p[1]]*(v.x-u[0].x)+(float)img->val[p[0]]*(u[1].x-v.x);
  value_array[1] = (float)img->val[p[3]]*(v.x-u[2].x)+(float)img->val[p[2]]*(u[3].x-v.x);
  value_array[2] = (float)img->val[p[5]]*(v.x-u[4].x)+(float)img->val[p[4]]*(u[5].x-v.x);
  value_array[3] = (float)img->val[p[7]]*(v.x-u[6].x)+(float)img->val[p[6]]*(u[7].x-v.x);
  value_array[4] = value_array[1]*(v.y-u[0].y) + value_array[0]*(u[2].y-v.y);
  value_array[5] = value_array[3]*(v.y-u[0].y) + value_array[2]*(u[2].y-v.y);
  value  = (int)(value_array[5]*(v.z-u[0].z) + value_array[4]*(u[4].z-v.z) + 0.5);

  return value;
}
Voxel LinearInterpolationCoord(Image *img, FVoxel v)
{
  Voxel u;
  u.x = ROUND(v.x);
  u.y = ROUND(v.y);
  u.z = ROUND(v.z);
  return u;
}
FloatList * IntensityProfile(Image *img, Voxel p1, Voxel pn)
{
  FVoxelList *line = DDAAlgorithm(p1, pn);
  FloatList *list = CreateFloatList(line->n);
  int i = 0;

  for(i=0; i<line->n; i++)
    list->val[i] = LinearInterpolationValue(img, line->val[i]);

  DestroyFVoxelList(line);
  return list;
}
void DrawLine(Image *img, Voxel p1, Voxel pn, int color)
{
  FVoxelList *line = DDAAlgorithm(p1, pn);
  AdjRel *A = Spheric(3);
  Voxel v, u;
  int i=0, j=0, q=0;

  for(i=0; i<line->n; i++)
  {
    u = LinearInterpolationCoord(img, line->val[i]);
    for (j=0; j < A->n; j++) 
    {
      v = GetAdjacentVoxel(A, u, j);
      if (ValidVoxel(img, v))
      {
        q = GetVoxelIndex(img, v);
        img->val[q] = color;
      }
    }
  }
  DestroyFVoxelList(line);
}

CubeFaces* LoadCubeFaces(Image *I) 
{
  int Nx = I->xsize;
  int Ny = I->ysize;
  int Nz = I->zsize;
  int i;

  CubeFaces *cf = (CubeFaces*) malloc(sizeof(CubeFaces)*6);
  for (i=0; i<6; i++)
  {
    cf[i].orthogonal = CreateMatrix(1, 4);
    cf[i].center = CreateMatrix(1, 4);
  }

  // Face of Plane XY
  cf[0].orthogonal->val[AXIS_X] = 0;
  cf[0].orthogonal->val[AXIS_Y] = 0;
  cf[0].orthogonal->val[AXIS_Z] = -1;
  cf[0].orthogonal->val[AXIS_H] = 1;

  cf[0].center->val[AXIS_X] = Nx/2;
  cf[0].center->val[AXIS_Y] = Ny/2;
  cf[0].center->val[AXIS_Z] = 0;
  cf[0].center->val[AXIS_H] = 1;

  // Face of Plane XZ
  cf[1].orthogonal->val[AXIS_X] = 0;
  cf[1].orthogonal->val[AXIS_Y] = -1;
  cf[1].orthogonal->val[AXIS_Z] = 0;
  cf[1].orthogonal->val[AXIS_H] = 1;

  cf[1].center->val[AXIS_X] = Nx/2;
  cf[1].center->val[AXIS_Y] = 0;
  cf[1].center->val[AXIS_Z] = Nz/2;
  cf[1].center->val[AXIS_H] = 1;

  // Face of Plane YZ
  cf[2].orthogonal->val[AXIS_X] = -1;
  cf[2].orthogonal->val[AXIS_Y] = 0;
  cf[2].orthogonal->val[AXIS_Z] = 0;
  cf[2].orthogonal->val[AXIS_H] = 1;

  cf[2].center->val[AXIS_X] = 0;
  cf[2].center->val[AXIS_Y] = Ny/2;
  cf[2].center->val[AXIS_Z] = Nz/2;
  cf[2].center->val[AXIS_H] = 1;

  // Face of Opposite Plane XY
  cf[3].orthogonal->val[AXIS_X] = 0;
  cf[3].orthogonal->val[AXIS_Y] = 0;
  cf[3].orthogonal->val[AXIS_Z] = 1;
  cf[3].orthogonal->val[AXIS_H] = 1;

  cf[3].center->val[AXIS_X] = Nx/2;
  cf[3].center->val[AXIS_Y] = Ny/2;
  cf[3].center->val[AXIS_Z] = Nz-1;
  cf[3].center->val[AXIS_H] = 1;

  // Face of Opposite Plane XZ
  cf[4].orthogonal->val[AXIS_X] = 0;
  cf[4].orthogonal->val[AXIS_Y] = 1;
  cf[4].orthogonal->val[AXIS_Z] = 0;
  cf[4].orthogonal->val[AXIS_H] = 1;

  cf[4].center->val[AXIS_X] = Nx/2;
  cf[4].center->val[AXIS_Y] = Ny-1;
  cf[4].center->val[AXIS_Z] = Nz/2;
  cf[4].center->val[AXIS_H] = 1;

  // Face of Opposite Plane YZ
  cf[5].orthogonal->val[AXIS_X] = 1;
  cf[5].orthogonal->val[AXIS_Y] = 0;
  cf[5].orthogonal->val[AXIS_Z] = 0;
  cf[5].orthogonal->val[AXIS_H] = 1;

  cf[5].center->val[AXIS_X] = Nx-1;
  cf[5].center->val[AXIS_Y] = Ny/2;
  cf[5].center->val[AXIS_Z] = Nz/2;
  cf[5].center->val[AXIS_H] = 1;

  return cf;
}

void ComputeIntersection(Matrix *Tpo, Image *img, Matrix *Tn, CubeFaces *cf, int *p1, int *pn)
{
  int i;
  Matrix *Fi, *CiTpo;
  float lambda[6] = {-1};
  float FiDotTn=0, FiDotCiTpo=0;

  Fi = CreateMatrix(1, 3);
  CiTpo = CreateMatrix(1, 3);
  // Computing the lambda for each face
  for (i = 0; i < 6; i++) {
    Fi->val[AXIS_X] = cf[i].orthogonal->val[AXIS_X];
    Fi->val[AXIS_Y] = cf[i].orthogonal->val[AXIS_Y];
    Fi->val[AXIS_Z] = cf[i].orthogonal->val[AXIS_Z];

    FiDotTn = MatrixDot(Fi, Tn);

    if (FiDotTn != 0) // Fi not orthogonal to Tn
    {
      CiTpo->val[AXIS_X] = cf[i].center->val[AXIS_X] - Tpo->val[AXIS_X];
      CiTpo->val[AXIS_Y] = cf[i].center->val[AXIS_Y] - Tpo->val[AXIS_Y];
      CiTpo->val[AXIS_Z] = cf[i].center->val[AXIS_Z] - Tpo->val[AXIS_Z];

      FiDotCiTpo = MatrixDot(Fi, CiTpo);

      lambda[i] = (float) FiDotCiTpo/FiDotTn;
    }
  }

  *p1 = -1;
  *pn = -1;
  for (int i = 0; i < 6; i++) 
  {
    if (lambda[i] != -1) 
    {
      // Find the voxel of the face[i]
      Voxel v;
      v.x = ROUND(Tpo->val[AXIS_X] + lambda[i] * Tn->val[AXIS_X]);
      v.y = ROUND(Tpo->val[AXIS_Y] + lambda[i] * Tn->val[AXIS_Y]);
      v.z = ROUND(Tpo->val[AXIS_Z] + lambda[i] * Tn->val[AXIS_Z]);

      if (ValidVoxel(img, v)) 
      {
        if (*p1 == -1)
          *p1 = GetVoxelIndex(img, v);
        else if (*p1 != -1) 
        {
          // if p1 and pn is not a vertex
          if (*p1 != GetVoxelIndex(img, v)) 
          {
            *pn = GetVoxelIndex(img, v);
          }
        }
      }
    }
  }

  // if p1 and pn belong to a same vertex
  if ((*p1 != -1) && (*pn == 1)) 
  {
    *pn = *p1;
  }

  if (*p1 > *pn) 
  {
    int p_aux = *p1;
    *p1 = *pn;
    *pn = p_aux;
  }

  DestroyMatrix(Fi);
  DestroyMatrix(CiTpo);
}


Image* MaximumIntensityProfile(Image *img, float xtheta, float ytheta, float ztheta)
{
  float intensity;
  int diagonal=0, p=0, i=0;
  int Nu, Nv;
  int p1, pn;

  CubeFaces *cf;
  Voxel p0, v1, vn;
  FloatList *LineValues;
  Matrix *Mt1, *Mt2, *Mrx, *Mry, *Mrz, *Mtemp, *T;
  Matrix *Norigin, *Nend, *Tnorigin, *Tnend, *Tn;
  Matrix *Tpo;

  //The worst case is when the plane is in the diagonal of the scene.
  diagonal = ROUND(sqrt((double) (img->xsize*img->xsize)+(img->ysize*img->ysize)+(img->zsize*img->zsize)));
  Nu = Nv = diagonal;

  //Creating the output image of the largest size possible, avoids clipping.
  Image *output = CreateImage(Nu, Nv, 1);

  //Compute the final transformation matrix T by multiplying the other matrices
  Mt1 = TranslationMatrix(-Nu/2, -Nv/2, -diagonal);  
  Mt2 = TranslationMatrix(img->xsize/2, img->ysize/2, img->zsize/2);
  Mrx = RotationMatrix('X', xtheta);
  Mry = RotationMatrix('Y', ytheta);
  Mrz = RotationMatrix('Z', ztheta);

  //Computation of the final transformation matrix
  Mtemp = MatrixMultiply(Mrx, Mt1);
  T     = MatrixMultiply(Mry, Mtemp);
  DestroyMatrix(Mtemp);
  Mtemp = T;
  T     = MatrixMultiply(Mrz, Mtemp);
  DestroyMatrix(Mtemp);
  Mtemp = T;
  T     = MatrixMultiply(Mt2, Mtemp);
  DestroyMatrix(Mtemp);
  //End of the computation of the final transformation matrix

  //Normal vector is multiplied by T.
  Norigin =  CreateMatrix(4, 1);
  Norigin->val[0] = Norigin->val[1] = Norigin->val[2] = 0;
  Norigin->val[3] = 1.0;
  Tnorigin = MatrixMultiply(T, Norigin);

  //Normal vector at the end is multiplied by T.
  Nend    =  CreateMatrix(4, 1);
  Nend->val[0] = Nend->val[1] = 0;
  Nend->val[2] = Nend->val[3] = 1;
  Tnend    = MatrixMultiply(T, Nend);

  //A matrix of the differences of Tnend and Tnorigin
  Tn = CreateMatrix(1, 3);
  Tn->val[0] = Tnend->val[0] - Tnorigin->val[0];
  Tn->val[1] = Tnend->val[1] - Tnorigin->val[1];
  Tn->val[2] = Tnend->val[2] - Tnorigin->val[2];

  //Compute a orthogonal vector for each face of the cube and a point in the center of this vector.
  cf = LoadCubeFaces(img);
  for(p=0; p<output->n; p++)
  {
    p1 = pn = -1;
    p0 = GetVoxelCoord(output, p);
    Mtemp = VoxelToMatrix(p0);
    Tpo = MatrixMultiply(T, Mtemp);

    ComputeIntersection(Tpo, img, Tn, cf, &p1, &pn);

    if ((p1 != -1) && (pn != -1)) 
    { 
      v1 = GetVoxelCoord(img, p1);
      vn = GetVoxelCoord(img, pn);

      LineValues = IntensityProfile(img, v1, vn);
      intensity = LineValues->val[0];
      for (i=1; i<LineValues->n; i++)
        if (LineValues->val[i] > intensity)
          intensity = LineValues->val[i];
      output->val[p] = intensity;
    }
    else
      output->val[p] = 0;

    DestroyMatrix(Mtemp);
    DestroyMatrix(Tpo);
  }

  DestroyMatrix(Mt1);
  DestroyMatrix(Mt2);
  DestroyMatrix(Mrx);
  DestroyMatrix(Mry);
  DestroyMatrix(Mrz);
  DestroyMatrix(T);
  DestroyMatrix(Norigin);
  DestroyMatrix(Nend);
  DestroyMatrix(Tnorigin);
  DestroyMatrix(Tnend);
  DestroyMatrix(Tn);
  DestroyCubeFaces(cf);
  return output;
}