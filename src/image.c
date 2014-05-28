#include "image.h"

/*******************


PRIVATE FUNCTIONS


******************/

VolumeFaces* CreateVolumeFaces(Image *I)
{
  int Nx = I->xsize;
  int Ny = I->ysize;
  int Nz = I->zsize;
  int i;

  VolumeFaces *vf = (VolumeFaces*) malloc(sizeof(VolumeFaces)*6);
  
  for (i=0; i<6; i++)
  {
    vf[i].orthogonal = CreateMatrix(1, 4);
    vf[i].center = CreateMatrix(1, 4);
  }

  // Face of Plane XY
  vf[0].orthogonal->val[AXIS_X] = 0;
  vf[0].orthogonal->val[AXIS_Y] = 0;
  vf[0].orthogonal->val[AXIS_Z] = -1;
  vf[0].orthogonal->val[AXIS_H] = 1;

  vf[0].center->val[AXIS_X] = Nx/2;
  vf[0].center->val[AXIS_Y] = Ny/2;
  vf[0].center->val[AXIS_Z] = 0;
  vf[0].center->val[AXIS_H] = 1;

  // Face of Plane XZ
  vf[1].orthogonal->val[AXIS_X] = 0;
  vf[1].orthogonal->val[AXIS_Y] = -1;
  vf[1].orthogonal->val[AXIS_Z] = 0;
  vf[1].orthogonal->val[AXIS_H] = 1;

  vf[1].center->val[AXIS_X] = Nx/2;
  vf[1].center->val[AXIS_Y] = 0;
  vf[1].center->val[AXIS_Z] = Nz/2;
  vf[1].center->val[AXIS_H] = 1;

  // Face of Plane YZ
  vf[2].orthogonal->val[AXIS_X] = -1;
  vf[2].orthogonal->val[AXIS_Y] = 0;
  vf[2].orthogonal->val[AXIS_Z] = 0;
  vf[2].orthogonal->val[AXIS_H] = 1;

  vf[2].center->val[AXIS_X] = 0;
  vf[2].center->val[AXIS_Y] = Ny/2;
  vf[2].center->val[AXIS_Z] = Nz/2;
  vf[2].center->val[AXIS_H] = 1;

  // Face of Opposite Plane XY
  vf[3].orthogonal->val[AXIS_X] = 0;
  vf[3].orthogonal->val[AXIS_Y] = 0;
  vf[3].orthogonal->val[AXIS_Z] = 1;
  vf[3].orthogonal->val[AXIS_H] = 1;

  vf[3].center->val[AXIS_X] = Nx/2;
  vf[3].center->val[AXIS_Y] = Ny/2;
  vf[3].center->val[AXIS_Z] = Nz-1;
  vf[3].center->val[AXIS_H] = 1;

  // Face of Opposite Plane XZ
  vf[4].orthogonal->val[AXIS_X] = 0;
  vf[4].orthogonal->val[AXIS_Y] = 1;
  vf[4].orthogonal->val[AXIS_Z] = 0;
  vf[4].orthogonal->val[AXIS_H] = 1;

  vf[4].center->val[AXIS_X] = Nx/2;
  vf[4].center->val[AXIS_Y] = Ny-1;
  vf[4].center->val[AXIS_Z] = Nz/2;
  vf[4].center->val[AXIS_H] = 1;

  // Face of Opposite Plane YZ
  vf[5].orthogonal->val[AXIS_X] = 1;
  vf[5].orthogonal->val[AXIS_Y] = 0;
  vf[5].orthogonal->val[AXIS_Z] = 0;
  vf[5].orthogonal->val[AXIS_H] = 1;

  vf[5].center->val[AXIS_X] = Nx-1;
  vf[5].center->val[AXIS_Y] = Ny/2;
  vf[5].center->val[AXIS_Z] = Nz/2;
  vf[5].center->val[AXIS_H] = 1;

  return vf;
}
void       DestroyVolumeFaces(VolumeFaces *cf)
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

/*******************


 PUBLIC FUNCTIONS


******************/
void CopyVoxelSize(Image *img1, Image *img2)
{
  img2->dx = img1->dx;
  img2->dy = img1->dy;
  img2->dz = img1->dz;
}

char ValidVoxel(Image *img, Voxel u)
{
  if ((u.x >= 0)&&(u.x < img->xsize)&&
      (u.y >= 0)&&(u.y < img->ysize)&&
      (u.z >= 0)&&(u.z < img->zsize) )
    return(1);
  else
    return(0);
}
char FValidVoxel(Image *img, FVoxel u)
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
  img->Cb = img->Cr = NULL;
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
Image *CreateColorImage(int xsize, int ysize, int zsize)
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
  img->Cb    = AllocIntArray(img->n);
  img->Cr    = AllocIntArray(img->n);
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
FImage *CreateFImage(int xsize, int ysize, int zsize)
{
  FImage *img=NULL;
  int    i;

  img = (FImage *) calloc(1,sizeof(FImage));
  img->n     = xsize*ysize*zsize;
  img->xsize = xsize; 
  img->ysize = ysize; 
  img->zsize = zsize; 
  img->dx    = img->dy = img->dz = 1.0;
  img->val   = AllocFloatArray(img->n);
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
void DestroyFImage(FImage *img)
{

  if(img != NULL){
    free(img->val); 
    free(img->tby);
    free(img->tbz);
    free(img);
  }

}

void DestroyImage(Image *img)
{

  if(img != NULL){
    free(img->val); 
    if (img->Cb != NULL) free(img->Cb);
    if (img->Cr != NULL) free(img->Cr);
    free(img->tby);
    free(img->tbz);
    free(img);
  }

}
Color RGBtoYCbCr(Color cin)
{
  Color cout;

  cout.val[0]=(int)(0.257*(float)cin.val[0]+
        0.504*(float)cin.val[1]+
        0.098*(float)cin.val[2]+16.0);
  cout.val[1]=(int)(-0.148*(float)cin.val[0]+
        -0.291*(float)cin.val[1]+
        0.439*(float)cin.val[2]+128.0);
  cout.val[2]=(int)(0.439*(float)cin.val[0]+
        -0.368*(float)cin.val[1]+
        -0.071*(float)cin.val[2]+128.0);
   
  for(int i=0; i < 3; i++) {
    if (cout.val[i]<0)   cout.val[i]=0; 
    if (cout.val[i]>255) cout.val[i]=255; 
  }

  return(cout);
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
// void WriteImageP6(Image *img, char *filename) 
// {
//   FILE *fp=NULL;
//   int   p;
//   Color YCbCr,RGB;

//   fp = fopen(filename,"w"); 
//   if (fp == NULL) 
//     Error("Cannot open file in the specified path","WriteImageP6");

//   fprintf(fp,"P6\n");
//   fprintf(fp,"%d %d\n",img->xsize,img->ysize);

//   img->maxval = MaximumValue(img);
//   img->minval = MinimumValue(img);

//   if ((img->maxval < 256)&&(img->minval>=0)) {
//     fprintf(fp,"%d\n",255);
//     for (p=0; p < img->n; p++) {
//       YCbCr.val[0] = img->val[p];
//       YCbCr.val[1] = img->Cb[p];
//       YCbCr.val[2] = img->Cr[p];

//       RGB = YCbCrtoRGB(YCbCr);

//       fputc(((uchar)RGB.val[0]),fp);
//       fputc(((uchar)RGB.val[1]),fp);
//       fputc(((uchar)RGB.val[2]),fp);        
//     }
//   }else {
//     Error("Cannot write image as P6","WriteImageP6");
//   }
//   fclose(fp);
// }

/* ----------------- Voxel-based algorithms ---------------------*/
Image *Normalize(Image *img, float minval, float maxval)
{
  Image *nimg = CreateImage(img->xsize, img->ysize, img->zsize);
  int p;

  // if (img->Cb != NULL) 
  //   CopyCbCr(img,nimg);
  CopyVoxelSize(img,nimg);

  img->minval  = MinimumValue(img);
  img->maxval  = MaximumValue(img);
  nimg->maxval = maxval;
  nimg->minval = minval;

  if (img->minval < img->maxval){
    for (p=0; p < img->n; p++) 
      nimg->val[p] = (int)((maxval-minval)*((float)img->val[p]-(float)img->minval)/((float)img->maxval-(float)img->minval) + minval);
  }else{
    Error("Cannot normalize empty image","Normalize");
  }
  return(nimg);
}

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



/* ----------------- Task 3 ---------------------*/
//Given a transformed ray (Tpo), a volume, a normal vector transformed and vectores orthogonal to each volume face, return the points of intersection.
int ComputeIntersection(Matrix *Tpo, Image *img, Matrix *Tn, VolumeFaces *vf, int *p1, int *pn)
{
  int i, p_aux;
  Matrix *Fi, *CiTpo;
  float lambda[6] = {-1};
  float FiDotTn=0, FiDotCiTpo=0;
  Voxel v;
  *p1 = *pn = -1;

  Fi = CreateMatrix(1, 3);
  CiTpo = CreateMatrix(1, 3);
  // Computing the lambda for each face
  for (i = 0; i < 6; i++) 
  {
    Fi->val[AXIS_X] = vf[i].orthogonal->val[AXIS_X];
    Fi->val[AXIS_Y] = vf[i].orthogonal->val[AXIS_Y];
    Fi->val[AXIS_Z] = vf[i].orthogonal->val[AXIS_Z];

    FiDotTn = MatrixInnerProduct(Fi, Tn);

    if (FiDotTn != 0) // vf[i].orthogonal not orthogonal to Tn
    {
      CiTpo->val[AXIS_X] = vf[i].center->val[AXIS_X] - Tpo->val[AXIS_X];
      CiTpo->val[AXIS_Y] = vf[i].center->val[AXIS_Y] - Tpo->val[AXIS_Y];
      CiTpo->val[AXIS_Z] = vf[i].center->val[AXIS_Z] - Tpo->val[AXIS_Z];

      FiDotCiTpo = MatrixInnerProduct(Fi, CiTpo);

      lambda[i] = (float) FiDotCiTpo/FiDotTn;

      v.x = ROUND(Tpo->val[AXIS_X] + lambda[i] * Tn->val[AXIS_X]);
      v.y = ROUND(Tpo->val[AXIS_Y] + lambda[i] * Tn->val[AXIS_Y]);
      v.z = ROUND(Tpo->val[AXIS_Z] + lambda[i] * Tn->val[AXIS_Z]);

      if (ValidVoxel(img, v)) 
      {
        if (*p1 == -1)
          *p1 = GetVoxelIndex(img, v);
        else if (*p1 != -1) // if p1 and pn is not a vertex
        {
          if (*p1 != GetVoxelIndex(img, v)) 
            *pn = GetVoxelIndex(img, v);
        }
      }
    }
  }

  // if p1 and pn belong to a same vertex
  if ((*p1 != -1) && (*pn == 1)) 
    *pn = *p1;

  if (*p1 > *pn)
  {
    p_aux = *p1;
    *p1 = *pn;
    *pn = p_aux;
  }

  DestroyMatrix(Fi);
  DestroyMatrix(CiTpo);
  if ((*p1 != -1) && (*pn != -1))
    return 1;
  else
    return 0;
}


Image* MaximumIntensityProjection(Image *img, float xtheta, float ytheta, float ztheta)
{
  float intensity;
  int diagonal=0, p=0, i=0;
  int Nu, Nv;
  int p1, pn;

  VolumeFaces *vf;
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

  T = ComputeTransformation(5, Mt1, Mrx, Mry, Mrz, Mt2);
  //Computation of the final transformation matrix
  // Mtemp = MatrixMultiply(Mrx, Mt1);
  // T     = MatrixMultiply(Mry, Mtemp);
  // DestroyMatrix(Mtemp);
  // Mtemp = T;
  // T     = MatrixMultiply(Mrz, Mtemp);
  // DestroyMatrix(Mtemp);
  // Mtemp = T;
  // T     = MatrixMultiply(Mt2, Mtemp);
  // DestroyMatrix(Mtemp);
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
  vf = CreateVolumeFaces(img);
  for(p=0; p<output->n; p++)
  {
    p1 = pn = -1;
    p0 = GetVoxelCoord(output, p);
    Mtemp = VoxelToMatrix(p0);
    Tpo = MatrixMultiply(T, Mtemp);

    

    if (ComputeIntersection(Tpo, img, Tn, vf, &p1, &pn))
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
  DestroyVolumeFaces(vf);
  return output;
}









/* ----------------- Task 4 ---------------------*/
float PhongShadingTask4(int p, float distance, float diagonal, Matrix *N, Matrix *ObserverVector)
{
  float cos_theta;
  float cos_2theta, power, phong_val=0.0;
  float ka, kd, ks, ra, ns;
  float D = (int)(2 * diagonal + 1);
  ka = 0.1;
  kd = 0.7;
  ks = 0.2;
  ra = 255;
  ns = 5.0;

  cos_theta  = MatrixInnerProduct(ObserverVector, N);
  if (cos_theta > 1.0) 
    cos_theta=1.0;
  
  if (cos_theta > Epsilon) 
  {  
    cos_2theta = 2*cos_theta*cos_theta - 1;

    if (cos_2theta <= Epsilon)
      power = 0.;
    else 
    {
      power = 1.;
      for (int k=0; k < ns; k++)
        power = power*cos_2theta;
    }
    phong_val = ka*ra + 255*((distance-D)/(0-D))*(kd*cos_theta + ks * power);
  }
  return phong_val;
}

void ImageGradientMagnitudeAndIndex(Image *img, Image *gradientImg, Image *normalIndexImg, AdjRel *A)
{
  float   dist, gx , gy, gz, gModule, alfa, gamma;
  int     i, p, q;
  Voxel   u,v;
  Vector  unitNormalVector;
  float   *mag = AllocFloatArray(A->n);

  CopyVoxelSize(img, gradientImg);
  CopyVoxelSize(img, normalIndexImg);

  for (i=0; i < A->n; i++)
      mag[i]=sqrtf(A->adj[i].dx*A->adj[i].dx+A->adj[i].dy*A->adj[i].dy+A->adj[i].dz*A->adj[i].dz);

   for (u.z=0; u.z < img->zsize; u.z++)
      for (u.y=0; u.y < img->ysize; u.y++)
         for (u.x=0; u.x < img->xsize; u.x++) {
      p = GetVoxelIndex(img,u);
      gx = gy = gz = 0.0;
        for (i=1; i < A->n; i++) {
         v.x = u.x + A->adj[i].dx;
         v.y = u.y + A->adj[i].dy;
         v.z = u.z + A->adj[i].dz;
         if (ValidVoxel(img,v)){
            q = GetVoxelIndex(img,v);
            dist = img->val[q]-img->val[p];
                  gx  += dist*A->adj[i].dx/mag[i];
            gy  += dist*A->adj[i].dy/mag[i];
            gz  += dist*A->adj[i].dz/mag[i];
         }              
      }
      gModule = sqrtf(gx*gx + gy*gy + gz*gz);
      gradientImg->val[p]= ROUND(gModule);   
      if (gModule > 0.0)     
      {
         unitNormalVector.x = (-1) * gx / gModule;
         unitNormalVector.y = (-1) * gy / gModule;
         unitNormalVector.z = (-1) * gz / gModule;
         gamma = asin(unitNormalVector.z) * 180 / PI;
         if (unitNormalVector.x != 0.0)
            alfa = atan(unitNormalVector.y / unitNormalVector.x) * 180 / PI;
         else
           alfa = 0.0;
         if (alfa < 0)
            alfa = alfa + 360; 
         normalIndexImg->val[p] = ROUND(360 * (gamma + 90) + alfa + 1);
      }
      else
      {
         normalIndexImg->val[p] = 0;
      }
   }

  free(mag);
}

FImage *CreateOpacityImage(Image *img, Image *gradient)
{
  FImage *output = CreateFImage(img->xsize, img->ysize, img->zsize);

  int p;
  int g = 1000;
  float l1s =  900, l2s = 950, l3s = 1050, l4s = 1100, alfaTs = 2.0, gs = 1000;
  float l1b = 1101, l2b = 1200, l3b = 2400, l4b = 2500, alfaTb = 1.0, gb = 1000;
  int intensity;
  float ogp;
  float oip;

  for (p=0; p<output->n; p++)
  {
    intensity = img->val[p];
    
    ogp = 1.0 / (1+exp(-gradient->val[p] + g));

    //Skin
    if (intensity <= l1s)
      oip = 0.0;
    if ((intensity > l1s) && (intensity <= l2s))
      oip = (intensity - l1s) / (l2s - l1s);
    if ((intensity > l2s) && (intensity <= l3s))
      oip = alfaTs;
    if ((intensity > l3s) && (intensity <= l4s))
      oip = (l4s - intensity) / (l4s - l3s);

    //Bone
    if ((intensity > l4s) && (intensity <= l1b))
      oip = 0.0;
    if ((intensity > l1b) && (intensity <= l2b))
      oip = (intensity - l1b) / (l2b - l1b);
    if ((intensity > l2b) && (intensity <= l3b))
      oip = alfaTb;
    if ((intensity > l3b) && (intensity <= l4b))
      oip = (l4b - intensity) / (l4b - l3b);

    if (intensity > l4s)
      g = gb;
    else
      g = gs;

    output->val[p] = ogp * oip;
  }
  return output;
}
Matrix *CreateNormalLookUpMatrix()
{
  //65160 = 361*180. Number of iterations below
  Matrix *table = CreateMatrix(65160, 3);

  int row = 0;
  int gamma, alpha;
  float gamma_rad, alpha_rad;
  for (gamma=-90; gamma <= 90; gamma++)
  {
    gamma_rad = (PI*gamma)/180.0;
    for (alpha=0; alpha <= 359; alpha++)
    {
      alpha_rad = (PI * alpha)/180.0;
      table->val[GetMatrixIndex(table, row, AXIS_X)] = cos(gamma_rad)*cos(alpha_rad);
      table->val[GetMatrixIndex(table, row, AXIS_Y)] = cos(gamma_rad)*sin(alpha_rad);
      table->val[GetMatrixIndex(table, row, AXIS_Z)] = sin(gamma_rad);
      row++;
    }
  }
  //printf("Number of rows: %d\n", row);
  return table;
}
int VolumeRenderValue(Voxel p0, Voxel p1, Voxel pn, Image *scene, Image *normalIndexImg, Matrix *normalTable, Matrix*ObserverVector, FImage *opacity)
{
  FVoxelList *FVoxels = DDAAlgorithm(p1, pn);
  Voxel v;
  int i=0, p=0;
  float acc_opacity = 1.0, distance=0;
  float opac=0, phong_val=0;
  float intensity = 0;
  float diagonal = ROUND(sqrt((double) (scene->xsize*scene->xsize)+(scene->ysize*scene->ysize)+(scene->zsize*scene->zsize)));
  Matrix *N = CreateMatrix(1, 3);

  for (i=0; (i < FVoxels->n) && (acc_opacity > Epsilon); i++)
  {
    if (FValidVoxel(scene, FVoxels->val[i]))
    {
      //Zero Degree Interpolation
      v.x = ROUND(FVoxels->val[i].x); 
      v.y = ROUND(FVoxels->val[i].y);
      v.z = ROUND(FVoxels->val[i].z);

      p = GetVoxelIndex(scene, v);
      if (opacity->val[p] > 0)
      {
        distance = VoxelDistance(v, p0);
        N->val[AXIS_X] = normalTable->val[GetMatrixIndex(normalTable, normalIndexImg->val[p], AXIS_X)];
        N->val[AXIS_Y] = normalTable->val[GetMatrixIndex(normalTable, normalIndexImg->val[p], AXIS_Y)];
        N->val[AXIS_Z] = normalTable->val[GetMatrixIndex(normalTable, normalIndexImg->val[p], AXIS_Z)];
        opac = opacity->val[p];
        phong_val = opac * PhongShadingTask4(p, distance, diagonal, N, ObserverVector) * acc_opacity;
        intensity = phong_val * scene->val[p];
        acc_opacity = acc_opacity * (1.0 -  opac);
      }

    }
  }
  return ROUND(intensity);
}
Image *RayCasting(Image *img, float xtheta, float ytheta, float ztheta)
{
  int diagonal=0, p=0, Nu, Nv;
  int p1, pn;

  AdjRel *radius1 = Spheric(1);
  FImage *opacity;
  Matrix *normalTable;

  VolumeFaces *vf;
  Voxel p0, v1, vn;
  Matrix *Mt1, *Mt2, *Mrx, *Mry, *Mrz, *Mtemp, *T;
  Matrix *Norigin, *Nend, *Tnorigin, *Tnend, *Tn;
  Matrix *Tpo;
  Matrix *Result;
  Matrix *ObserverVector = CreateMatrix(1, 3);


  //The worst case is when the plane is in the diagonal of the scene.
  diagonal = ROUND(sqrt((double) (img->xsize*img->xsize)+(img->ysize*img->ysize)+(img->zsize*img->zsize)));
  Nu = Nv = diagonal;
  Image *output = CreateImage(Nu, Nv, 1);
  Image *normalized;

  //Compute the final transformation matrix T by multiplying the other matrices
  Mt1 = TranslationMatrix(-Nu/2, -Nv/2, -diagonal);  
  Mt2 = TranslationMatrix(img->xsize/2, img->ysize/2, img->zsize/2);
  Mrx = RotationMatrix('X', xtheta);
  Mry = RotationMatrix('Y', ytheta);
  Mrz = RotationMatrix('Z', ztheta);

  T = ComputeTransformation(5, Mt1, Mrx, Mry, Mrz, Mt2);

  //Computing the observer vector
  Mtemp = CreateMatrix(1, 4);
  Mtemp->val[AXIS_X] = Mtemp->val[AXIS_Y] = 0;
  Mtemp->val[AXIS_Z] = -1;
  Mtemp->val[AXIS_H] = 1;
  Result = ComputeTransformation(3, Mrx, Mry, Mtemp);
  ObserverVector->val[AXIS_X] = Result->val[AXIS_X];
  ObserverVector->val[AXIS_Y] = Result->val[AXIS_Y];
  ObserverVector->val[AXIS_Z] = Result->val[AXIS_Z];
  DestroyMatrix(Mtemp);

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
  Image *gradientImg = CreateImage(img->xsize,img->ysize,img->zsize);
  Image *normalIndexImg = CreateImage(img->xsize,img->ysize,img->zsize);
  ImageGradientMagnitudeAndIndex(img, gradientImg, normalIndexImg, radius1);
  vf = CreateVolumeFaces(img);
  opacity = CreateOpacityImage(img, gradientImg);
  normalTable = CreateNormalLookUpMatrix();
  for(p=0; p<output->n; p++)
  {
    p1 = pn = -1;

    p0 = GetVoxelCoord(output, p); //Since the plane is 2d, Z axis will always be zero.
    p0.z = -diagonal;              //Put z in -diagonal.
    Mtemp = VoxelToMatrix(p0);
    Tpo = MatrixMultiply(T, Mtemp);

    if (ComputeIntersection(Tpo, img, Tn, vf, &p1, &pn)) //If the ray intersects the image
    { 
      v1 = GetVoxelCoord(img, p1);
      vn = GetVoxelCoord(img, pn);
      
      int val = VolumeRenderValue(p0, v1, vn, img, normalIndexImg, normalTable, ObserverVector, opacity);
      output->val[p] = val;
    }
    else
      output->val[p] = 0;

    DestroyMatrix(Mtemp);
    DestroyMatrix(Tpo);
  }

  normalized = Normalize(output, 0, 255);
  // for (p=0; p<output->n; p++)
  //   if (output->val[p])
  //     printf("%d\n", output->val[p]);

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
  DestroyMatrix(Result);
  DestroyMatrix(ObserverVector);
  DestroyMatrix(normalTable);
  DestroyAdjRel(radius1);
  DestroyFImage(opacity);
  DestroyVolumeFaces(vf);
  DestroyImage(output);
  return normalized;
}