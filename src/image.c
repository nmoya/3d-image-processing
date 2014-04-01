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
FVoxel * DDAAlgorithm(Voxel p1, Voxel pn)
{
  int n, k;
  float dx=0, dy=0, dz=0, DX, DY, DZ;
  FVoxel *output = NULL;

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
  output = (FVoxel*)malloc((n+1)*sizeof(FVoxel));
  output[0].x = p1.x;
  output[0].y = p1.y;
  output[0].z = p1.z;

  for (k=1; k<n; k++)
  {
    output[k].x = output[k-1].x + dx;
    output[k].y = output[k-1].y + dy;
    output[k].z = output[k-1].z + dz;
    output[k].terminus = 0;
  }
  output[k].terminus = 1;
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
float * IntensityProfile(Image *img, Voxel p1, Voxel pn)
{
  FVoxel *queue = DDAAlgorithm(p1, pn);
  float *output = NULL;
  int i = 0;

  //Count the number of points in the line
  output = AllocFloatArray(FVoxelSize(queue));

  i=0;
  while(!queue[i].terminus)
  {
    output[i] = LinearInterpolationValue(img, queue[i]);
    printf("%g\n", output[i]);
    i++;
  }

  free(queue);
  return output;
}
void DrawLine(Image *img, Voxel p1, Voxel pn, int color)
{
  FVoxel *queue = DDAAlgorithm(p1, pn);
  Voxel v;
  int i=0;
  while(!queue[i].terminus)
  {
    v = LinearInterpolationCoord(img, queue[i++]);
    img->val[GetVoxelIndex(img, v)] = color;
  }
}