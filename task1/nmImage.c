#include "nmImage.h"


nmImage *nmCreateImage(int xsize,int ysize,int zsize)
{
    nmImage *img = NULL;

    img = (nmImage *) malloc(sizeof(nmImage));
    if (!img)
        nmError("Error creating and image", "nmCreateImage");

    img->val = NULL;
    img->val = nmAllocIntArray(xsize * ysize * zsize);
    img->xsize = xsize;
    img->ysize = ysize;
    img->zsize = zsize;
    img->dx      = 1.0;
    img->dy      = 1.0;
    img->dz      = 1.0;
    img->maxval  = 0;
    img->minval  = 0;
    img->n       = xsize*ysize*zsize;

    if (!img->val) 
        nmError("Error creating values array", "nmCreateImage");
    
    return img;
}
void nmDestroyImage(nmImage **img)
{
    nmImage *aux = *img;
    if (aux != NULL)
    {
        free(aux->val);
        aux->val = NULL;
        free(aux);
        aux = NULL;
    }
}
nmImage* nmReadSCNImage(char filename[])
{
    nmImage  *img=NULL;
    FILE    *fp=NULL;
    unsigned char   *data8=NULL;
    unsigned short  *data16=NULL;
    int     *data32=NULL;
    char    type[10];
    int     p,v,xsize,ysize,zsize;

    fp = fopen(filename, "r");
    if (!fp)
        nmError("Error opening file","nmReadSCNmage");

    if (fscanf(fp,"%s\n",type) != 1) 
        nmError("Reading type error","nmReadSCNmage");

    if (strcmp(type,"SCN")!=0) 
        nmError("Invalid file type","nmReadSCNImage");

    if (fscanf(fp,"%d %d %d\n",&xsize,&ysize,&zsize)!=3) 
        nmError("Reading xsize, ysize, zsize error","nmReadSCNImage");

    img = nmCreateImage(xsize, ysize, zsize);
    if (fscanf(fp,"%f %f %f\n",&img->dx,&img->dy,&img->dz)!=3) 
        nmError("Reading dx, dy, dz error","nmReadSCNImage");

    if (fscanf(fp,"%d\n",&v)!=1) 
        nmError("Reading maxval error","nmReadSCNImage");

    if (v==8)
    {
        data8 = (unsigned char *) malloc(img->n * sizeof(unsigned char));
        if (fread(data8, sizeof(unsigned char), img->n, fp) != img->n)
            nmError("Reading data8 error","nmReadSCNImage");
        
        for (p=0; p < img->n; p++)
            img->val[p] = (int) data8[p];
        free(data8);
    }
    else if (v == 16)
    {
        data16 = (unsigned short *) malloc(img->n * sizeof(unsigned short));
        if (fread(data16, sizeof(unsigned short), img->n, fp) != img->n)
            nmError("Reading data16 error","nmReadSCNImage");
        
        for (p=0; p < img->n; p++)
            img->val[p] = (int) data16[p];
        free(data16);
    }
    else if (v == 32)
    {
        data32 = nmAllocIntArray(img->n);
        if (fread(data32, sizeof(int), img->n, fp) != img->n)
            nmError("Reading data32 error","nmReadSCNImage");
        
        for (p=0; p < img->n; p++)
            img->val[p] = data32[p];
        free(data32);
    }
    else
        nmError("Image depth is bigger then expected", "nmReadSCNImage");

    img->minval = nmMinimumValue(img);
    img->maxval = nmMaximumValue(img);
    
    fclose(fp);
    return img;
}

void nmWriteImageP2(nmImage *img, char filename[])
{
    FILE *fp=NULL;
    int   p;

    fp = fopen(filename,"w"); 
    if (!fp) 
        nmError("Failed to open P2 image","nmWriteImageP2");

    fprintf(fp,"P2\n");
    fprintf(fp,"%d %d\n",img->xsize,img->ysize);

    img->maxval = nmMaximumValue(img);

    fprintf(fp,"%d\n",img->maxval);
    for (p=0; p < img->n; p++) 
        fprintf(fp,"%d ",img->val[p]);
    fclose(fp);
}
nmImage * nmGetAxialSlice(nmImage *img, int slice)
{
    int i, index;
    nmImage *tempImg = nmCreateImage(img->xsize, img->ysize, 1);
    index = 0;
    for (i= (img->xsize*img->ysize)*slice; i < (img->xsize*img->ysize)*(slice+1); i++)
    {
        tempImg->val[index] = img->val[i];
        index++;
    }
    return tempImg;
}
nmImage * nmGetCoronalSlice(nmImage *img, int slice)
{
    int i, index;
    nmImage *tempImg = nmCreateImage(img->xsize, img->ysize, 1);
    index = 0;
    for (i= (img->xsize*img->zsize)*slice; i < (img->xsize*img->zsize)*(slice+1); i++)
    {
        tempImg->val[index] = img->val[i];
        index++;
    }
    return tempImg;
}
nmImage * nmGetSagitalSlice(nmImage *img, int slice)
{
    int i, index;
    nmImage *tempImg = nmCreateImage(img->xsize, img->ysize, 1);
    index = 0;
    for (i= (img->ysize*img->zsize)*slice; i < (img->ysize*img->zsize)*(slice+1); i++)
    {
        tempImg->val[index] = img->val[i];
        index++;
    }
    return tempImg;
}

void nmError(char *msg, char *func)
{ 
  fprintf(stderr,"Error in %s: \n%s. \n",func, msg);
  exit(1);
}
int nmMinimumValue(nmImage *img)
{
  int p;       

  img->minval = INT_MAX;
  for (p=0; p < img->n; p++) 
    if (img->minval > img->val[p])
      img->minval = img->val[p];

  return img->minval;
}
int nmMaximumValue(nmImage *img)
{
  int p;       

  img->maxval = -INT_MAX;
  for (p=0; p < img->n; p++) 
    if (img->maxval < img->val[p])
      img->maxval = img->val[p];

  return img->maxval;
}

int *nmAllocIntArray(int n)
{
    int *arr=NULL;
    arr = (int *) malloc(n * sizeof(int));
    if (arr == NULL)
        nmError("Error allocing int array" ,"nmAllocIntArray");
    memset(arr, 0, n*sizeof(int));
    return arr;
}
float *nmAllocFloatArray(int n)
{
    float *arr=NULL;
    arr = (float *) malloc(n * sizeof(float));
    if (arr == NULL)
        nmError("Error allocing float array" ,"nmAllocFloatArray");
    memset(arr, 0, n*sizeof(float));
    return arr;
}