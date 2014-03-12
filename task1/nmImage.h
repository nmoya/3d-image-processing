#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

typedef struct nm_image {
  int    *val;
  int xsize,ysize,zsize;
  float dx,dy,dz;
  int maxval, minval, n; // minimum and maximum values, and number of voxels 
} nmImage;

typedef struct nm_voxel {
	int x, y, z;
}nmVoxel;

nmImage *nmCreateImage(int xsize,int ysize,int zsize);
nmImage *nmReadSCNImage(char filename[]);
nmImage * nmGetAxialSlice(nmImage *img, int slice);
nmImage * nmGetSagitalSlice(nmImage *img, int slice);
nmImage * nmGetCoronalSlice(nmImage *img, int slice);
nmImage * nmLinearStretching(nmImage *img, float level, float width);

void nmDestroyImage(nmImage **img);
void nmWriteImageP2(nmImage *img, char filename[]);
void nmError(char message[], char function[]);

int *nmAllocIntArray(int n); 
int nmVoxelToIndex(nmImage *img, nmVoxel v);
int nmMinimumValue(nmImage *img);
int nmMaximumValue(nmImage *img);

float nmMeanValue(nmImage *img);
float nmStdDevValue(nmImage *img);
float *nmAllocFloatArray(int n);

nmVoxel nmIndexToVoxel(nmImage *img, int index);