#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


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
void nmDestroyImage(nmImage **img);

nmImage* nmReadSCNImage(char filename[]);
void nmWriteImageP2(nmImage *img, char filename[]);

void nmError(char message[], char function[]);
int *nmAllocIntArray(int n); 
float *nmAllocFloatArray(int n);

int nmMinimumValue(nmImage *img);
int nmMaximumValue(nmImage *img);

nmImage * nmGetAxialSlice(nmImage *img, int slice);
nmImage * nmGetSagitalSlice(nmImage *img, int slice);
nmImage * nmGetCoronalSlice(nmImage *img, int slice);

nmVoxel nmIndexToVoxel(nmImage *img, int index);
int nmVoxelToIndex(nmImage *img, nmVoxel v);