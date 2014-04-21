#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "common.h"
#include "adjacency.h"
#include "kernel.h"

/* 3D grayscale image */

typedef struct _image {
  int    xsize,ysize,zsize; /* number of voxels along x, y, and z */
  int    n, maxval, minval; /* number of voxels */
  int   *val;               /* list of voxel intensities */ 
  float  dx,dy,dz;          /* voxel dimensions in mm along x, y, and z */
  int   *tby, *tbz;         /* look-up tables to speed up conversions between 
			       voxel and its coordinates */
} Image;

/* Voxel and coordinate conversions */

#define GetXCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define GetYCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define GetZCoord(s,p) ((p) / (((s)->xsize)*((s)->ysize)))
#define GetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])

Voxel   GetVoxelCoord(Image *img, int p); 
char    ValidVoxel(Image *img, Voxel u);
int 	VoxelValue(Image *img, Voxel v);

/* Allocate memory to store image */ 
Image        *CreateImage(int xsize, int ysize, int zsize); 
void          DestroyImage(Image *img); /* Free memory */
Image        *ReadImage(char *filename); /* Read image from a file */
void          WriteImage(Image *img, char *filename); /* Write image into a file */ 


/*--------------------- Voxel-based methods---------------------------- */

int           MaximumValue(Image *img); /* Maximum voxel intensity */
int           MinimumValue(Image *img); /* Minimum voxel intensity */
Image        *WindowAndLevel(Image *img, int window, int level, int H); /* Intensity transformation */

/*-------------------- Adjacency-based methods ------------------------- */

Image        *LinearFilter(Image *img, Kernel *K); /* Filter image by convolution */

/*------------------- Connectivity-based methods ------------------------*/

Image        *LabelBinaryImage(Image *img, AdjRel *A); /* Label foreground components */


/*--------------------- Task 1---------------------------- */
Image 		* GetSagitalSlice(Image *img, int slice);
Image 		* GetAxialSlice(Image *img, int slice);
Image 		* GetCoronalSlice(Image *img, int slice);
void 		WriteImageP2(Image *img, char filename[]);

/*--------------------- Task 2---------------------------- */
int 		CompareVoxels(Voxel v1, Voxel v2);
FloatList 	* IntensityProfile(Image *img, Voxel p1, Voxel pn);
FVoxelList  * DDAAlgorithm(Voxel p1, Voxel pn);
int 		LinearInterpolationValue(Image *img, FVoxel v);
Voxel 		LinearInterpolationCoord(Image *img, FVoxel v);
void 		DrawLine(Image *img, Voxel p1, Voxel pn, int color);

/*--------------------- Task 3---------------------------- */
Image* MaximumIntensityProfile(Image *img, float xtheta, float ytheta, float ztheta);

#endif
