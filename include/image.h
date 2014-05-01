#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "common.h"
#include "adjacency.h"
#include "kernel.h"
#include "matrix.h"

/* 3D grayscale image */

typedef struct _image {
  int    xsize,ysize,zsize; /* number of voxels along x, y, and z */
  int    n, maxval, minval; /* number of voxels */
  int   *val;               /* list of voxel intensities */ 
  float  dx,dy,dz;          /* voxel dimensions in mm along x, y, and z */
  int   *tby, *tbz;         /* look-up tables to speed up conversions between 
			       voxel and its coordinates */
} Image;

typedef struct _fimage {
  int    xsize,ysize,zsize; /* number of voxels along x, y, and z */
  int    n, maxval, minval; /* number of voxels */
  float   *val;               /* list of voxel intensities */ 
  float  dx,dy,dz;          /* voxel dimensions in mm along x, y, and z */
  int   *tby, *tbz;         /* look-up tables to speed up conversions between 
			       voxel and its coordinates */
} FImage;

/* Structure that stores for each face of the volume a vector that is orthogonal and a central point in this vector */
typedef struct _volumefaces {
	Matrix *orthogonal;
	Matrix *center;
} VolumeFaces;


/* Voxel and coordinate conversions */

#define GetXCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) % (s)->xsize)
#define GetYCoord(s,p) (((p) % (((s)->xsize)*((s)->ysize))) / (s)->xsize)
#define GetZCoord(s,p) ((p) / (((s)->xsize)*((s)->ysize)))
#define GetVoxelIndex(s,v) ((v.x)+(s)->tby[(v.y)]+(s)->tbz[(v.z)])

Voxel   GetVoxelCoord(Image *img, int p); 
char    ValidVoxel(Image *img, Voxel u);
int 	VoxelValue(Image *img, Voxel v);
void 	CopyVoxelSize(Image *img1, Image *img2);

/* Allocate memory to store image */ 
Image        *CreateImage(int xsize, int ysize, int zsize); 
void          DestroyImage(Image *img); /* Free memory */
FImage        *CreateFImage(int xsize, int ysize, int zsize); 
void          DestroyFImage(FImage *img); /* Free memory */
Image        *ReadImage(char *filename); /* Read image from a file */
void          WriteImage(Image *img, char *filename); /* Write image into a file */ 


/*--------------------- Voxel-based methods---------------------------- */

int           MaximumValue(Image *img); /* Maximum voxel intensity */
int           MinimumValue(Image *img); /* Minimum voxel intensity */
Image        *WindowAndLevel(Image *img, int window, int level, int H); /* Intensity transformation */
Image        *Normalize(Image *img, float minval, float maxval);

/*-------------------- Adjacency-based methods ------------------------- */

Image        *LinearFilter(Image *img, Kernel *K); /* Filter image by convolution */

/*------------------- Connectivity-based methods ------------------------*/

Image        *LabelBinaryImage(Image *img, AdjRel *A); /* Label foreground components */


/*--------------------- Task 1---------------------------- */
Image 		* GetSagitalSlice(Image *img, int slice);
Image 		* GetAxialSlice(Image *img, int slice);
Image 		* GetCoronalSlice(Image *img, int slice);
void 		WriteImageP2(Image *img, char filename[]);


/*Move functions below this line to Visualization files */
/* --------------------- Task 2---------------------------- */
int 		CompareVoxels(Voxel v1, Voxel v2);
FloatList 	* IntensityProfile(Image *img, Voxel p1, Voxel pn);
FVoxelList  * DDAAlgorithm(Voxel p1, Voxel pn);
int 		LinearInterpolationValue(Image *img, FVoxel v);
Voxel 		LinearInterpolationCoord(Image *img, FVoxel v);
void 		DrawLine(Image *img, Voxel p1, Voxel pn, int color);

/*--------------------- Task 3---------------------------- */
VolumeFaces *CreateVolumeFaces(Image *I);
void		 DestroyVolumeFaces(VolumeFaces *vf);
int 		 ComputeIntersection(Matrix *Tpo, Image *img, Matrix *Tn, VolumeFaces *cf, int *p1, int *pn);
Image 		*MaximumIntensityProjection(Image *img, float xtheta, float ytheta, float ztheta);


/*---------------------- Task 4---------------------------*/
Image 		*RayCasting(Image *img, float xtheta, float ytheta, float ztheta);
int       VolumeRenderValue(Voxel p0, Voxel p1, Voxel pn, Image *scene, Image *normalIndexImg, Matrix *normalTable, Matrix*ObserverVector, FImage *opacity);
Matrix    *CreateNormalLookUpMatrix();
FImage    *CreateOpacityImage(Image *img, Image *gradient);
void      ImageGradientMagnitudeAndIndex(Image *img, Image *gradientImg, Image *normalIndexImg, AdjRel *A);
float     PhongShading(int p, float distance, float diagonal, Matrix *N, Matrix *ObserverVector);


#endif
