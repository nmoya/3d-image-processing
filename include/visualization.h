#ifndef _VISUALIZATION_H
#define _VISUALIZATION_H

#include "mo815-3dvis.h"

#define  NUM_OF_NORMALS    65161
#define  SCENE_NORMAL      0  
#define  OBJECT_NORMAL     1
#define  RAYCASTING        0
#define  SPLATTING         1

/* Graphical context with the visualization attributes for ray casting
   and voxel splatting, including current viewing direction
   parameters.
     
       Octants of the Scene
 z
  /
 /
 ---- x
 |
y|
       4 ------- 5 
      /        /|
     /        / |
    /        /  | 
  0 --------1   |  
   |  6     |  / 7  
   |        | / 
   |        |/  
  2 -------- 3 


The origin of the image to be rendered on the viewing plane is
initially at (0,0,-diag/2). The observer is initially at (diag/2,
diag/2, -infinity) with a single white light source at the same
position. The transformation T is applied to the scene for voxel
splatting or the inverse of T is applied to the vieweing plane for ray
casting. The front-to-back (FTB) voxel access order for splatting
depends on the pair (closest octant, principal axis). The FTB
possibilities (dx,dy,dz,xo,yo,zo,xf,yf,zf) are stored in a 2D array
whose column is the closest octant and line is the principal axis. Any
change in viewing direction requires to update all viewing parameters below.

*/

/* Extra surface rendering buffers, especially required for voxel
   splatting. They all have the same dimensions of the resulting
   rendered image (diag x diag pixels) to avoid clipping. */

typedef struct _plane {
  Point  pos;    // reference position
  Vector normal; // normal vector which gives its orientation
} Plane;

typedef struct _surface_rendering_buffers {
  float  depth;         /* distance to the viewing plane, whose center
			   is initially positioned at
			   (diag/2,diag/2,-diag/2) */ 
  float  opacity;       /* accumulated opacity of the voxels projected
			   onto the viewing plane */
  int    voxel;         /* voxels in the scene that have been projected
			   onto the viewing plane */
  int    object;        /* object in the scene whose voxel has been
			   projected onto the viewing plane */
} SRBuffers;
  
/* Object opacity, visibility, and color attributes in the scene */

typedef struct _object_attributs {
  float opacity;          /* opacity value in [0,1] */
  float red, green, blue; /* proportions in [0,1] of red, green and
                             blue components in its color */
  char visibility;        /* visibility flag (0/1 for invisible/visible) */
} ObjectAttributes; 

/* Parameters of the Phong's illumination model */

typedef struct _phong_model {
  float      ka;     /* ambient reflection constant  (0.1 by default) */
  float      kd;     /* diffuse reflection constant  (0.7 by default) */
  float      ks;     /* specular refrection constant (0.2 by default) */
  float      ns;     /* exponent of the specular reflection component (5 by default) */
  Vector *normal; /* normal table used to speed up shading */ 
  float     *Idist;  /* distance table used to convert voxel's distance to the viewing plane into intensity in [0,1] (also known as depth shading) */ 
  int        ndists; /* number of distance values, which is the size of the scene's diagonal. */
} PhongModel; 

/* Viewing direction and FTB voxel access order */

typedef struct _ftb_access {
  int dx, dy, dz; /* -1 or 1. Default is 1,1,1 */
  int xo, yo, zo; /* the x,y,z original coordinates of the closest  octant: the initial FTB voxel access coordinates */
  int xf, yf, zf; /* the x,y,z original coordinates of the farthest octant: the final FTB voxel access coordinates */
} FTBaccess;

typedef struct _viewing_direction {
  Matrix           *T;            /* scene's transformation matrix: default is spin=tilt=0.0 */ 
  Matrix           *Tinv;         /* inverse of the scene's transformation matrix */ 
  Matrix           *R;            /* rotation matrix, only for vectors */
  Matrix           *Rinv;         /* inverse of the rotation matrix, only for vectors */
  char                 paxis;        /* principal axis, 0 for 'x', 1 for 'y', or 2 for 'z', which is the last one to be visited in the voxel splatting loop. This is the most orthogonal axis to the viewing plane. */
  int                  octant;       /* closest octant, 0, 1, ..., 7, to the vieweing plane, which together with the principal axis indicate the front-to-back voxel access order for splatting. */
  FTBaccess        *ftb;          /* gives the front-to-back access order from the closest octant to the farthest octant. Used for voxel splatting only. */ 
  Vector            V;            /* visualization direction --- vector that points to the observer. */
} ViewDir; 


/* Graphical context of the scene for 3D visualization */

typedef struct _graphic_context {
  ObjectAttributes *object;       /* list of attributes per object
					in the scene */
  int                  nobjs;        /* number of objects in the scene */  
  float                overall_opac; /* it indicates the overall opacity of the objects in the scene */ 
  char                 proj_mode;    /* voxel projection mode: ray casting is 0 and voxel splatting is 1 */ 
  PhongModel       *phong;        /* rendering parameters of the Phong's illumination model */
  ViewDir          *viewdir;      /* viewing direction parameters */
  SRBuffers        *surf_render;  /* extra surface rendering buffers, especially required for voxel splatting */
  FImage           *scene;        /* the 3D image property for visualization (intensity, BMD, etc.) */
  Image            *label;        /* the 3D label image with objects for visualization, used only for surface rendering */
  Plane            *face;         /* planes of the 6 faces of the scene  */
  Image            *normal;       /* normal indices of the scene's voxels */
  FImage           *opacity;      /* opacity scene used for volume rendering only */
} GraphicalContext;


/* 
   Graphical Context for 3D visualization of the image properties in
   scene. The object labels i > 0 (0 is background) are used only for
   surface rendering (i.e., label==NULL for volume rendering). The
   normal vectors may be estimated as the gradient of the scene
   (normal_type=SCENE_NORMAL), as the gradient of the EDT of the
   objects (normal_type=OBJECT_NORMAL) or on-the-fly from the index
   buffer (normal_type=NIL).

*/

void                 SetSceneNormal(GraphicalContext *gc);
void                 SetProjectionMode(GraphicalContext *gc, char proj_mode);
void                 SetObjectNormal(GraphicalContext *gc);
void                 SetObjectColor(GraphicalContext *gc, int object, float red, float green, float blue);
void                 SetSceneOpacity(GraphicalContext *gc, float min_val, float max_val, Image *grad, int grad_thres, float max_opac);
void                 SetObjectOpacity(GraphicalContext *gc, int object, float opacity);
void                 SetObjectVisibility(GraphicalContext *gc, int object, char visibility);
void                 SetViewDir(GraphicalContext *gc, float tilt, float spin);
char                 IntersectionPoints(GraphicalContext *gc, Point P0, Vector n, Point *P1, Point *Pn);
Image            *SurfaceRender(GraphicalContext *gc);
Image            *VolumeRender(GraphicalContext *gc);
GraphicalContext *CreateGraphicalContext(FImage *scene, Image *label); 
void                 DestroyGraphicalContext(GraphicalContext *gc);

// void      DrawPoint(Image *img, Voxel u, Color YCbCr, AdjRel *B);
// void      DrawPoints(Image *img, Set *S, Color YCbCr, AdjRel *B);
// void      DrawLine(Image *img, Voxel u1, Voxel u2, Color YCbCr, AdjRel *B);
// void      DrawBorders(Image *img, Image *label, AdjRel *A, Color YCbCr, AdjRel *B);
// void      DrawObject(Image *img, Image *bin, Color YCbCr, AdjRel *B);
// Image *Draw2DFeatureSpace(DataSet *Z, uchar opt, uchar status);
// Image *DrawVoxelSamples(DataSet *Z, uchar opt, char *ref_data_type);
// Image *ColorizeComp(Image *label);
// Image *ColorizeCompOverImage(Image *orig, Image *label);
// Image  *ProjectMaxValue(Image *img,  Plane *cutplane, int xviewsize, int yviewsize, int slabsize); 
// Image  *ProjectMinValue(Image *img,  Plane *cutplane, int xviewsize, int yviewsize, int slabsize);
// Image  *ProjectMeanValue(Image *img, Plane *cutplane, int xviewsize, int yviewsize, int slabsize);
// Image  *ProjectMaxObjectValue(Image *img,  Image *obj, Plane *cutplane, int xviewsize, int yviewsize, int slabsize); 
// Image  *ProjectMinObjectValue(Image *img,  Image *obj, Plane *cutplane, int xviewsize, int yviewsize, int slabsize);
// Image  *ProjectMeanObjectValue(Image *img, Image *obj, Plane *cutplane, int xviewsize, int yviewsize, int slabsize);
// Image  *CurvilinearSplatting(Image *img, FImage *dist, Plane *cutplane, int viewsize, float depth);


#endif

