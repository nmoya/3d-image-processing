#include "mo815-3dvis.h"

Image * CreateSphereImage(int number_of_spheres)
{
    AdjRel *A=NULL;
    int i, j, q, img_position=0;
    Voxel u, v;
    int radius=0;


    Image *img = CreateImage(100, 100, 100);
    srand(time(NULL));
    for (i=0; i< number_of_spheres; i++)
    {
        radius = rand()%20+5;
        printf("Creating Spehere: %d\n", radius);
        A = Spheric(radius);
        img_position = rand()%(img->n);
        u = GetVoxelCoord(img, img_position);
        for(j=0; j <A->n; j++)
        {
            v = GetAdjacentVoxel(A, u, j);
            if (ValidVoxel(img, v))
            {
              q = GetVoxelIndex(img,v);
              img->val[q] = 1;
            }
        }
        DestroyAdjRel(A);
    }
    return img;
}

int main(int argc, char *argv[]) 
{
  Image  *img=NULL, *labeled_image=NULL; 
  AdjRel *A=NULL;

  int number_of_spheres=5;

  if(argc != 2)
    Error("Input must be: <output.scn>", "main");

  A = Spheric(1.0);
  img = CreateSphereImage(number_of_spheres);
  labeled_image = LabelBinaryImage(img, A);
  
  WriteImage(img, "binary-image.scn");
  WriteImage(labeled_image, argv[1]);
  DestroyImage(img); 
  DestroyImage(labeled_image);


  return(0); 
}

