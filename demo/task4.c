#include "mo815-3dvis.h"
/*
Student: Nikolas Moya
Number: 144678
Email: nikolasmoya@gmail.com

Instructions:
make task4       //Will generate the volume rendering binary
./task4 <filename> <output.pgm> <xtheta> <ytheta> <ztheta>

Example:
./task4 ../../base/small-foot.scn volumerendering.pgm 0 0 
./task4 ../../base/small-foot.scn volumerendering.pgm 30 0 
./task4 ../../base/small-foot.scn volumerendering.pgm 45 0 

*/

int main(int argc, char *argv[])
{
    if (argc != 5)
        Error("Run: ./main <filename> <output> <tilt> <spin>", "main");
    
    char buffer[512];
    float tilt, spin;
    tilt = atof(argv[3]);
    spin = atof(argv[4]);
    Image *output = NULL, *rend = NULL;
    Image *img = ReadImage(argv[1]);
    Image *grad = NULL, *index = NULL;
    AdjRel *A = Spheric(sqrtf(3));
    Image *scene  = ImageToFImage(img);
    GraphicalContext *gc;
    
    grad = CreateImage(img->xsize, img->ysize, img->zsize);
    index = CreateImage(img->xsize, img->ysize, img->zsize);
    ImageGradientMagnitudeAndIndex(img, grad, index, A);

    gc     = CreateGraphicalContext(scene, NULL);
    SetSceneOpacity(gc, 0 /*min value*/, 255 /*max value*/, grad, 50 /*gradient threshold*/, 1 /*opacity*/);
    SetViewDir(gc, tilt, spin);
    
    rend = VolumeRender(gc);

    output = Normalize(rend, 0, 255);
    sprintf(buffer, "%.0f%0.f%s", tilt, spin, argv[2]);
    WriteImageP2(output, buffer);

    DestroyImage(img);
    DestroyImage(output);
    DestroyImage(rend);
    DestroyFImage(scene);
    return 0;
}