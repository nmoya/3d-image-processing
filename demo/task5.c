#include "mo815-3dvis.h"
/*
Student: Nikolas Moya
Number: 144678
Email: nikolasmoya@gmail.com

Instructions:
make task5       //Will generate the mip binary
./task5 <filename> <label> <output.pgm> <tilt> <spin>

Example:
./task5 ../../base/small-foot.scn surface.pgm 0 0

*/

int main(int argc, char *argv[])
{
    Image  *img, *label, *output;
    FImage *scene;
    GraphicalContext *gc;
    timer     *t1 = NULL, *t2 = NULL;
    char buffer[120];
    float tilt, spin;

    if (argc != 5)
        Error("Usage: SurfaceRender <image.scn> <label.scn> <tilt> <spin>", "main");

    img    = ReadImage(argv[1]);
    label  = ReadImage(argv[2]);
    scene  = ImageToFImage(img);
    DestroyImage(img);

    gc     = CreateGraphicalContext(scene, label);

    tilt = atof(argv[3]);
    spin = atof(argv[4]);
    SetViewDir(gc, tilt, spin);

    //SetObjectNormal(gc);
    SetObjectColor(gc, 1, 1.0, 1.0, 0.0);
    SetObjectVisibility(gc, 1, 1);
    SetObjectOpacity(gc, 1, 1);

    SetProjectionMode(gc, RAYCASTING);
    
    output   = SurfaceRender(gc);

    DestroyImage(label);
    DestroyFImage(scene);
    DestroyGraphicalContext(gc);

    img = Normalize(output, 0, 255);

    sprintf(buffer, "surfacer%.0f%.0f.ppm", tilt, spin);
    WriteImageP6(img, buffer);
    DestroyImage(output);
    DestroyImage(img);


    return (0);
}
