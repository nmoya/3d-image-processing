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

    if (argc != 5)
        Error("Usage: SurfaceRender <image.scn> <label.scn> <tilt> <spin>", "main");

    img    = ReadImage(argv[1]);
    label  = ReadImage(argv[2]);
    scene  = ImageToFImage(img);
    DestroyImage(img);

    gc     = CreateGraphicalContext(scene, label);

    SetViewDir(gc, atof(argv[3]), atof(argv[4]));
    SetObjectNormal(gc);

    SetObjectColor(gc, 1, 1.0, 1.0, 0.0);
    SetObjectVisibility(gc, 1, 1);
    SetObjectOpacity(gc, 1, 0.1);

    SetObjectColor(gc, 2, 0.0, 1.0, 1.0);
    SetObjectVisibility(gc, 2, 1);
    SetObjectOpacity(gc, 2, 1);

    SetObjectColor(gc, 3, 0.5, 0.5, 1.0);
    SetObjectVisibility(gc, 3, 1);
    SetObjectOpacity(gc, 3, 1);

    SetObjectColor(gc, 4, 0.5, 1.0, 0.5);
    SetObjectVisibility(gc, 4, 1);
    SetObjectOpacity(gc, 4, 0.1);

    SetProjectionMode(gc, RAYCASTING);

    output   = SurfaceRender(gc);


    DestroyImage(label);
    DestroyFImage(scene);
    DestroyGraphicalContext(gc);

    img = Normalize(output, 0, 255);
    WriteImageP6(img, "result.ppm");
    DestroyImage(output);
    DestroyImage(img);


    return (0);
}
