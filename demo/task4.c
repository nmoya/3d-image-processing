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
    //char buffer2[512];

    float tx, ty, tz;
    tx = atof(argv[3]);
    ty = atof(argv[4]);
    tz = 0;
    Image *img = ReadImage(argv[1]);
    Image *output = NULL;

    output = RayCasting(img, tx, ty, tz);
    sprintf(buffer, "%.1f%.1f%.1f%s", tx, ty, tz, argv[2]);
    //sprintf(buffer2, "convert %s %s.png", buffer, buffer);
    //system(buffer2);

    WriteImageP2(output, buffer);
    DestroyImage(img);
    DestroyImage(output);
    return 0;
}