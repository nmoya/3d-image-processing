#include "mo815-3dvis.h"
/*
Student: Nikolas Moya
Number: 144678
Email: nikolasmoya@gmail.com

Instructions:
make task3       //Will generate the mip binary
./task3 <filename> <output.pgm> <xtheta> <ytheta> <ztheta>

Example:
./task3 ../../base/small-foot.scn mip.pgm 0 0 0
./task3 ../../base/small-foot.scn mip.pgm 30 0 0
./task3 ../../base/small-foot.scn mip.pgm 45 0 0

*/

int main(int argc, char *argv[])
{
    if (argc != 6)
        Error("Run: ./main <filename> <output> <xtheta> <ytheta> <ztheta>", "main");
    
    char buffer[512];

    float tx, ty, tz;
    tx = atof(argv[3]);
    ty = atof(argv[4]);
    tz = atof(argv[5]);
    Image *img = ReadImage(argv[1]);
    Image *output = NULL;

    output = RayCasting(img, tx, ty, tz);
    sprintf(buffer, "../data/%.1f%.1f%.1f%s", tx, ty, tz, argv[2]);
   
    WriteImageP2(output, buffer);
    DestroyImage(img);
    DestroyImage(output);
    return 0;
}