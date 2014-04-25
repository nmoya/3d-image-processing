#include "mo815-3dvis.h"
/*
Student: Nikolas Moya
Number: 144678
Email: nikolasmoya@gmail.com

Instructions:
make mip       //Will generate the mip binary
./mip <filename> <output.pgm>

Example:
./mip ../../base/small-foot.scn mip1.pgm
./mip ../../base/small-foot.scn mip2.pgm
./mip ../../base/small-foot.scn mip3.pgm

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


    output = MaximumIntensityProfile(img, tx, ty, tz);
    sprintf(buffer, "../data/%s", argv[2]);

    WriteImage(output, buffer);
    DestroyImage(img);
    DestroyImage(output);
    return 0;
}