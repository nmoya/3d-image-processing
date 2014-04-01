#include "mo815-3dvis.h"
/*
Student: Nikolas Moya
Number: 144678
Email: nikolasmoya@gmail.com

Instructions:
make main       //Will generate the main binary
./main <filename> <p1.x p1.y p1.z p2.x p2.y p2.z> 

Example:
./main ../../base/small-foot.scn 0 0 0 10 10 10
./main ../../base/small-foot.scn 0 0 0 100 100 100
./main ../../base/small-foot.scn 0 0 0 0 0 190

*/


int main(int argc, char *argv[])
{
    if (argc != 8)
        Error("Run: ./main <filename> <p1.x> <p1.y> <p1.z> <p2.x> <p2.y> <p2.z> ", "main");
    
    Voxel p1, p2;
    int i;
    Image *img = ReadImage(argv[1]);
    FloatList *output = NULL;
    p1.x = atoi(argv[2]);
    p1.y = atoi(argv[3]);
    p1.z = atoi(argv[4]);

    p2.x = atoi(argv[5]);
    p2.y = atoi(argv[6]);
    p2.z = atoi(argv[7]);

    output = IntensityProfile(img, p1, p2);
    DrawLine(img, p1, p2, 4095);
    
    for(i=0; i< output->n; i++)
        printf("%g ", output->val[i]);
    printf("\n");

    WriteImage(img, "../data/line.scn");
    DestroyFloatList(output);
    DestroyImage(img);
    return 0;
}