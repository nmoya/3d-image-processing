/*
Student: Nikolas Moya
Number: 144678
Email: nikolasmoya@gmail.com

Instructions:
make main       //Will generate the main binary
./main <filename> <0=axial, 1=coronal, 2=sagital> 

Example:
./main ../../base/small-foot.scn 0
./main ../../base/small-foot.scn 1
./main ../../base/small-foot.scn 2

*/

#include "nmImage.h"

int main(int argc, char *argv[])
{
    if (argc != 3)
        nmError("Run: ./main <filename> <0=axial, 1=coronal, 2=sagital>", "main");
    
    char view = atoi(argv[2]);
    int x, y, z;
    char filename[120];
    nmImage *nm = nmReadSCNImage(argv[1]);
    nmImage *tempImg = NULL;

    float mean = nmMeanValue(nm);
    float stdev = nmStdDevValue(nm);
    printf("Mean: %.2f\n", mean);
    printf("StdDev: %.2f\n", stdev);

    nmImage *adjusted_img = NULL;
    adjusted_img = nmLinearStretching(nm, mean, stdev);

    if (view == 0) //Axial
    {
        for (z=0; z<nm->zsize; z++)
        {
            tempImg = nmGetAxialSlice(nm, z);
            sprintf(filename, "axial%d.pgm", z);
            nmWriteImageP2(tempImg, filename);
            filename[0] = '\0';
            nmDestroyImage(&tempImg);
        }
    }
    else if (view == 1) // Coronal
    {
        for (y=0; y<nm->ysize; y++)
        {
            tempImg = nmGetCoronalSlice(nm, y);
            sprintf(filename, "coronal%d.pgm", y);
            nmWriteImageP2(tempImg, filename);
            filename[0] = '\0';
            nmDestroyImage(&tempImg);
        }
    }
    else if (view == 2) //Sagital
    {
        for (x=0; x<nm->xsize; x++)
        {
            tempImg = nmGetSagitalSlice(nm, x);
            sprintf(filename, "sagital%d.pgm", x);
            nmWriteImageP2(tempImg, filename);
            filename[0] = '\0';
            nmDestroyImage(&tempImg);
        }
    }
    else
        nmError("Wrong view option", "main");

    nmDestroyImage(&nm);
    nmDestroyImage(&adjusted_img);

    return 0;
}