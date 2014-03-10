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

    return 0;
}