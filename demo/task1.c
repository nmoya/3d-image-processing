#include "mo815-3dvis.h"
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

float MeanValue(Image *img)
{
    float number_of_pixels = 0;
    float sum = 0;
    int i =0;
    for (i=0; i< img->n; i++)
    {
        if (img->val[i] > 100)
        {
            sum += img->val[i];
            number_of_pixels++;
        }
    }
    return sum / number_of_pixels * 1.0;
}
float StdDevValue(Image *img)
{
    int i=0;
    float average = MeanValue(img);
    float accumulator = 0;
    float variance = 0;
    float number_of_pixels = 0;

    for (i=0; i<img->n; i++)
    {
        if (img->val[i] > 100)
        {
            accumulator += pow((img->val[i] - average), 2);
            number_of_pixels++;
        }
    }
    variance = (float) (accumulator*1.0) / (number_of_pixels * 1.0);
    return (float) sqrt(variance);
}

int main(int argc, char *argv[])
{
    if (argc != 3)
        Error("Run: ./main <filename> <0=axial, 1=coronal, 2=sagital>", "main");
    
    char view = atoi(argv[2]);
    int x, y, z;
    char filename[120];
    Image *img = ReadImage(argv[1]);
    Image *tempImg = NULL;

    float mean = MeanValue(img);
    float stdev = StdDevValue(img);
    printf("Mean: %.2f\n", mean);
    printf("StdDev: %.2f\n", stdev);

    Image *adjusted_img = NULL;
    adjusted_img = WindowAndLevel(img, mean, stdev*2, 4095);
    //adjusted_img = img;

    if (view == 0) //Axial
    {
        for (z=0; z<adjusted_img->zsize; z++)
        {
            tempImg = GetAxialSlice(adjusted_img, z);
            sprintf(filename, "axial%d.pgm", z);
            WriteImageP2(tempImg, filename);
            filename[0] = '\0';
            DestroyImage(tempImg);
        }
    }
    else if (view == 1) // Coronal
    {
        for (y=0; y<adjusted_img->ysize; y++)
        {
            tempImg = GetCoronalSlice(adjusted_img, y);
            sprintf(filename, "coronal%d.pgm", y);
            WriteImageP2(tempImg, filename);
            filename[0] = '\0';
            DestroyImage(tempImg);
        }
    }
    else if (view == 2) //Sagital
    {
        for (x=0; x<adjusted_img->xsize; x++)
        {
            tempImg = GetSagitalSlice(adjusted_img, x);
            sprintf(filename, "sagital%d.pgm", x);
            WriteImageP2(tempImg, filename);
            filename[0] = '\0';
            DestroyImage(tempImg);
        }
    }
    else
        Error("Wrong view option", "main");

    DestroyImage(img);
    DestroyImage(adjusted_img);

    return 0;
}