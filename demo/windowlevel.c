#include "mo815-3dvis.h"

void EstimateWindowAndLevel(Image *img, int *window, int *level, int minimum)
{
  int p, n; 
  float stdev;

  *level = 0; n = 0;  
  for (p=0; p < img->n; p++) 
    if (img->val[p] > minimum){ 
      *level += img->val[p];
      n++;
    }
  *level /= n; 
  
  stdev = 0.0; 
  for (p=0; p < img->n; p++) 
    if (img->val[p] > minimum){ 
      stdev += (img->val[p] - *level)*(img->val[p] - *level); 
    }
  stdev = sqrtf(stdev/n); 

  *window = (int) 2*stdev; 
}

int main(int argc, char *argv[]) 
{
  Image *img[2]; 
  int    window, level;
  timer *t1, *t2; 

  /*--------------------------------------------------------*/

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;

  /*--------------------------------------------------------*/

  if (argc != 3)
    Error("Usage: windowlevel <input.scn> <output.scn>","windowlevel.c");


  t1 = Tic();

  img[0] = ReadImage(argv[1]); 
  EstimateWindowAndLevel(img[0],&window,&level,100); 
  printf("window %d level %d\n",window,level);
  img[1] = WindowAndLevel(img[0],window,level,4095); 
  WriteImage(img[1],argv[2]); 
  DestroyImage(img[0]); 
  DestroyImage(img[1]); 

  t2 = Toc();
  fprintf(stdout,"Intensity correction in %f ms\n",CompTime(t1,t2)); 


  /* ------------------------------------------------------ */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   

  return(0); 
}

