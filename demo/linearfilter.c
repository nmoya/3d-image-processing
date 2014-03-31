#include "mo815-3dvis.h"

int main(int argc, char *argv[]) 
{
  Image  *img[2]; 
  Kernel *K;
  timer  *t1, *t2; 

  /*--------------------------------------------------------*/

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;

  /*--------------------------------------------------------*/

  if (argc != 4)
    Error("Usage: linearfilter <input.scn> <output.scn> <adj. radius>","windowlevel.c");


  t1 = Tic();

  img[0] = ReadImage(argv[1]); 
  K      = GaussianKernel(atof(argv[3]));
  img[1] = LinearFilter(img[0],K);
  WriteImage(img[1],argv[2]); 
  DestroyImage(img[0]); 
  DestroyImage(img[1]); 
  DestroyKernel(K);

  t2 = Toc();
  fprintf(stdout,"Linear filtering in %f ms\n",CompTime(t1,t2)); 


  /* ------------------------------------------------------ */

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   

  return(0); 
}

