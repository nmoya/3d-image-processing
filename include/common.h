#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdlib.h>
#include <stdio.h>
#if !defined(__APPLE__)
	#include <malloc.h>
#endif
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include <cblas.h>
#include <dirent.h>
#include <regex.h>

/* 
 * Common data types
 */


#define INFINITY_INT  INT_MAX
#define INFINITY_FLT  FLT_MAX
#define INFINITY_DBL  DBL_MAX
#define INFINITY_LDBL LDBL_MAX

typedef struct timeval timer;

typedef unsigned char  uchar;
typedef unsigned short ushort;
typedef unsigned int   uint;
typedef unsigned long long ullong;

typedef struct _vector {
  float x,y,z;
} Vector, Point;

typedef struct _voxel {
  int x,y,z;
} Voxel;

typedef struct _dcomplex
{
  double r;
  double i;
} Complex;

/** 
 * Common definitions 
 */

#define AXIS_X 0
#define AXIS_Y 1
#define AXIS_Z 2
#define PI          3.1415926536
#define INTERIOR    0
#define EXTERIOR    1
#define BOTH        2
#define WHITE       0
#define GRAY        1
#define BLACK       2
#define NIL        -1
#define INCREASING  1
#define DECREASING  0
#define Epsilon     1E-05

/** 
 * Common operations 
 */

#ifndef MAX
#define MAX(x,y) (((x) > (y))?(x):(y))
#endif

#ifndef MIN
#define MIN(x,y) (((x) < (y))?(x):(y))
#endif

#define ROUND(x) ((x < 0)?(int)(x-0.5):(int)(x+0.5))
#define SIGN(x) ((x >= 0)?1:-1)


/** 
 * Common functions to allocate memory
 */

char        *AllocCharArray(int n);  
uchar       *AllocUCharArray(int n); 
short       *AllocShortArray(int n);
ushort      *AllocUShortArray(int n);
uint        *AllocUIntArray(int n); 
ullong      *AllocULLongArray(int n); 
int         *AllocIntArray(int n);  
float       *AllocFloatArray(int n);
double      *AllocDoubleArray(int n);
Complex     *AllocComplexArray(int n);
long double *AllocLongDoubleArray(int n);

/** 
 * Error messages 
 */

#define MSG1  "Cannot allocate memory space"
#define MSG2  "Cannot open file"

/**
 *  Error message msg is printed in function func and the program
 *  exits abnormally.
 */
                     
void Error(char *msg,char *func); 

/**
 *  Warning message msg is printed in function func and the program
 *  continues.
 */
                     
void Warning(char *msg,char *func); 

timer *Tic(); /* It marks the initial time */

timer *Toc(); /* It marks the final time */

float CompTime(timer *tic, timer *toc); /* It computes the time difference */

#endif
