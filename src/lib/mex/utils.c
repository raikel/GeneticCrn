#include "utils.h"

double xMax(double a, double b)
{
	if(a > b) {
    return a;
  }	else {
    return b;
  }
}

double xMin(double a, double b)
{
	if(a < b) {
    return a;
  }	else {
    return b;
  }
}

// Uniform integer random number generator
int randui(int xMax, int xMin)
{
    return ( rand() % (xMax + 1) + xMin );
}

// Uniform double random number generator
double randud(double xMin, double xMax)
{
	double u = ( (double)( rand() ) )/(RAND_MAX + 1);
    return ( u*(xMax - xMin) + xMin );
}

// Normal double random number generator
double randnd(double xMean, double xStd)
{
  static double V1, V2, S;
  double X;
	
  do {
    double U1 = (double)rand() / RAND_MAX;
    double U2 = (double)rand() / RAND_MAX;

    V1 = 2 * U1 - 1;
    V2 = 2 * U2 - 1;
    S = V1 * V1 + V2 * V2;
  } while(S >= 1 || S == 0);
  
  X = V1 * sqrt(-2 * log(S) / S);
  
  return (X*xStd + xMean);
}


// Uniform double random number generator
void copyDoubleData(double *y, double *x, mwSize n)
{
  mwIndex i;
  
	for(i = 0; i < n; i++) 
  {
    y[i] = x[i];
  }
}

// Print matrix
void mxPrintMatrix(double *x, mwSize M, mwSize N) {
  mwIndex i, j;
  mexPrintf("\n*********************************\n");
  for (i =0; i<M; i++) {
    for (j =0; j<N; j++) {
      mexPrintf("%e ", x[j*M + i]);
    }
    mexPrintf("\n");
  }
  mexPrintf("*********************************\n");
}



