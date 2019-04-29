#ifndef utils_h
#define utils_h

#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <string.h>

double xMin(double, double);
double xMax(double, double);
int randui(int, int);
double randud(double, double);
double randnd(double, double);
void copyDoubleData(double *, double *, mwSize);
void mxPrintMatrix(double *, mwSize , mwSize );

#endif
