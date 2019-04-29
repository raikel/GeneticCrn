/*=================================================================
 *
 * DISTMATRIX.C	Sample .MEX file corresponding to DISTMATRIX.M
 *	        Compute distance matrix. 
 *
 * The calling syntax is:
 *
 *		d = distmatrix(x1, x2)
 * In D each colum J contains the distances form point J in X1 to every 
 * points in X2. In X1 and X2 each colum represents a point.
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2006 The MathWorks, Inc.
 *
 *=================================================================*/
/* $Revision: 1.10.6.4 $ */
#include <math.h>
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *x1; 
    double *x2;
    double *y;
    double d;
    mwIndex i, j, k;
    mwSize nPoints, nDims; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 2) { 
	mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    }

    /* Check the dimensions of input parameters */ 
    if ( ( mxGetN(prhs[0]) != mxGetN(prhs[1]) ) ||
         ( mxGetM(prhs[0]) != mxGetM(prhs[1]) ) ) { 
        mexErrMsgTxt("DISTMATRIX requires that input parameters have the same size."); 
    }  
    
    nPoints = mxGetN(prhs[0]);
    nDims   = mxGetM(prhs[0]);

    /* Assign pointers to input parameters */ 
    x1 = mxGetPr(prhs[0]); 
    x2 = mxGetPr(prhs[1]);
 
    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(nPoints, nPoints, mxREAL); 
    
    y = mxGetPr(plhs[0]);
    
    for (i = 0; i < nPoints; i++)
    {
        for (j = 0; j < nPoints; j++)
        {
            d = 0;
            for (k = 0; k < nDims; k++)
            {
                d = d + pow(x1[i*nDims + k] -  x2[j*nDims + k], 2);
            }
            d = sqrt(d);
            y[i*nPoints + j] = d;
        }
    }        
    
    return;
    
}