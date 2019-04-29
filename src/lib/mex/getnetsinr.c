/*=================================================================
 *
 * NETSINR.C Sample .MEX file corresponding to DISTMATRIX.M
 *	        Compute distance matrix. 
 *
 * The calling syntax is:
 *
 *		sinr = distmatrix(channels, power, noise)
 *
 * CHANNELS is a square matrix where the element in colunm I, row J 
 * represents the channel cofficient from tramsmitter I to receiver J.
 * POWER is a vector of the network power where each element I represents the
 * transmission power of transmitter I.
 * NOISE is a vector of the network noise where each element I represents the
 * noise power at receiver I.
 * SINR is a vector of the network sinr where each element I represents the
 * sinr at receiver I.
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
    double *channels; 
    double *power;
    double *noise;
    double *sinr;    
    double signalx, interfx;
    mwSize n;
    mwIndex i, j;
        
    /* Check for proper number of arguments */
    
    if (nrhs != 3) { 
	   mexErrMsgTxt("Three input arguments required."); 
    } else if (nlhs > 1) {
	   mexErrMsgTxt("Too many output arguments."); 
    }

    /* Check the dimensions of input parameters */ 
    if ( mxGetN(prhs[0]) != mxGetM(prhs[0]) )  { 
        mexErrMsgTxt("First input argument must be a square matrix."); 
    }  

    /* Get number of links */
    n = mxGetN(prhs[0]);

    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
    
    /* Assign pointers to input parameters */
    channels = mxGetPr(prhs[0]); 
    power = mxGetPr(prhs[1]);
    noise = mxGetPr(prhs[2]);    

    /* Assign a pointer to the return argument */
    sinr = mxGetPr(plhs[0]);
    
    for (i = 0; i < n; i++)
    {        
        if (power[i] == 0)
        {
            sinr[i] = 0;
        }
        else 
        {
            signalx = power[i]*channels[i*n + i];
            interfx = noise[i];
            for (j = 0; j < n; j++)
            {
                interfx += power[j]*channels[j*n + i];
            }            
            sinr[i] = signalx/(interfx - signalx);
        }
    }        
    
    return;
    
}