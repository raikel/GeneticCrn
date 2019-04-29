#include <math.h>
#include "mex.h"
#include "ga.h"
#include "fitness.h"

void mexFunction( int nlhs, mxArray *plhs[], 
				  int nrhs, const mxArray*prhs[] )
     
{
	mwIndex i,j;
  mwSize linksQ, genomeLength;
  mwSize nWorlds;
  mxArray *channels;
  Ga      ga; 
	Fitness	fitness;
  double *power;	

	// Check proper number of input arguments
	if(nrhs != 2) {
		mexErrMsgTxt("Two input arguments required.");
	}  
  // Set fitness parameters  
	fitParamSet(&fitness, prhs[1]);
  // Total number of links
  linksQ = fitness.pLinksQ + fitness.sLinksQ;
  // Set ga parameters
	gaParamSet(&ga, prhs[0]);
  // Genome lenght
  genomeLength = ga.genomeLength;
  // Numero de escenarios
  channels = mxGetField(prhs[1], 0, "channels");
  if ( mxGetNumberOfDimensions(channels) == 2 ) {
    nWorlds = 1;
  } else if ( mxGetNumberOfDimensions(channels) == 3 ) {
    nWorlds = mxGetDimensions(channels)[2];
  } else {
    mexErrMsgTxt("Parameter CHANNELS must be a bi-dimensional or tri-dimensional array.");
  }
  // Create an array for the return argument
  plhs[0] = mxCreateDoubleMatrix(linksQ, nWorlds, mxREAL);
  // Get pointer to return argument data array
  power = mxGetPr(plhs[0]);
  
 // mexPrintf("\n");
//   for(i = 0; i < linksQ; i++) {
//     for(j = 0; j < linksQ; j++) {
//       mexPrintf("%e ", fitness.channels[j*linksQ + i]);
//     }
//     mexPrintf("\n"); 
//   }
//     
  // Iterate throw the fitness worlds
  for(i = 0; i < nWorlds; i++) {
    // Copy solution to the output array
    decode(&power[i*linksQ], genomeLength, evolution(&ga, &fitness) );
    // Next channel matrix
    fitness.channels += linksQ*linksQ;
  }  
	return;
}
