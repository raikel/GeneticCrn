#include "fitness.h"

void fitParamSet(Fitness *fitness, mxArray *fitnessParams)
{
	mxArray *field;
  mwSize  linksQ;
	// Check that second argument is an struct array
  if( !mxIsStruct(fitnessParams) ) {
    mexErrMsgTxt("Second input argument must be a structure array.");
  }
  // Set the the fitness function
  fitness->fitnessFcn = fitnessFcn;
  // Get number of primary links
  field = mxGetField(fitnessParams, 0, "pLinksQ");
  if( field != NULL ) {
    if ( mxIsUint32(field) && (mxGetNumberOfElements(field) == 1) ) {
      fitness->pLinksQ = ( (mwIndex*)mxGetData(field) )[0];
    } else {
      mexErrMsgTxt("Field PLINKSQ in FITNESSPARAMS struct must be an INT32 type scalar.");
    }           
  } else {
      mexErrMsgTxt("Must specified field PLINKSQ in FITNESSPARAMS struct.");
  }
  // Get number of secondary links
  field = mxGetField(fitnessParams, 0, "sLinksQ");
  if( field != NULL ) {
    if ( mxIsUint32(field) && (mxGetNumberOfElements(field) == 1) ) {
      fitness->sLinksQ = ( (mwIndex*)mxGetData(field) )[0];
    } else {
      mexErrMsgTxt("Field SLINKSQ in FITNESSPARAMS struct must be an INT32 type scalar.");
    }           
  } else {
      
  }
  // Total number of users
  linksQ = fitness->pLinksQ + fitness->sLinksQ;
  // Get channel coeficients matrix
	field = mxGetField(fitnessParams, 0, "channels");
	if( field != NULL  ) {
    if ( mxIsDouble(field) && (mxGetDimensions(field)[0] == linksQ) && (mxGetDimensions(field)[1] == linksQ) ) {
      fitness->channels = mxGetPr(field);
    } else {
       mexErrMsgTxt("Field CHANNELS must be a double square matrix.");
    } 
  } else {  
		mexErrMsgTxt("Must specified field CHANNELS in FITNESSPARAMS struct.");
	}
	// Get SINR thresold vector
	field = mxGetField(fitnessParams, 0, "sinrThr");
	if( field != NULL) {
		if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == linksQ) ) {
			fitness->sinrThr = mxGetPr(field);						
		} else {
			mexErrMsgTxt("Field SINRTHR in second input argument must be a double vector.");
		}		
	} else {
		mexErrMsgTxt("Must specified field SINRTHR in FITNESSPARAMS struct.");
	}
	// Get SINR thresold vector
	field = mxGetField(fitnessParams, 0, "sinr");
	if( field != NULL) {
		if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == linksQ) ) {
			fitness->sinr = mxGetPr(field);						
		} else {
			mexErrMsgTxt("Field SINR in second input argument must be a double vector.");
		}		
	} else {
		mexErrMsgTxt("Must specified field SINR in FITNESSPARAMS struct.");
	}
	// Get power vector
	field = mxGetField(fitnessParams, 0, "power");
	if( field != NULL) {
		if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == linksQ) ) {
			fitness->power = mxGetPr(field);						
		} else {
			mexErrMsgTxt("Field POWER in second input argument must be a double vector.");
		}		
	} else {
		mexErrMsgTxt("Must specified field POWER in FITNESSPARAMS struct.");
	}
  // Get power vector
	field = mxGetField(fitnessParams, 0, "noise");
	if( field != NULL) {
		if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == linksQ) ) {
			fitness->noise = mxGetPr(field);						
		} else {
			mexErrMsgTxt("Field NOISE in second input argument must be a double vector.");
		}		
	} else {
		mexErrMsgTxt("Must specified field NOISE in FITNESSPARAMS struct.");
	}
}

/*=================================================================
 *
 * The calling syntax is:
 *
 *		getNetSinr(sinr, channels, power, noise, linksQ)
 *
 * Inputs:
 *
 * CHANNELS is a pointer to a N-by-N square matrix where the element in colunm I, row J 
 * represents the channel cofficient from tramsmitter I to receiver J.
 * POWER is a pointer to a N-elements vector of the network power where each element J represents the
 * transmission power of transmitter J.
 * NOISE is a pointer to a N-elements vector of the network noise where each element I represents the
 * noise power at receiver I.
 * LINKSQ is the numbre of links in the system
 *
 * Returns:
 *
 * SINR is a pointer to a 1-by-N vector of the network sinr where each element I represents the
 * sinr at receiver I.
 *
 *=================================================================*/

void getNetSinr(double *sinr, double *channels, double *power, double *noise, mwSize n) {
  double signalx, interfx;
	mwIndex i, j; 
    
  for (i = 0; i < n; i++) {
   if (power[i] == 0){
     sinr[i] = 0;
   } else {
     signalx = power[i]*channels[i*n + i];
     interfx = noise[i];
     for (j = 0; j < n; j++) {
       interfx += power[j]*channels[j*n + i];
     }
     sinr[i] = signalx/(interfx - signalx);
   }
  }
}

void fitnessFcn(mwSize genomeLength, struct Gen *gen, struct Fitness *fitness)
{
	double  *sinrThr = fitness->sinrThr;
  double  *channels = fitness->channels;
  mwSize  linksQ = fitness->pLinksQ + fitness->sLinksQ;
  mwIndex i;
	double  *sinr = fitness->sinr;
  
  // Decode the chromosome
	for (i = 0; i < linksQ; i++) {		    
		fitness->power[i] = gen->chromosomes[i];
	}
  // Compute network SINR
  getNetSinr(sinr, channels, fitness->power, fitness->noise, linksQ);
  // Compute the fitness value
  gen->score[0] = 0;
  gen->score[1] = 0;
  gen->score[2] = 0;
  // The first entry of the fitness vector represents the number of primary constraints that are meet
  for (i = 0; i < fitness->pLinksQ; i++) {
    if (sinr[i] >= sinrThr[i]) (gen->score[0])++;    
	}
  // The second entry of the fitness vector represents the number of secondary constraints that are meet	
  for (i = fitness->pLinksQ; i < linksQ; i++) {
    if (sinr[i] >= sinrThr[i]) (gen->score[1])++;    
	}
  // The third entry of the fitness vector represents the total power for the secondary links	
  for (i = fitness->pLinksQ; i < linksQ; i++) {
     gen->score[2] += fitness->power[i];
	}
  gen->score[2] = -gen->score[2];
}

void decode(double *x, mwSize genomeLength, struct Gen *gen) {
  mwIndex i;
  // Decode the chromosome
	for (i = 0; i < genomeLength; i++) {		    
		x[i] = gen->chromosomes[i];
	}
}
