/*=================================================================
 *
 * DISTMATRIX.C	Sample .MEX file corresponding to DISTMATRIX.M
 *	        Compute distance matrix. 
 *
 * The calling syntax is:
 *
 *		d = distmatrix(x1, x2)
 * In X1 and X2 each column represents a point.
 * In D each column J contains distances form point J in X1 to every points in X2
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2006 The MathWorks, Inc.
 *
 *=================================================================*/
/* $Revision: 1.10.6.4 $ */
#include "ga.h"

void gaParamSet(Ga *ga, mxArray *gaOptions)
{
    
    char    string[STRLEN];
    mxArray *field;
    mwIndex i;
    double xOverFraction;

    ga->genomeLength = 0; 
    ga->popSize = 20;
    ga->nEliteKids = 2;
    xOverFraction = 0.8;    

    
    // Set default fitness vector lenght
    ga->nScores = 1;
    
    // Set default stoping criteria
    (ga->stopping).nGenerations = 50;

    // Set default scaling schema parameters
    (ga->fitScaling).scalingFcn = fitScalingRank;
    (ga->fitScaling).quantity = 0.4;
    (ga->fitScaling).maximumSurvivalRate = 2;

    // Set default selection schema parameters
    (ga->selection).selectionFcn = selectionRoulette;
    (ga->selection).nTournamentPlayers = 2;

    // Set crossover schema parameters
    (ga->crossover).crossoverFcn = crossoverScattered;
    
    // Set mutation schema parameters
    (ga->mutation).mutationFcn = mutationGaussian;
    (ga->mutation).rate = 0.05;    
    (ga->mutation).shrink = 1;
    (ga->mutation).scale = 1;    
    (ga->mutation).span = NULL;

    ga->classMask = NULL;
    ga->lowerBound = NULL;
    ga->upperBound = NULL;

    // Check that second argument is an struct array
    if( !mxIsStruct(gaOptions) )
    {
        mexErrMsgTxt("GAOPTIONS must be a structure array.");
    }

    // Get genomeLength
    field = mxGetField(gaOptions, 0, "genomeLength");
    if( field != NULL ) {
      if ( mxIsUint32(field) && (mxGetNumberOfElements(field) == 1) )
      {
        ga->genomeLength = ( (mwIndex*)mxGetData(field) )[0];
      } else {
        mexErrMsgTxt("Field GENOMELENGTH in GAOPTIONS struct must be an INT32 type scalar.");
      }           
    } else {
        mexErrMsgTxt("Must specified field genomeLength in GAOPTIONS struct.");
    }
    // Get nScores
    field = mxGetField(gaOptions, 0, "nScores");
    if( field != NULL ) {
      if ( mxIsUint32(field) && (mxGetNumberOfElements(field) == 1) )
      {
        ga->nScores = ( (mwIndex*)mxGetData(field) )[0];
      } else {
        mexErrMsgTxt("Field NELITEKIDS in GAOPTIONS struct must be INT32 type scalar.");
      }
    }
    // Get nEliteKids
    field = mxGetField(gaOptions, 0, "nEliteKids");
    if( field != NULL )
    {
        if ( mxIsUint32(field) && (mxGetNumberOfElements(field) == 1) )
        {
            ga->nEliteKids = ( (mwIndex*)mxGetData(field) )[0];   
        }
        else
        {
            mexErrMsgTxt("Field NELITEKIDS in GAOPTIONS struct must be INT32 type scalar.");
        }           
    }
    // Get popSize
    field = mxGetField(gaOptions, 0, "popSize");
    if( field != NULL )
    {
        if ( mxIsUint32(field) && (mxGetNumberOfElements(field) == 1)  )
        {
            ga->popSize = ( (mwIndex*)mxGetData(field) )[0]; 
        }
        else
        {
            mexErrMsgTxt("Field POPSIZE in GAOPTIONS struct must be INT32 type scalar.");
        }           
    }
    // Get xOverFraction
    field = mxGetField(gaOptions, 0, "xOverFraction");
    if( field != NULL )
    {
        if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == 1) )
        {
            xOverFraction = mxGetScalar(field); 
        }
        else
        {
            mexErrMsgTxt("Field XOVERFRACTION in GAOPTIONS struct must be INT16 type.");
        }           
    }
    // Get scalingQuantity
    field = mxGetField(gaOptions, 0, "scalingQuantity");
    if( field != NULL )
    {
        if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == 1) )
        {
            (ga->fitScaling).quantity = mxGetScalar(field);  
        }
        else
        {
            mexErrMsgTxt("Field SCALINGQUANTITY in GAOPTIONS struct must be DOUBLE type.");
        }           
    }
    // Get scalingMaximumSurvivalRate
    field = mxGetField(gaOptions, 0, "scalingMaximumSurvivalRate");
    if( field != NULL )
    {
        if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == 1) )
        {
            (ga->fitScaling).maximumSurvivalRate = mxGetScalar(field);
        }
        else
        {
            mexErrMsgTxt("Field SCALINGMAXIMUMSURVIVALRATE in GAOPTIONS struct must be DOUBLE type.");
        }           
    }
    // Get nTournamentPlayers 
    field = mxGetField(gaOptions, 0, "nTournamentPlayers");
    if( field != NULL )
    {
        if ( mxIsUint32(field) && (mxGetNumberOfElements(field) == 1) )
        {
            (ga->selection).nTournamentPlayers = ( (mwIndex*)mxGetData(field) )[0];  
        }
        else
        {
            mexErrMsgTxt("Field NTOURNAMENTPLAYERS in GAOPTIONS struct must be INT16 type.");
        }           
    }
    // Get mutationRate
    field = mxGetField(gaOptions, 0, "mutationRate");
    if( field != NULL )
    {
        if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == 1) )
        {
            (ga->mutation).rate = mxGetScalar(field); 
        }
        else
        {
            mexErrMsgTxt("Field SCALINGMAXIMUMSURVIVALRATE in GAOPTIONS struct must be DOUBLE type.");
        }           
    }
    // Get mutationShrink
    field = mxGetField(gaOptions, 0, "mutationShrink");
    if( field != NULL )
    {
        if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == 1) )
        {
            (ga->mutation).shrink = mxGetScalar(field);   
        }
        else
        {
            mexErrMsgTxt("Field MUTATIONSHRINK in GAOPTIONS struct must be DOUBLE type.");
        }           
    }
    // Get mutationScale
    field = mxGetField(gaOptions, 0, "mutationScale");
    if( field != NULL )
    {
        if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == 1) )
        {
            (ga->mutation).scale = mxGetScalar(field); 
        }
        else
        {
            mexErrMsgTxt("Field MUTATIONSCALE in GAOPTIONS struct must be DOUBLE type.");
        }           
    }
    // Get nGenerations 
    field = mxGetField(gaOptions, 0, "nGenerations");
    if( field != NULL )
    {
        if ( mxIsUint32(field) && (mxGetNumberOfElements(field) == 1) )
        {
           (ga->stopping).nGenerations = ( (mwIndex*)mxGetData(field) )[0];  
        }
        else
        {
            mexErrMsgTxt("Field NGENERATIONS in GAOPTIONS struct must be INT16 type.");
        }           
    }
    // Get lowerBound
    field = mxGetField(gaOptions, 0, "lowerBound");
    if( field != NULL )
    {
        if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == ga->genomeLength) )
        {
          ga->lowerBound = mxGetPr(field);                     
        }
        else
        {
            mexErrMsgTxt("Field LOWERBOUND in GAOPTIONS struct must be DOUBLE type.");
        }           
    }
    // Get upperBound
    field = mxGetField(gaOptions, 0, "upperBound");
    if( field != NULL )
    {
        if ( mxIsDouble(field) && (mxGetNumberOfElements(field) == ga->genomeLength) )
        {
          ga->upperBound = mxGetPr(field); 
        }
        else
        {
            mexErrMsgTxt("Field UPPERBOUND in GAOPTIONS struct must be DOUBLE type.");
        }           
    }
    // Get classMask
    field = mxGetField(gaOptions, 0, "classMask");
    if( field != NULL )
    {
        if ( mxIsUint32(field) && (mxGetNumberOfElements(field) == ga->genomeLength) )
        {
          ga->classMask = (mwIndex*)mxGetData(field); 
        }
        else
        {
            mexErrMsgTxt("Field CLASSMASK in GAOPTIONS struct must be INT16 type.");
        }           
    }
    // Get fitScalingFcn
    field = mxGetField(gaOptions, 0, "fitScalingFcn");
    if( field != NULL )
    {
        if ( mxIsChar(field) )
        {
            mxGetString(field, string, STRLEN);
            if      ( !strcmp(string, "fitScalingRank") )
            {
                (ga->fitScaling).scalingFcn = fitScalingRank;
            }
            else if ( !strcmp(string, "fitScalingTop" ) )
            {
                (ga->fitScaling).scalingFcn = fitScalingTop;
            }
            else if ( !strcmp(string, "fitScalingShiftLinear") )
            {
                //(ga->fitScaling).scalingFcn = fitScalingShiftLinear;
            }
            else
            {
                mexErrMsgTxt("Unknown fitness scaling funcion.");
            }                        
        }
        else
        {
            mexErrMsgTxt("Field FITSCALINGFCN in GAOPTIONS struct must be CHAR type.");
        }           
    }
    // Get selectionFcn
    field = mxGetField(gaOptions, 0, "selectionFcn");
    if( field != NULL )
    {
        if ( mxIsChar(field) )
        {
            mxGetString(field, string, STRLEN);
            if      ( !strcmp(string, "selectionRoulette" ) )
            {
                (ga->selection).selectionFcn = selectionRoulette;
            }
            else if ( !strcmp(string, "selectionTournament") )
            {
                (ga->selection).selectionFcn = selectionTournament;
            }
            else
            {
                mexErrMsgTxt("Unknown selection funcion.");
            }                
        }
        else
        {
            mexErrMsgTxt("Field SELECTIONFCN in GAOPTIONS struct must be CHAR type.");
        }           
    }
    // Get crossoverFcn
    field = mxGetField(gaOptions, 0, "crossoverFcn");
    if( field != NULL )
    {
        if ( mxIsChar(field) )
        {
            mxGetString(field, string, STRLEN);
            if      ( !strcmp(string, "crossoverScattered") )
            {
                (ga->crossover).crossoverFcn = crossoverScattered;
            }
            else
            {
                mexErrMsgTxt("Unknown crossover funcion.");
            }                
        }
        else
        {
            mexErrMsgTxt("Field CROSSOVERFCN in GAOPTIONS struct must be CHAR type.");
        }           
    }
    // Get mutationFcn
    field = mxGetField(gaOptions, 0, "mutationFcn");
    if( field != NULL )
    {
        if ( mxIsChar(field) )
        {
            mxGetString(field, string, STRLEN);            
            if      ( !strcmp(string, "mutationGaussian" ) )
            {
                (ga->mutation).mutationFcn = mutationGaussian;
            }
            else
            {
                mexErrMsgTxt("Unknown mutation funcion.");
            }                
        }
        else
        {
            mexErrMsgTxt("Field MUTATIONFCN in GAOPTIONS struct must be CHAR type.");
        }           
    }
    // Set ga size
    ga->nXoverKids = xOverFraction*(ga->popSize - ga->nEliteKids);
    ga->nMutateKids = ga->popSize - (ga->nEliteKids + ga->nXoverKids);
    //
    (ga->mutation).generations = (ga->stopping).nGenerations;
    // Allocate memory for the ga
    ga->thisPopulation = mxCalloc( ga->popSize, sizeof(Gen) );
    ga->nextPopulation = mxCalloc( ga->popSize, sizeof(Gen) );
    // Allocate memory for genes in each chromosome
    for (i = 0; i < ga->popSize; i++)
    {
        (ga->thisPopulation[i]).chromosomes = mxCalloc( ga->genomeLength, sizeof(double) );
        (ga->thisPopulation[i]).score = mxCalloc( ga->nScores, sizeof(double) );
        (ga->nextPopulation[i]).chromosomes = mxCalloc( ga->genomeLength, sizeof(double) );
        (ga->nextPopulation[i]).score = mxCalloc( ga->nScores, sizeof(double) );
    }
    // Set the default mask class
    if(ga->classMask == NULL) 
    {
       ga->classMask = mxCalloc( ga->genomeLength, sizeof(mwIndex) );
        for (i = 0; i < ga->genomeLength; i++)
        {
            ga->classMask[i] = CLASS_DOUBLE;
        } 
    }

    // Set the default upper and lower bounds
    if( (ga->lowerBound == NULL) || (ga->upperBound == NULL) ) 
    {
        ga->lowerBound = mxCalloc( ga->genomeLength, sizeof(double) );
        ga->upperBound = mxCalloc( ga->genomeLength, sizeof(double) );
        for (i = 0; i < ga->genomeLength; i++)
        {
            ga->lowerBound[i] = 0;
            ga->upperBound[i] = 1;
        } 
    }

    // Set the default mutation span
    if( (ga->mutation).span == NULL ) 
    {
        (ga->mutation).span = mxCalloc( ga->genomeLength, sizeof(double) );
        for (i = 0; i < ga->genomeLength; i++)
        {
           (ga->mutation).span[i] = ga->upperBound[i] - ga->lowerBound[i];            
        } 
    }
}

void moveGen(Gen *genX, Gen *genY, mwSize genomeLength)
{
    mwIndex i;
    for (i = 0; i < genomeLength; i++)
    {
        genX->chromosomes[i] = genY->chromosomes[i];
    }
}

void printPopulation(Ga *ga)
{
	mwIndex i, j;  
  mexPrintf("\n ************** Next Genome ************** \n");
	for (i = 0; i < ga->popSize; i++) {     
      for (j = 0; j < ga->genomeLength; j++) {
          mexPrintf("%f ", (ga->nextPopulation[i]).chromosomes[j]);
      }
      mexPrintf("   |   [");
      for (j = 0; j < ga->nScores; j++) {
          mexPrintf("%e ", (ga->nextPopulation[i]).score[j]);
      }
      mexPrintf("]   |   %e \n", (ga->nextPopulation[i]).expectation);      
  }
  
  mexPrintf("\n ************** This Genome ************** \n");
	for (i = 0; i < ga->popSize; i++) {     
      for (j = 0; j < ga->genomeLength; j++) {
          mexPrintf("%f ", (ga->thisPopulation[i]).chromosomes[j]);
      }
      mexPrintf("   |   [");
      for (j = 0; j < ga->nScores; j++) {
          mexPrintf("%e ", (ga->thisPopulation[i]).score[j]);
      }
      mexPrintf("]   |   %e \n", (ga->thisPopulation[i]).expectation);
  }
}

void initPopulation(Ga *ga)
{
	mwIndex i, j;
	
	for (i = 0; i < ga->popSize; i++)
    {      
        for (j = 0; j < ga->genomeLength; j++)
        {
            switch(ga->classMask[j])
            {
                case CLASS_DOUBLE:                
                    (ga->nextPopulation[i]).chromosomes[j] = 
                    randud(ga->lowerBound[j], ga->upperBound[j]);
                    break;
                
                case CLASS_BITSTR:                
                    if ( randud(0, 1) > 0.5 )
                    {
                        (ga->nextPopulation[i]).chromosomes[j] = 1;              
                    }
                    else
                    {
                        (ga->nextPopulation[i]).chromosomes[j] = 0;
                    } 
                    break;
                default:
                    mexErrMsgTxt("Unkwon gen class.");
            }
        }		
    }
}

void boundGen(Gen *gen, double *lowerBound, double *upperBound, mwSize genomeLength)
{
    mwIndex i;    
    for (i = 0; i < genomeLength; i++)
    {
       gen->chromosomes[i] = xMin(upperBound[i], gen->chromosomes[i]);       
       gen->chromosomes[i] = xMax(lowerBound[i], gen->chromosomes[i]);
    }
}

/* Prototipo de la funcion de comparar */
int compareGen(const void *genX, const void *genY)
{
  int i;
  for(i = 0; i < 3; i++) {
    if      ( ( (Gen*)genX )->score[i] >  ( (Gen*)genY )->score[i] ) return -1;
    else if ( ( (Gen*)genX )->score[i] <  ( (Gen*)genY )->score[i] ) return  1;
    else if ( ( (Gen*)genX )->score[i] == ( (Gen*)genY )->score[i] ) continue;    
  }
  return  0;
}

//  Fitness scaling is the process of mapping raw fitness scores
//  as returned by the fitness function, to an expected number of 
//  children for each individual. It is important that the expected 
//  number of children for each individual be in a "good" range. 
//       
//  If the expectations have too large a range, then the individuals with
//  the highest expectation will reproduce too rapidly, taking over the 
//  ga  gene pool too quickly, and preventing the Genetic Algorithm 
//  from searching other areas of the solution space.   On the other hand, if 
//  the expectations do not vary enough, then all individuals have about the 
//  same chance of reproduction and the search will progress very slowly. 
//       
//  The fitness scaling function attempts to map an arbitrary fitness range 
//  into a good expectation range. Experimentally, it has been found that if 
//  the most fit individual reproduces about twice as often as the average 
//  individual the Genetic Algorithm progresses about as quickly as possible.
//  
//  The scaling function expects that genes in ga be ordered in 
//  asscendig order according its scores. The function computes the expected 
//  fraction of offspring for each member of the ga over the total 
//  number of parents. The sum of the expectations should be equal to 1.

/* Rank based fitness scaling.
*   This function calculates the fitness of each gen in 
*   POPULATION using its scores. This relationship can be linear or 
*   nonlinear. The genes in POPULATION must be oredered 
*   in asscendig order according its scores */
void fitScalingRank(Gen *population, FitScaling *fitScaling, mwSize popSize)
{
    double expectation, expectationTotal;
    mwIndex i;    
    
    expectationTotal = 0;
    for (i = 0; i < popSize; i++)
    {
      expectation = 1/sqrt(i + 1);
      (population[i]).expectation = expectation;
      expectationTotal += expectation;
    }

    for (i = 0; i < popSize; i++)
    {
        (population[i]).expectation = (population[i]).expectation/expectationTotal;
    }
}

/* Top individuals reproduce equally.
*   This function calculates the fitness of each gen in 
*   POPULATION using its scores and QUANTITY in
*   SELECTION. QUANTITY represents 
*   the number of expectations. The genes in 
*   POPULATION must be oredered in asscendig order 
*   according its scores */
void fitScalingTop(Gen *population, FitScaling *fitScaling, mwSize popSize)
{
    mwIndex i;
    double nTopGenes = floor( popSize*(fitScaling->quantity) );
    
    if ( (nTopGenes > popSize) || (nTopGenes < 0) ) 
    {
        mexErrMsgTxt("QUANTITY must be between 0 and 1.");
    }

    for (i = 0; i < nTopGenes; i++)
    {
        (population[i]).expectation = 1/nTopGenes;
    }

    for (i = nTopGenes; i < popSize; i++)
    {
        (population[i]).expectation = 0.0;
    }
}

/* Offset and scale fitness to desired range.
*   This function calculates FITNESS of each gen in 
*   struct POPULATION using its SCORE, and NPARENTS and 
*   MAXIMUMSURVIVALRATE in struct FITSCALING. 
*   MAXIMUMSURVIVALRATE is the ratio of the expectation of the best
*   individual to the average expectation of the population. 
*   Values near 2 have been found to work well. The genes in 
*   POPULATION must be oredered in asscendig order according its SCORE 
void fitScalingShiftLinear(Gen *population, FitScaling *fitScaling, mwSize popSize)
{
    double minScore, maxScore, meanScore;
    double desiredMean, scale, offset;    
    mwIndex i;

    // Compute mean score
    meanScore = 0;
    for (i = 0; i < popSize; i++)
    {
        meanScore = meanScore + (population[i]).score;
    }
    meanScore = meanScore/popSize;

    // We're MINIMIZING here 
    maxScore  = -(population[0]).score;
    meanScore = -meanScore;
    minScore  = -(population[popSize - 1]).score;

    // Take care of the degenerate case where all scores are the same
    if (maxScore == minScore)
    {
        for (i = 0; i < popSize; i++)
        {
            (population[i]).expectation = 1/popSize;
        }
    }
    
    // Since we must sum to nParents, our mean must be this:
    desiredMean = 1/popSize;

    // We want to find a scale and an offset so that:
    // 1. scale * max + offset = MaximumSurvivalRate * desiredMean
    // and
    // 2.  Scale * mean + offset = desiredMean 
    // Subtracting 2 from 1, Factoring out scale and desiredMean, and dividing
    // both sides by max - mean gives:
    scale = desiredMean * (fitScaling->maximumSurvivalRate - 1) / (maxScore - meanScore);

    // offset so that the mean is desiredMean
    offset = desiredMean - (scale * meanScore);

    // if the above causes the least fitness to go negative,
    // change our goal to have a min of zero & mean of nParents/length(scores)
    if (offset + scale * minScore < 0)
    {
        scale = desiredMean / (meanScore - minScore);
        offset = desiredMean - (scale * meanScore);
    }

    for (i = 0; i < popSize; i++)
    {
      (population[i]).expectation = offset + scale*(population[i]).score;
    }    
}
 */

Gen *selectionRoulette(Gen *population, Selection *selection, mwSize popSize)
{
    mwIndex i;
    double position = randud(0, 1);
    double slotPosition = 0.0;

    for (i = 0; i < popSize; i++)
    {
        slotPosition = slotPosition + (population[i]).expectation;
        
        if (position < slotPosition) 
        {
            return &(population[i]);
        }
    }
    
    mexErrMsgTxt("Roulette selection method failed.");
}

Gen *selectionTournament(Gen *population, Selection *selection, mwSize popSize)
{
    mwIndex i;
    Gen  *winner, *player;

    // Assume an initial random winner player
    winner = &(population[ randui(0, popSize - 1) ]);
    // Loop for the rest of the players 
    for (i = 0; i < selection->nTournamentPlayers - 1; ++i)
    {
        // Get a new random player
        player = &(population[ randui(0, popSize - 1) ]);
        // If this player have a better fitness, 
        // replace winner
        if ( player->expectation > winner->expectation )
        {
            winner = player;
        }        
    }

    return winner;

}


void mutationGaussian( Gen *mutationChild, 
                       Gen *mutationParent, 
                       Mutation *mutation,
                       mwIndex generation,                        
                       mwIndex *classMask,                       
                       mwSize genomeLength) 
{
    mwIndex i;
    double scaleMax;
    
    scaleMax = (mutation->scale)*(1 - (mutation->shrink)*generation/(mutation->generations) );
    for (i = 0; i < genomeLength; i++){
        switch(classMask[i]) {
            case CLASS_DOUBLE: {                               
                mutationChild->chromosomes[i] = 
                mutationParent->chromosomes[i] + randnd(0.0, scaleMax*(mutation->span[i]));
                break;
            }
            case CLASS_BITSTR: {
                if (randud(0, 1) < mutation->rate){
                    mutationChild->chromosomes[i] = !mutationParent->chromosomes[i];
                } else {
                    mutationChild->chromosomes[i] = mutationParent->chromosomes[i];
                }
                break;
            }
            default:{
                mexErrMsgTxt("Unknow gen class.");
                break;
            }
        }        
    }
}

void crossoverScattered(Gen *xOverChild, Gen *xOverParentX, Gen *xOverParentY, mwSize genomeLength)
{
    mwIndex i;

    for (i = 0; i < genomeLength; i++)
    {
        if (randud(0, 1) < 0.5)
        {
            xOverChild->chromosomes[i] = xOverParentX->chromosomes[i];
        } 
        else 
        {
            xOverChild->chromosomes[i] = xOverParentY->chromosomes[i];
        }
    }
}



Gen *evolution(Ga *ga, Fitness *fitness)
{
  mwIndex i, j, k;
  mwSize     popSize = ga->popSize;    
  Gen        *mutationParent, *xOverParentY, *xOverParentX;
  Gen        *mutationChild, *xOverChild;
  mwIndex    genomeLength = ga->genomeLength;    
  Gen        *thisPopulation  = ga->thisPopulation;
  Gen        *nextPopulation  = ga->nextPopulation;
  FitScaling *fitScaling  = &(ga->fitScaling);
  Selection  *selection   = &(ga->selection);
  Mutation   *mutation    = &(ga->mutation);
  Crossover  *crossover   = &(ga->crossover);
  Stopping   *stopping    = &(ga->stopping);
  // Create a random initial population
  initPopulation(ga);
  // Main loop
  for (i = 0; i < stopping->nGenerations; i++)
  { //mexPrintf("\n generation = %d \n", i);
    // Update population         
    for (j = 0; j < popSize; j++)
    {
      moveGen( &(thisPopulation[j]), &(nextPopulation[j]), genomeLength);           
    }
    // Score population
    for (j = 0; j < popSize; j++)
    {          
      boundGen(&(thisPopulation[j]), ga->lowerBound, ga->upperBound, genomeLength);            
      fitness->fitnessFcn(genomeLength, &(thisPopulation[j]), fitness);            
    }
    // Sort the genes of ga in asscending order 
    // (Best genes first, we are minimizing)
    qsort(thisPopulation, popSize, sizeof(Gen), compareGen);
    // Scale ga score
    fitScaling->scalingFcn(thisPopulation, fitScaling, popSize);

    // Initialize ga iterator
    j = 0;

    // Initialize elite itertor
    k = 0;
    while (k < ga->nEliteKids)
    {
      moveGen( &(nextPopulation[j]), &(thisPopulation[j]), genomeLength);
      // Increment ga iterator  
      j++;
      // Increment elite iterator
      k++;
    }
    // Initialize crossover itertor
    k = 0;
    while (k < ga->nXoverKids)
    {
      // Get a crossover parent
      xOverParentX = selection->selectionFcn(thisPopulation, selection, popSize);
      // Get another crossover parent
      xOverParentY = selection->selectionFcn(thisPopulation, selection, popSize);
      // Select the crossover child
      xOverChild = &(nextPopulation[j]);
      // Perform the crossover
      crossover->crossoverFcn(xOverChild, xOverParentX, xOverParentY, genomeLength);
      // Increment ga iterator  
      j++;
      // Increment crossover iterator
      k++;
    }
    // Initialize mutation iterator
    k = 0;
    while (k < ga->nMutateKids)
    { 
      // Get a mutation parent
      mutationParent = selection->selectionFcn(thisPopulation, selection, popSize);
      // select the mutation child 
      mutationChild = &(nextPopulation[j]);
      // Perform the mutation
      mutation->mutationFcn(mutationChild, mutationParent, mutation, i, ga->classMask, genomeLength);
      // Increment ga iterator  
      j++;
      // Increment mutation iterator
      k++;
    } //printPopulation(ga);     
  }
  // Return
  return &nextPopulation[0];
}
