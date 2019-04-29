#ifndef ga_h
#define ga_h

#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include "fitness.h"
#include "utils.h"


#define STRLEN 50
#define CLASS_DOUBLE 0
#define CLASS_BITSTR 1

typedef struct Gen
{
    double  *chromosomes;
    double  *score;
    double	expectation;
    bool	  feasible;    
} Gen;

typedef struct FitScaling
{
	// Pointer to scaling fucntion
    void (*scalingFcn) (struct Gen *, struct FitScaling *, mwSize);	
    // Used by FITSCALINGTOP scaling function
	double  quantity;
	// Used by FITSCALINGSHIFTLINEAR scaling function
	double	maximumSurvivalRate;

} FitScaling;

typedef struct Selection
{
	// Pointer to selection fucntion
    Gen * (*selectionFcn) (struct Gen *, struct Selection *, mwSize);	
    // Used by SELECTIONTOURNAMENT scaling function
    mwSize nTournamentPlayers;
} Selection;

typedef struct Crossover
{
	// Pointer to crossover fucntion
    void (*crossoverFcn) (struct Gen *, struct Gen *, struct Gen *, mwSize );    
} Crossover;

typedef struct Mutation
{
	// Pointer to mutation fucntion
    void (*mutationFcn) (struct Gen *, struct  Gen *, struct Mutation *, mwIndex, mwIndex *, mwSize );
    // Used by MUTATIONGAUSSIAN mutation function
	double rate;	
	double shrink;
	double scale;
	mwSize generations;
	double *span; 

} Mutation;

typedef struct Stopping
{
	mwSize nGenerations;
} Stopping;

typedef struct Ga
{
	struct Gen         *thisPopulation;
	struct Gen         *nextPopulation;
	struct FitScaling  fitScaling;
	struct Mutation    mutation;
	struct Crossover   crossover;
	struct Selection   selection;
	struct Stopping    stopping;
	mwIndex			   *classMask;
	double             *lowerBound;
	double			   *upperBound;
	mwSize             genomeLength;
  mwSize            nScores;
	mwSize   	       popSize;
	mwSize		       nEliteKids;
	mwSize             nXoverKids;
	mwSize             nMutateKids;
	
} Ga;

void gaParamSet(struct Ga *, mxArray *);
int compareGen(const void *, const void *);
int compareGen2(const void *, const void *);
void sortPopulation(struct Ga *);
void fitScalingRank(struct Gen *, struct FitScaling *, mwSize);
void fitScalingTop(struct Gen *, struct FitScaling *, mwSize);
//void fitScalingShiftLinear(struct Gen *, struct FitScaling *, mwSize);
Gen *selectionRoulette(struct Gen *, struct Selection *, mwSize);
Gen *selectionTournament(struct Gen *, struct Selection *, mwSize);
void mutationGaussian(struct Gen *, struct  Gen *, struct Mutation *, mwIndex, mwIndex *, mwSize );
void crossoverScattered(struct Gen *, struct  Gen *, struct  Gen *, mwSize );
Gen *evolution(struct Ga *, struct Fitness *);
void moveGen(Gen *, Gen *, mwSize );

#endif
