#ifndef fitness_h
#define fitness_h

#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "ga.h"
#include "utils.h"

typedef struct Fitness
{
	void (*fitnessFcn)(mwSize, struct Gen *, struct Fitness *);
	mwSize  pLinksQ;   // Cantidad de enlaces primarios
  mwSize  sLinksQ;   // Cantidad de enlaces secundarios
	double  *channels; // Matriz de coeficientes del canal
	double  *sinrThr;  // SINR requerida por cada enlace  
	double  *sinr;     // SINR de cada enlace
	double  *power;    // Potencia de transmision de cada enlace
	double  *noise;    // Potencia de ruido de cada enlace
} Fitness;

void decode(double *, mwSize , struct Gen *);
void fitnessFcn(mwSize , struct Gen *, struct Fitness *);
void fitParamSet(struct Fitness *, mxArray *);

#endif
