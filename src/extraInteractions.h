#ifndef _EXTRA_INTERACTIONS_H_
#define _EXTRA_INTERACTIONS_H_

#include "world.h"

/* For umbrella sampling. The extra potential energy is K/2*(R - Rref)^2, 
 * with R the end-to-end distance of the strand */
typedef struct {
	double K;
	double Rref;
	Strand *s;
} HarmonicEndToEndInt;

void registerHarmonicEndToEndInt(HarmonicEndToEndInt *conf);

/* This is a suitable header to put at the top of an end-to-end 
 * measurement. It dumps the config parameters as non-comment in two 
 * columns. You need to free the returned pointer afterwards. */
char *harmonicEndToEndIntHeader(HarmonicEndToEndInt *conf);

#endif
