#ifndef _EXTRA_INTERACTIONS_H_
#define _EXTRA_INTERACTIONS_H_

#include "world.h"

/* For umbrella sampling. The extra potential energy is K/2*(R - Rref)^2, 
 * with R the end-to-end distance of the strand. */
typedef struct {
	double K;
	double Rref;
	Strand *s;
} HarmonicEndToEndInt;

void registerHarmonicEndToEndInt(HarmonicEndToEndInt *conf);

/* You need to free the returned pointer afterwards. */
char *harmonicEndToEndIntHeader(HarmonicEndToEndInt *conf);



/* Apply a constant end-to-end force. */
typedef struct {
	Vec3 F; /* Force F at first monomer, -F at last monomer */
	Strand *s;
} EndToEndForceInt;

void registerEndToEndForceInt(EndToEndForceInt *conf);

/* You need to free the returned pointer afterwards. */
char *endToEndForceIntHeader(EndToEndForceInt *conf);

#endif
