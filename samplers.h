#ifndef _SAMPLERS_H_
#define _SAMPLERS_H_

#include "measure.h"

/* Sampler that averages the temperature (in a naive way, so don't use it 
 * on a huge amount of samples or you will loose precision). Prints the 
 * result at the end. */
Sampler averageTemperatureSampler(void);

/* Sampler that dumps physics stats. */
Sampler dumpStatsSampler(void);

/* Sampler that dumps the position (in fm) of the particle. */
Sampler particlePositionSampler(Particle *p);
/* Sampler that dumps the Center Of Mass position (in fm) of the strand. */
Sampler strandCOMSampler(Strand *s);

/* Sampler that dumps the squared displacement (in fm^2) of the particle. */
Sampler particleSquaredDisplacementSampler(Particle *p);
/* Sampler that dumps the squared displacement of Center Of Mass (in fm^2) 
 * the strand. */
Sampler strandCOMSquaredDisplacementSampler(Strand *s);

/* Sampler that counts the number of base pair bindings. A pair is 'bound' 
 * if its base pair potential is lower than the given threshold. */
Sampler basePairingSampler(double energyThreshold);
#endif
