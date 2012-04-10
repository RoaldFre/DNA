#ifndef _SAMPLERS_H_
#define _SAMPLERS_H_

#include "measure.h"

/* Sampler that averages the temperature (in a naive way, so don't use it 
 * on a huge amount of samples or you will loose precision). Prints the 
 * result at the end. */
Sampler averageTemperatureSampler(void);

/* Sampler that dumps physics stats. */
Sampler dumpStatsSampler(void);

/* Sampler that dumps the position of the particle. */
Sampler particlePositionSampler(Particle *p);

/* Sampler that dumps the Center Of Mass position of the strand. */
Sampler strandCOMSampler(Strand *s);

#endif
