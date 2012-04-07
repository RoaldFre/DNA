#ifndef _SAMPLERS_H_
#define _SAMPLERS_H_

#include "measure.h"

/* Simple sampler that averages the temperature (in a naive way, so don't 
 * use it on a huge amount of samples or you will loose precision). Prints 
 * the result at the end. */
Sampler averageTemperatureSampler(void);

#endif
