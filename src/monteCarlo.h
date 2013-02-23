#ifndef _MONTE_CARLO_H_
#define _MONTE_CARLO_H_

#include "physics.h"

/* Request a negative number of monteCarloSweeps to keep on going 
 * indefinitely. */
Task makeMonteCarloTask(int monteCarloSweeps);

#endif
