#ifndef _MONTE_CARLO_H_
#define _MONTE_CARLO_H_

#include "physics.h"

typedef struct {
	int sweeps; /* Negative to go on indefinitely */
	bool verbose;
} MonteCarloConfig;
Task makeMonteCarloTask(MonteCarloConfig *config);

#endif
