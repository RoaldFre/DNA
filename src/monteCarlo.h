#ifndef _MONTE_CARLO_H_
#define _MONTE_CARLO_H_

#include "system.h"

typedef struct monteCarloMover MonteCarloMover;

typedef struct MonteCarloMove {
	/* The mover for this type of move. */
	MonteCarloMover *m;

	/* Weight of this move: how likely is it to be selected, compared 
	 * to other moves. Only relative values w.r.t other weights in the 
	 * same list of MonteCarloMoves are relevant. */
	double weight;
} MonteCarloMove;

typedef struct {
	MonteCarloMove *moves;
	int numMoves;
} MonteCarloMoves;


/* Pivot mover */
extern MonteCarloMover pivotMover;



typedef struct {
	int sweeps; /* Negative to go on indefinitely */
	bool verbose;
	MonteCarloMoves moves;
} MonteCarloConfig;

Task makeMonteCarloTask(MonteCarloConfig *config);

#endif
