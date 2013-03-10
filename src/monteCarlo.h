#ifndef _MONTE_CARLO_H_
#define _MONTE_CARLO_H_

#include "system.h"

typedef struct monteCarloMover MonteCarloMover;

typedef struct MonteCarloMove {
	/* Weight of this move: how likely is it to be selected, compared 
	 * to other moves. Only relative values w.r.t other weights in the 
	 * same list of MonteCarloMoves are relevant. */
	double weight;

	/* The mover for this type of move. */
	MonteCarloMover *m;
} MonteCarloMove;

typedef struct {
	MonteCarloMove *moves;
	int numMoves;
} MonteCarloMoves;


/* Pivot mover: Randomly pivot the strand around random positions. */
extern MonteCarloMover pivotMover;

/* Jiggle mover: Randomly perturb the positions of the particles in the strand. */
extern MonteCarloMover jiggleMover;



typedef struct {
	int sweeps; /* Negative to go on indefinitely */
	bool verbose;
	MonteCarloMoves moves;
} MonteCarloConfig;

Task makeMonteCarloTask(MonteCarloConfig *config);

#endif
