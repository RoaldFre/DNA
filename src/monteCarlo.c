#include "monteCarlo.h"
#include "spgrid.h"
#include "render.h"
#include <string.h>


/* TODO make part of MonteCarloState, also for all other tasks that use a 
 * renderString. */
#define RENDER_STRING_CHARS 80
static char renderString[RENDER_STRING_CHARS];


#define TARGET_ACCEPTANCE_RATIO 0.3 /* Aim for this acceptance ratio */
#define ACCEPTANCE_BUNCH 5000 /* update acceptance after this many moves */

/* Acceptance stats of the current and previous period. A period is defined 
 * as ACCEPTANCE_BUNCH moves */
typedef struct {
	long attempted; /* Number of attempted moves in this period */
	long accepted;  /* Number of accepted moves in this period */
	float prevAcceptance; /* The acceptance of the previous period */
} Acceptance;

struct monteCarloMover {
	/* Short description of the type of Monte Carlo move. */
	const char *description;

	/* Accumulate acceptance data here */
	Acceptance acceptance;

	/* Data for the Monte Carlo move (e.g. parameters) */
	void *data;
	
	/* Set up and return the data pointer */
	void *(*init)(void);

	/* Perform a move. */
	void (*doMove)(void *data);

	/* Undo the last move. */
	void (*undoMove)(void *data);

	/* Update the parameters of this type of Monte Carlo move to reach 
	 * the target acceptance ratio. */
	void (*updateParameters)(double acceptance, double targetAcceptance,
								void *data);

	/* De-initialize */
	void (*exit)(void *data);
};


/* PIVOT MOVE */

#define INITIAL_NUM_PIVOT_POINTS_FRACTION (0.05)
#define INITIAL_PIVOT_ANGLE_STDDEV        (30 * DEGREE)

typedef struct {
	Strand *s; /* The strand to perturb */
	double maxNumPivPts; /* up to floor(maxNumPivPts) pivot points */
	double angleStdDev; /* std deviation of the generated pivot angle */
} PivotConfig;

typedef enum {
	PIVOT_PHOSPHATE = 0, /* pivot around the phosphate */
	PIVOT_SUGAR,         /* pivot around the sugar */
	PIVOT_BASE,          /* pivot only the base */
} PivotType;
#define NUM_PIVOT_TYPES 3

typedef struct {
	PivotType type;
	int pivotIndex;
	Mat4 transfo;
} PivotPoint;

/* Multiple pivot points for a single strand. Pivot points need to be 
 * sorted with ascending pivotIndex.
 *
 * There can only be at most one pivotPoint for a specific index!
 * TODO generalize */
typedef struct {
	int numPivotPoints;
	PivotPoint* pivotPoints;
} PivotChain;

/* XXX WARNING: due to numerical errors, the bond spacing at the end of 
 * long chains will start to vary significantly. This will get corrected by 
 * the acceptance probability when bond interactions are still enabled. It 
 * could give a slightly biassed sampling, though.
 * TODO this only is the case for long pivot chains? */
static void applyPivotChain(Strand *s, PivotChain chain)
{
	int n = s->numMonomers;
	int npts = chain.numPivotPoints;
	PivotPoint *pts = chain.pivotPoints;

	if (npts == 0)
		return;

	/* We can't deal with periodic boundary conditions in our 
	 * transformations! */
	undoPeriodicBoundaryConditions(s);

	Mat4 transfo = mat4identity();
	int pivPt = 0;

	for (int i = 0; i < n; i++) {
		/* Save original positions */
		s->Ss[i].prevPos = s->Ss[i].pos;
		s->Bs[i].prevPos = s->Bs[i].pos;
		s->Ps[i].prevPos = s->Ps[i].pos;

		/* Apply transformation */
		if (pivPt >= npts || pts[pivPt].pivotIndex != i) {
			/* Current monomer is not a special pivot monomer */
			s->Ss[i].pos = mat4point(transfo, s->Ss[i].pos);
			s->Bs[i].pos = mat4point(transfo, s->Bs[i].pos);
			s->Ps[i].pos = mat4point(transfo, s->Ps[i].pos);
		} else {
			/* Current monomer *is* a pivot monomer */
			switch (pts[pivPt].type) {
			case PIVOT_PHOSPHATE:
				/* Pivot entire monomer with old transfo */
				s->Ss[i].pos = mat4point(transfo, s->Ss[i].pos);
				s->Bs[i].pos = mat4point(transfo, s->Bs[i].pos);
				s->Ps[i].pos = mat4point(transfo, s->Ps[i].pos);
				transfo = mat4multiply(pts[pivPt].transfo, transfo);
				break;
			case PIVOT_SUGAR:
				/* TODO 'half pivot' the base? (rotate 
				 * 'half' with old transfo and 'half' with 
				 * new transfo?) */
				s->Ss[i].pos = mat4point(transfo, s->Ss[i].pos);
				s->Bs[i].pos = mat4point(transfo, s->Bs[i].pos);
				transfo = mat4multiply(pts[pivPt].transfo, transfo);
				s->Ps[i].pos = mat4point(transfo, s->Ps[i].pos);
				break;
			case PIVOT_BASE:
				s->Ss[i].pos = mat4point(transfo, s->Ss[i].pos);
				Mat4 baseTransfo = mat4multiply(pts[pivPt].transfo, transfo);
				s->Bs[i].pos = mat4point(baseTransfo, s->Bs[i].pos);
				s->Ps[i].pos = mat4point(transfo, s->Ps[i].pos);
				break;
			default:
				die("Internal error: unknown pivot type\n");
			}
			pivPt++;
		}
		reboxParticle(&s->Ss[i]);
		reboxParticle(&s->Bs[i]);
		reboxParticle(&s->Ps[i]);
	}
}

static PivotPoint generatePivotPoint(Strand *s, int pivotIndex)
{
	Vec3 origin;
	Vec3 axis = randomDirection();	
	//double theta = 60*DEGREE * (2*rand01() - 1);
	double theta = 10*DEGREE * randNorm();

	PivotType type = NUM_PIVOT_TYPES * rand01();

	switch (type) {
	case PIVOT_PHOSPHATE:
		/* Pivot around the phosphate */
		origin = s->Ps[pivotIndex].pos;
		break;
	case PIVOT_SUGAR:
		/* Pivot around the sugar */
		origin = s->Ss[pivotIndex].pos;
		break;
	case PIVOT_BASE:
		/* Pivot the base *around* the sugar! */
		origin = s->Ss[pivotIndex].pos;
		break;
	default:
		die("Internal error: unknown pivot type\n");
		origin = vec3(0, 0, 0); /* avoid compiler warning */
	}

	PivotPoint res;
	res.type = type;
	res.pivotIndex = pivotIndex;
	res.transfo = mat4rotationAround(origin, axis, theta);
	return res;
}


/* Generates a new pivot chain (with chain->numPivotPoints points). The 
 * given pivotChain must have enough space allocated for these points. */
static void generatePivotChain(PivotChain *chain, Strand *s)
{
	//TODO don't realloc all the time?
	int *indices = calloc(chain->numPivotPoints, sizeof(*indices));
	uniformSortedIndices(s->numMonomers, chain->numPivotPoints, indices);
	for (int i = 0; i < chain->numPivotPoints; i++)
		chain->pivotPoints[i] = generatePivotPoint(s, indices[i]);
	free(indices);
}

static void *initPivotMove(void)
{
	if (world.numStrands != 1)
		die("Expected a single strand in the world!\n");

	PivotConfig *cfg = malloc(sizeof(*cfg));
	cfg->s = &world.strands[0];
	int n = cfg->s->numMonomers;
	cfg->maxNumPivPts = MAX(1, n * INITIAL_NUM_PIVOT_POINTS_FRACTION);
	cfg->angleStdDev = INITIAL_PIVOT_ANGLE_STDDEV;

	return cfg;
}

static void pivotMove(void *data)
{
	PivotConfig *cfg = (PivotConfig *) data;
	PivotChain chain;
	int n = 1 + randIndex((int) cfg->maxNumPivPts);
	chain.numPivotPoints = n;
	chain.pivotPoints = calloc(n, sizeof(*chain.pivotPoints));

	generatePivotChain(&chain, cfg->s);
	applyPivotChain(cfg->s, chain);

	free(chain.pivotPoints); //TODO don't realloc all the time
}

static void undoPivotMove(void *data)
{
	Strand *s = ((PivotConfig*) data)->s;
	for (int i = 0; i < 3*s->numMonomers; i++) {
		s->all[i].pos = s->all[i].prevPos;
		reboxParticle(&s->all[i]);
	}
}

MonteCarloMover pivotMover = {
	.description = "Pivot",
	.init = &initPivotMove,
	.doMove = &pivotMove,
	.undoMove = &undoPivotMove,
	.updateParameters = NULL, //TODO
	.exit = &freePointer,
};





/* GENERIC MONTE CARLO MOVE DRIVER */

static bool acceptMove(double *previousPotentialEnergy)
{
	double V = getPotentialEnergy();
	double dE = V - *previousPotentialEnergy;
	double beta = 1.0/(BOLTZMANN_CONSTANT * getHeatBathTemperature());
	if (dE < 0 || rand01() < exp(-dE * beta)) {
		/* Move is accepted */
		*previousPotentialEnergy = V;
		return true;
	} else {
		/* Move is rejected */
		return false;
	}
}

/* Update acceptance stats. If we exceed ACCEPTANCE_BUNCH, then update acc 
 * accordingly and call the update function (if it isn't NULL) with the 
 * given data */
static void updateAcceptance(Acceptance *acc, bool accepted,
			void (*update)(double, double, void*), void *data)
{
	acc->attempted++;
	if (accepted)
		acc->accepted++;

	if (acc->attempted < ACCEPTANCE_BUNCH)
		return;

	acc->prevAcceptance = ((float) acc->accepted) / acc->attempted;
	acc->accepted = 0;
	acc->attempted = 0;

	if (update != NULL)
		update(acc->prevAcceptance, TARGET_ACCEPTANCE_RATIO, data);
}

typedef struct {
	long numMoves; /* Number of (attempted) moves up till now */
	long maxMoves; /* Maximum number of (attempted) moves */
	Acceptance totalAcceptance; /* Total acceptance of last few moves */
	double previousPotentialEnergy;
	MonteCarloMoves moves; /* List of possible moves */
	bool verbose; /* Make some noise? */
} MonteCarloState;

static void monteCarloMove(MonteCarloState *mcs)
{
	MonteCarloMoves *moves = &mcs->moves;
	
	/* Select a move based on their weights */
	double totWeight = 0;
	for (int i = 0; i < moves->numMoves; i++)
		totWeight += moves->moves[i].weight;
	double rnd = totWeight * rand01();
	double acc = 0;
	int i = 0;
	while (i < moves->numMoves  &&  rnd > acc)
		acc += moves->moves[i].weight;
	MonteCarloMover *m = moves->moves[i].m;

	m->doMove(m->data);
	bool accepted = acceptMove(&mcs->previousPotentialEnergy);
	if (!accepted)
		m->undoMove(m->data);
	updateAcceptance(&m->acceptance, accepted, m->updateParameters,
								m->data);
	updateAcceptance(&mcs->totalAcceptance, accepted, NULL, NULL);

	snprintf(renderString, RENDER_STRING_CHARS, "acceptance: %f",
					mcs->totalAcceptance.prevAcceptance);
}



/* TASK */

static void *monteCarloTaskStart(void *initialData)
{
	MonteCarloConfig *mcc = (MonteCarloConfig*) initialData;
	if (world.numStrands != 1)
		die("Expected a single strand in the world!\n");

	/* Set up state */
	MonteCarloState *state = malloc(sizeof(*state));
	state->numMoves = 0;
	state->maxMoves = mcc->sweeps * world.strands[0].numMonomers;
	state->verbose = mcc->verbose;
	state->previousPotentialEnergy = getPotentialEnergy();
	state->moves = mcc->moves;

	/* Initialize movers */
	MonteCarloMoves *m = &state->moves;
	for (int i = 0; i < m->numMoves; i++)
		m->moves[i].m->data = m->moves[i].m->init();

	/* Rendering stuff */
	renderString[0] = '\0';
	RenderStringConfig rsc;
	rsc.string = renderString;
	rsc.x = 10;
	rsc.y = 60;
	registerString(&rsc);

	if (mcc->verbose) {
		if (state->maxMoves > 0)
			printf("Monte carlo move 0 / %ld (0%%)", state->maxMoves);
		else
			printf("Monte carlo move 0");
		fflush(stdout);
	}

	free(initialData);
	return state;
}
static TaskSignal monteCarloTaskTick(void *state)
{
	MonteCarloState *mcs = (MonteCarloState*) state;
	mcs->numMoves++;
	if (mcs->maxMoves >= 0  &&  mcs->numMoves > mcs->maxMoves)
		return TASK_STOP;

	monteCarloMove(mcs);

	if (mcs->verbose && 0 == mcs->numMoves % MAX(1, mcs->maxMoves / 100)) {
		if (mcs->maxMoves > 0)
			printf("\rMonte carlo move %ld / %ld (%ld%%)",
					mcs->numMoves, mcs->maxMoves,
					100 * mcs->numMoves / mcs->maxMoves);
		else
			printf("\rMonte carlo move %ld", mcs->numMoves);
		fflush(stdout);
	}

	return TASK_OK;
}
static void monteCarloTaskStop(void *state)
{
	MonteCarloState *mcs = (MonteCarloState*) state;
	if (mcs->verbose)
		printf("\n");
	free(state);
}
Task makeMonteCarloTask(MonteCarloConfig *config)
{
	MonteCarloConfig *mccCopy = malloc(sizeof(*mccCopy));
	memcpy(mccCopy, config, sizeof(*mccCopy));

	Task task;
	task.initialData = mccCopy;
	task.start = &monteCarloTaskStart;
	task.tick  = &monteCarloTaskTick;
	task.stop  = &monteCarloTaskStop;
	return task;
}

