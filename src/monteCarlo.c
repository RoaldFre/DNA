#include "monteCarlo.h"
#include "spgrid.h"
#include "render.h"
#include <string.h>


/* TODO make part of MonteCarloState, also for all other tasks that use a 
 * renderString. */
#define RENDER_STRING_CHARS 80
static char renderString[RENDER_STRING_CHARS];


static void initializePreviousPosition(Particle *p)
{
	p->prevPos = p->pos;
}

/* One MC sweep 'equals' this time, so we can still use the measurement 
 * machinery, which relies on the simulated time to determine when to 
 * sample. */
#define SWEEP_TIME (1 * PICOSECONDS)

#define TARGET_ACCEPTANCE_RATIO 0.3 /* Aim for this acceptance ratio */
#define ACCEPTANCE_BUNCH 1000 /* update acceptance after this many moves */

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

#define INITIAL_NUM_PIVOT_POINTS_FRACTION (0.02)
#define INITIAL_PIVOT_ANGLE_STDDEV        (20 * DEGREE)
#define PIVOT_UPDATE_FRACTION 1.15

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
	Mat4 transfo;
} PivotPoint;

static PivotPoint generatePivotPoint(Strand *s, int pivotIndex, PivotConfig *cfg)
{
	Vec3 origin;
	Vec3 axis = randomDirection();	
	double theta = cfg->angleStdDev * randNorm();

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
	res.transfo = mat4rotationAround(origin, axis, theta);
	return res;
}


/* Indices: must be sorted ascendingly and not have repeating values. TODO: 
 * generalize to more than one (different) pivot move per monomer. */
static void generateAndApplyPivotMove(Strand *s, int *indices, int numIndices, 
		PivotConfig *cfg)
{
	int n = s->numMonomers;

	if (numIndices == 0)
		return;

	/* We can't deal with periodic boundary conditions in our 
	 * transformations! */
	undoPeriodicBoundaryConditions(s);

	Mat4 transfo = mat4identity();
	int pivPt = 0; /* index in the list of indices */

	for (int i = 0; i < n; i++) {
		/* Save original positions */
		s->Ss[i].prevPos = s->Ss[i].pos;
		s->Bs[i].prevPos = s->Bs[i].pos;
		s->Ps[i].prevPos = s->Ps[i].pos;

		/* Apply transformation */
		s->Ss[i].pos = mat4point(transfo, s->Ss[i].pos);
		s->Bs[i].pos = mat4point(transfo, s->Bs[i].pos);
		s->Ps[i].pos = mat4point(transfo, s->Ps[i].pos);

		/* If the current monomer is a pivot monomer, then same of 
		 * the transformations above will be redundant. 
		 * Nonetheless, it is important that we transformed the 
		 * pivot monomer itself with the 'old' transformation, 
		 * because otherwise the pivot origin will not be at the 
		 * correct position! */
		if (pivPt < numIndices && indices[pivPt] == i) {
			/* Current monomer *is* a pivot monomer! */
			PivotPoint pivotPoint = generatePivotPoint(s, indices[pivPt], cfg);
			switch (pivotPoint.type) {
			case PIVOT_PHOSPHATE:
				/* Entire monomer needed to be transformed 
				 * with the old transformation, which we 
				 * already did. Just update the accumulated 
				 * transfo. */
				transfo = mat4multiply(pivotPoint.transfo, transfo);
				break;
			case PIVOT_SUGAR:
				/* TODO 'half pivot' the base? (rotate 
				 * 'half' with old transfo and 'half' with 
				 * new transfo?) */
				s->Ps[i].pos = mat4point(pivotPoint.transfo, s->Ps[i].pos);
				transfo = mat4multiply(pivotPoint.transfo, transfo);
				break;
			case PIVOT_BASE:
				/* This transformation is only for the 
				 * base. Don't update the accumulated 
				 * transfo. */
				s->Bs[i].pos = mat4point(pivotPoint.transfo, s->Bs[i].pos);
				break;
			default:
				die("Internal error: unknown pivot type\n");
			}
			pivPt++;
		}
		/* Rebox after moving */
		reboxParticle(&s->Ss[i]);
		reboxParticle(&s->Bs[i]);
		reboxParticle(&s->Ps[i]);
	}
}

static void updatePivotParameters(double acc, double targetAcc, void *data)
{
	PivotConfig *cfg = (PivotConfig *) data;

	/* TODO do something more smart -- scale parameterFactor with targetAcc/acc? */
	double parameterFactor;
	if (acc > targetAcc)
		parameterFactor = PIVOT_UPDATE_FRACTION;
	else
		parameterFactor = 1.0 / PIVOT_UPDATE_FRACTION;

	cfg->angleStdDev *= parameterFactor;
	cfg->angleStdDev = MAX(cfg->angleStdDev, 1*DEGREE);

	cfg->maxNumPivPts *= parameterFactor;
	cfg->maxNumPivPts = MAX(cfg->maxNumPivPts, 1);
	cfg->maxNumPivPts = MIN(cfg->maxNumPivPts, cfg->s->numMonomers);

	//printf("Adjusted parameters: %f %f\n", cfg->angleStdDev / DEGREE, cfg->maxNumPivPts);
}

static void *initPivotMove(void)
{
	if (world.numStrands != 1)
		die("Expected a single strand in the world!\n");

	forEveryParticle(&initializePreviousPosition);

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
	int n = 1 + (int)(rand01() * cfg->maxNumPivPts);
	int *indices = calloc(n, sizeof(*indices));
	uniformSortedIndices(cfg->s->numMonomers, n, indices);

	generateAndApplyPivotMove(cfg->s, indices, n, cfg);

	//TODO don't realloc all the time
	free(indices);
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
	.updateParameters = &updatePivotParameters,
	.exit = &freePointer,
};





/* JIGGLE MOVE */
#define JIGGLE_UPDATE_FRACTION (1.2)
#define INITIAL_JIGGLE_LENGTH (0.01 * ANGSTROM)
static void jiggleWorker(Particle *p, void *data)
{
	double h = *(double*) data;
	p->prevPos = p->pos;
	p->pos = add(p->pos, randUniformVec(-h, h));
	reboxParticle(p);
}
static void undoJiggleWorker(Particle *p)
{
	p->pos = p->prevPos;
	reboxParticle(p);
}
static void jiggleMove(void *data)
{
	forEveryParticleD(&jiggleWorker, data);
}
static void undoJiggleMove(void *data)
{
	UNUSED(data);
	forEveryParticle(&undoJiggleWorker);
}
static void updateJiggleParameters(double acc, double targetAcc, void *data)
{
	double *h = (double*) data;

	/* TODO do something more smart -- scale parameterFactor with targetAcc/acc? */
	double parameterFactor;
	if (acc > targetAcc)
		parameterFactor = JIGGLE_UPDATE_FRACTION;
	else
		parameterFactor = 1.0 / JIGGLE_UPDATE_FRACTION;

	*h *= parameterFactor;

	//printf("Adjusted parameter: %f\n", *h / ANGSTROM);
}
static void *initJiggleMove(void)
{
	forEveryParticle(&initializePreviousPosition);
	double *h = malloc(sizeof(*h));
	*h = INITIAL_JIGGLE_LENGTH;
	return h;
}
MonteCarloMover jiggleMover = {
	.description = "Jiggle",
	.init = &initJiggleMove,
	.doMove = &jiggleMove,
	.undoMove = &undoJiggleMove,
	.updateParameters = &updateJiggleParameters,
	.exit = &freePointer,
};



/* JIGGLE SOME MOVE */
/* Looks like -- although this allows for larger jiggle lengths -- 
 * perturbing *all* particles is still more efficient due to the O(N) cost 
 * of getting the potential energy. */
#define JIGGLE_SOME_UPDATE_FRACTION (1.2)
#define JIGGLE_SOME_PARTICLE_FRACTION (0.1)
#define INITIAL_JIGGLE_SOME_LENGTH (0.05 * ANGSTROM)
typedef struct {
	double h; /* Jiggle length */
	Strand *s; /* Strand to jiggle */
	int n; /*  Number of particles to jiggle (must be <= s->numMonomers) */
} JiggleSomeConfig;
static void jiggleSomeMove(void *data)
{
	JiggleSomeConfig *jsc = (JiggleSomeConfig*) data;

	int *indices = calloc(jsc->n, sizeof(*indices));
	uniformSortedIndices(3 * jsc->s->numMonomers, jsc->n, indices);

	for (int i = 0; i < jsc->n; i++)
		jiggleWorker(&jsc->s->all[indices[i]], &jsc->h);
}
static void updateJiggleSomeParameters(double acc, double targetAcc, void *data)
{
	JiggleSomeConfig *jsc = (JiggleSomeConfig*) data;

	/* TODO do something more smart -- scale parameterFactor with targetAcc/acc? */
	double parameterFactor;
	if (acc > targetAcc)
		parameterFactor = JIGGLE_UPDATE_FRACTION;
	else
		parameterFactor = 1.0 / JIGGLE_UPDATE_FRACTION;

	jsc->h *= parameterFactor;

	//printf("Adjusted parameter: %f\n", jsc->h / ANGSTROM);
}
static void *initJiggleSomeMove(void)
{
	forEveryParticle(&initializePreviousPosition);
	if (world.numStrands != 1)
		die("Expected a single strand in the world!\n");

	JiggleSomeConfig *jsc = malloc(sizeof(*jsc));
	jsc->h = INITIAL_JIGGLE_SOME_LENGTH;
	jsc->s = &world.strands[0];
	jsc->n = MAX(1, 3*jsc->s->numMonomers * JIGGLE_SOME_PARTICLE_FRACTION);
	return jsc;
}
MonteCarloMover jiggleSomeMover = {
	.description = "JiggleSome",
	.init = &initJiggleSomeMove,
	.doMove = &jiggleSomeMove,
	.undoMove = &undoJiggleMove,
	.updateParameters = &updateJiggleSomeParameters,
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
	double timeStep; /* Time step for every MC move */
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
		acc += moves->moves[i++].weight;
	i--;
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

	advanceTimeBy(mcs->timeStep);
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
	state->timeStep = SWEEP_TIME / world.strands[0].numMonomers;

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

	/* De-initialize all movers */
	MonteCarloMoves *m = &mcs->moves;
	for (int i = 0; i < m->numMoves; i++)
		m->moves[i].m->exit(m->moves[i].m->data);

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

