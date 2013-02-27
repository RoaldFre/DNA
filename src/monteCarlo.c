#include "monteCarlo.h"
#include "spgrid.h"
#include "render.h"
#include <string.h>

#define PIVOT_SELECTION_PROBABILITY 0.05

#define RENDER_STRING_CHARS 80
static char renderString[RENDER_STRING_CHARS];
static double previousPotentialEnergy;

#define ACCEPTANCE_BUNCH 10000 /* update acceptance after this many moves */
static int totalAttemptedMoves = 0;
static int totalAcceptedMoves = 0;
static float acceptance = 0;

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


/* The given pivotChain must have enough allocated space for at least 
 * s->numMonomers pivot points! */
static void generatePivotChain(PivotChain *chain, Strand *s,
					double selectionProbability)
{
	int numPivPts = 0;
	for (int i = 0; i < s->numMonomers; i++) {
		/* TODO this is probably wasteful on random numbers. 
		 * Especially for low selectionProbability.
		 * Better: random number between 1 and numMonomers, then 
		 * generate that many pivot points with pivotIndices 
		 * sampled uniformly over subsets of permutations of 
		 * 1..numMonomers (how to do that fast?) */
		if (rand01() > selectionProbability)
			continue;

		/* Generate new pivot point */
		chain->pivotPoints[numPivPts] = generatePivotPoint(s, i);
		numPivPts++;
	}
	chain->numPivotPoints = numPivPts;
}

static void pivotMove(Strand *s)
{
	PivotChain chain;
	//TODO don't realloc all the time
	chain.pivotPoints = calloc(s->numMonomers, sizeof(*chain.pivotPoints));

	generatePivotChain(&chain, s, PIVOT_SELECTION_PROBABILITY);
	applyPivotChain(s, chain);

	free(chain.pivotPoints);
}

#if 0
/* Move a random single particle */
static void jiggleMove(void)
{
	Strand *s = &world.strands[0];
	int n = s->numMonomers;
	int i = randIndex(3*n);
	//Vec3 displacement = randUniformVec(-4*ANGSTROM, 4*ANGSTROM);
	Vec3 displacement = randNormVec(2*ANGSTROM);

	s->all[i].pos = add(s->all[i].pos, displacement);
	double V2 = getPotentialEnergy();

	double dE = V2 - previousPotentialEnergy;
	double beta = 1.0/(BOLTZMANN_CONSTANT * getHeatBathTemperature());
	if (dE < 0 || rand01() < exp(-dE * beta)) {
		/* Move is accepted -> update boxes of moved particle! */
		reboxParticle(&s->all[i]);
		previousPotentialEnergy = V2;
		totalAcceptedMoves++;
	} else {
		/* Move is rejected -> reset back to original configuration! */
		s->all[i].pos = sub(s->all[i].pos, displacement);
	}
}
#endif

static bool acceptMove(void)
{
	double V = getPotentialEnergy();
	double dE = V - previousPotentialEnergy;
	double beta = 1.0/(BOLTZMANN_CONSTANT * getHeatBathTemperature());
	if (dE < 0 || rand01() < exp(-dE * beta)) {
		/* Move is accepted */
		previousPotentialEnergy = V;
		return true;
	} else {
		/* Move is rejected */
		return false;
	}
}

static void undoPivotMove(Strand *s)
{
	for (int i = 0; i < 3*s->numMonomers; i++) {
		s->all[i].pos = s->all[i].prevPos;
		reboxParticle(&s->all[i]);
	}
}

static void monteCarloMove(void)
{
	totalAttemptedMoves++;
	if (totalAttemptedMoves > ACCEPTANCE_BUNCH) {
		acceptance = ((float) totalAcceptedMoves)
				/ (totalAttemptedMoves - 1);
		totalAttemptedMoves = 1;
		totalAcceptedMoves = 0;
	}

	Strand *s = &world.strands[0];
	pivotMove(s);
	if (acceptMove()) {
		totalAcceptedMoves++;
	} else {
		undoPivotMove(s);
	}

	snprintf(renderString, RENDER_STRING_CHARS, "acceptance: %f",
								acceptance);
}

typedef struct {
	long numMoves;
	long maxMoves;
	bool verbose;
} MonteCarloState;
static void *monteCarloTaskStart(void *initialData)
{
	MonteCarloConfig *mcc = (MonteCarloConfig*) initialData;
	if (world.numStrands != 1)
		die("Expected a single strand in the world!\n");

	previousPotentialEnergy = getPotentialEnergy();

	renderString[0] = '\0';
	RenderStringConfig rsc;
	rsc.string = renderString;
	rsc.x = 10;
	rsc.y = 60;
	registerString(&rsc);

	MonteCarloState *state = malloc(sizeof(*state));
	state->numMoves = 0;
	state->maxMoves = mcc->sweeps * world.strands[0].numMonomers;
	state->verbose = mcc->verbose;

	if (mcc->verbose) {
		printf("Monte carlo move 0 / %ld (0%%)", state->maxMoves);
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

	monteCarloMove();

	if (mcs->verbose && 0 == mcs->numMoves % MAX(1, mcs->maxMoves / 100)) {
		printf("\rMonte carlo move %ld / %ld (%ld%%)",
					mcs->numMoves, mcs->maxMoves,
					100 * mcs->numMoves / mcs->maxMoves);
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

