#include "monteCarlo.h"
#include "spgrid.h"

/* Rotates the monomers with index strictly higher than pivotMonomer of the 
 * strand over the given axis and angle. */
static void pivotStrandNextMonomers(Strand *s, int pivotMonomer,
		Vec3 origin, Vec3 axis, double theta)
{
	assert(s != NULL);
	assert(0 <= pivotMonomer && pivotMonomer < s->numMonomers);

	int n = s->numMonomers;

	for (int i = pivotMonomer + 1; i < n; i++) {
		s->Ps[i].pos = rotateAround(s->Ps[i].pos, origin, axis, theta);
		s->Ss[i].pos = rotateAround(s->Ss[i].pos, origin, axis, theta);
		s->Bs[i].pos = rotateAround(s->Bs[i].pos, origin, axis, theta);
	}
}
/* The pivot origin is the phosphate of the 'pivotMonomer' monomer of the 
 * given strand. */
static void pivotStrandPhosphate(Strand *s, int pivotMonomer,
						Vec3 axis, double theta)
{
	Vec3 origin = s->Ps[pivotMonomer].pos;
	pivotStrandNextMonomers(s, pivotMonomer, origin, axis, theta);
}
/* The pivot origin is the sugar of the 'pivotMonomer' monomer of the 
 * given strand. */
static void pivotStrandSugar(Strand *s, int pivotMonomer,
						Vec3 axis, double theta)
{
	Vec3 origin = s->Ss[pivotMonomer].pos;
	/* Also pivot the phosphate of the pivot monomer */
	int i = pivotMonomer;
	s->Ps[i].pos = rotateAround(s->Ps[i].pos, origin, axis, theta);

	pivotStrandNextMonomers(s, pivotMonomer, origin, axis, theta);
}
typedef enum {
	PIVOT_PHOSPHATE = 0, /* pivot around the phosphate */
	PIVOT_SUGAR,         /* pivot around the sugar */
} PivotType;
#define NUM_PIVOT_TYPES 2
static void pivotStrand(Strand *s, int pivotMonomer, PivotType type,
						Vec3 axis, double theta)
{
	switch (type) {
	case PIVOT_PHOSPHATE:
		pivotStrandPhosphate(s, pivotMonomer, axis, theta);
		break;
	case PIVOT_SUGAR:
		pivotStrandSugar(s, pivotMonomer, axis, theta);
		break;
	}
}
/* Also reboxes the pivot Monomer itself, so this can be used for pivots 
 * around phosphates and sugars. */
static void reboxAfterPivot(Strand *s, int pivotMonomer)
{
	assert(s != NULL);
	assert(0 <= pivotMonomer && pivotMonomer < s->numMonomers);

	int n = s->numMonomers;
	for (int i = pivotMonomer; i < n; i++) {
		reboxParticle(&s->Ps[i]);
		reboxParticle(&s->Ss[i]);
		reboxParticle(&s->Bs[i]);
	}
}

static void pivotMove(void)
{
	Strand *s = &world.strands[0];
	int pivot = (s->numMonomers - 1) * rand01();
	Vec3 axis = randomDirection();	
	//double theta = 60*DEGREE * (2*rand01() - 1);
	double theta = 30*DEGREE * randNorm();
	PivotType type = NUM_PIVOT_TYPES * rand01();

	double V1 = getPotentialEnergy();
	pivotStrand(s, pivot, type, axis, theta);
	double V2 = getPotentialEnergy();

	double dE = V2 - V1;
	double beta = 1.0/(BOLTZMANN_CONSTANT * getHeatBathTemperature());
	if (dE < 0 || rand01() < exp(-dE * beta)) {
		/* Move is accepted -> update boxes of moved particles! */
		reboxAfterPivot(s, pivot);
	} else {
		/* Move is rejected -> reset back to original configuration! */
		pivotStrand(s, pivot, type, axis, -theta);
	}
}

/* Move a random single particle */
static void jiggleMove(void)
{
	Strand *s = &world.strands[0];
	int n = s->numMonomers;
	int i = randIndex(3*n);
	//Vec3 displacement = randUniformVec(-4*ANGSTROM, 4*ANGSTROM);
	Vec3 displacement = randNormVec(2*ANGSTROM);

	double V1 = getPotentialEnergy();
	s->all[i].pos = add(s->all[i].pos, displacement);
	double V2 = getPotentialEnergy();

	double dE = V2 - V1;
	double beta = 1.0/(BOLTZMANN_CONSTANT * getHeatBathTemperature());
	if (dE < 0 || rand01() < exp(-dE * beta)) {
		/* Move is accepted -> update boxes of moved particle! */
		reboxParticle(&s->all[i]);
	} else {
		/* Move is rejected -> reset back to original configuration! */
		s->all[i].pos = sub(s->all[i].pos, displacement);
	}
}

static void monteCarloMove(void)
{
	// TODO: cache potential energy of previous move?
	if (rand01() < 0.2)
		pivotMove();
	else
		jiggleMove();
}

static TaskSignal monteCarloTaskTick(void *state)
{
	UNUSED(state);
	monteCarloMove();
	return TASK_OK;
}

Task makeMonteCarloTask(void)
{
	Task task;

	task.start = NULL;
	task.tick  = &monteCarloTaskTick;
	task.stop  = NULL;
	return task;
}

