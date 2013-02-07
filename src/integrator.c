#include <string.h>
#include <stdio.h>
#include "integrator.h"
#include "spgrid.h"

static double heatBathTemperature;

__inline__ void setHeatBathTemperature(double temperature)
{
	heatBathTemperature = temperature;
}
__inline__ double getHeatBathTemperature(void)
{
	return heatBathTemperature;
}

static double timeStep;

__inline__ double getTimeStep(void)
{
	return timeStep;
}
__inline__ void setTimeStep(double dt)
{
	timeStep = dt;
}

static double simulationTime;

double getTime(void)
{
	return simulationTime;
}


/* ===== TEMPERATURE TASK ===== */
typedef struct {
	TemperatureTable table;
	/* Index of the next entry in the table that still needs to 
	 * executed. */
	int index;
} TemperatureState;

static TaskSignal temperatureTaskTick(void *state)
{
	TemperatureState *ts = (TemperatureState*) state;
	if (ts->index >= ts->table.numSetpoints)
		return TASK_OK;

	TemperatureSetpoint setpoint = ts->table.setpoints[ts->index];
	if (setpoint.time < getTime())
		return TASK_OK;

	setHeatBathTemperature(setpoint.temperature);
	ts->index++;

	return TASK_OK;
}

Task makeTemperatureTask(TemperatureTable table)
{
	Task task;
	TemperatureState *state = malloc(sizeof(*state));
	state->table = table;
	state->index = 0;

	task.initialData = state;
	task.start = &passPointer;
	task.tick  = &temperatureTaskTick;
	task.stop  = &freePointer;

	return task;
}



/* ===== INTEGRATOR ===== */

static void thermostatHelper(Particle *p, void *data)
{
	double lambda = *(double*) data;
	p->vel = scale(p->vel, lambda);
}
/* Berendsen thermostat */
static void thermostat(double tau)
{
	if (tau <= 0)
		return;

	/* Mass and Boltzmann constant are 1 */ 
	double Tk = getKineticTemperature();
	assert(isSaneNumber(Tk));
	double T0 = getHeatBathTemperature();
	double dt = getTimeStep();
	double lambda2 = 1 + dt/tau * (T0/Tk - 1);
	double lambda;
	if (lambda2 >= 0)
		lambda = sqrt(lambda2);
	else
		lambda = 1e-20; //TODO sane?

	forEveryParticleD(&thermostatHelper, (void*) &lambda);
}

static void verletHelper1(Particle *p)
{
	double dt = getTimeStep();

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos, "start verletHelper1");
		debugVectorSanity(p->vel, "start verletHelper1");
		debugVectorSanity(p->F,   "start verletHelper1");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
	}

	/* vel(t + dt/2) = vel(t) + acc(t)*dt/2 */
	p->vel = add(p->vel, scale(p->F, dt / (2 * p->m)));

	/* pos(t + dt) = pos(t) + vel(t + dt/2)*dt */
	p->pos = add(p->pos, scale(p->vel, dt));

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos, "end verletHelper1");
		debugVectorSanity(p->vel, "end verletHelper1");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
	}
}
static void verletHelper2(Particle *p)
{
	double dt = getTimeStep();

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos, "start verletHelper2");
		debugVectorSanity(p->vel, "start verletHelper2");
		debugVectorSanity(p->F,   "start verletHelper2");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
	}

	/* vel(t + dt) = vel(t + dt/2) + acc(t + dt)*dt/2 */
	p->vel = add(p->vel, scale(p->F, dt / (2*p->m)));

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->vel, "end verletHelper2");
	} else {
		assert(isSaneVector(p->vel));
	}
}
static void verlet(void)
{
	forEveryParticle(&verletHelper1);
	calculateForces(); /* acc(t + dt) */
	forEveryParticle(&verletHelper2);
}


static void langevinBBKhelper(Particle *p, void *data)
{
	LangevinSettings *settings = (LangevinSettings*) data;

	double dt = getTimeStep();
	double g  = settings->gamma;
	double T  = getHeatBathTemperature();

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos, "start langevinBBKhelper");
		debugVectorSanity(p->vel, "start langevinBBKhelper");
		debugVectorSanity(p->F,   "start langevinBBKhelper");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
	}

	/* Regular forces have been calculated. Add the random force due 
	 * to collisions to the total force. The result is:
	 * p->F = F(t + dt) + R(t + dt) */
	double Rstddev = sqrt(2 * BOLTZMANN_CONSTANT * T * g * p->m / dt);
	Vec3 R = randNormVec(Rstddev);
	debugVectorSanity(R, "randNormVec in langevinBBKhelper");
	p->F = add(p->F, R);

	Vec3 tmp;
	tmp = scale(p->pos, 2);
	tmp = add(tmp, scale(p->prevPos, (g*dt/2 - 1)));
	tmp = add(tmp, scale(p->F, dt*dt / p->m));
	Vec3 newPos = scale(tmp, 1 / (1 + g*dt/2));

	p->vel = scale(sub(newPos, p->pos), 1/dt);
	p->prevPos = p->pos;
	p->pos = newPos;

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos, "end langevinBBKhelper");
		debugVectorSanity(p->vel, "end langevinBBKhelper");
		debugVectorSanity(p->F,   "end langevinBBKhelper");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
	}
}

/* BBK integrator for Langevin dynamics.
 * See http://localscf.com/LangevinDynamics.aspx 
 * (relocated to http://localscf.com/localscf.com/LangevinDynamics.aspx.html atm) */
static void langevinBBK(LangevinSettings *settings)
{
	calculateForces();
	forEveryParticleD(&langevinBBKhelper, settings);
}



static void stepPhysics(Integrator integrator)
{
	assert(worldSanityCheck());

	syncPhysics();

	switch(integrator.type) {
	case VERLET:
		verlet();
		assert(physicsCheck());
		thermostat(integrator.settings.verlet.tau);
		assert(physicsCheck());
		break;
	case LANGEVIN:
		langevinBBK(&integrator.settings.langevin);
		break;
	default:
		fprintf(stderr, "ERROR: Unknown integrator!\n");
		assert(false);
		break;
	}

	simulationTime += timeStep;
}

static void initializePreviousPosition(Particle *p)
{
	p->prevPos = sub(p->pos, scale(p->vel, timeStep));
}

typedef struct
{
	Integrator integrator;
	double reboxInterval;
	double lastReboxTime;
} IntegratorState;
/* The integrator task is responsible for handeling the space partition 
 * grid */
static void *integratorTaskStart(void *initialData)
{
	IntegratorConf *ic = (IntegratorConf*) initialData;

	forEveryParticle(&initializePreviousPosition);

	IntegratorState *state = malloc(sizeof(*state));
	state->integrator = ic->integrator;
	state->reboxInterval = ic->reboxInterval;
	state->lastReboxTime = 0;

	setTimeStep(ic->timeStep);
	simulationTime = 0;

	free(initialData);
	return state;
}
static TaskSignal integratorTaskTick(void *state)
{
	IntegratorState *is = (IntegratorState*) state;
	stepPhysics(is->integrator);

	if (getTime() > is->lastReboxTime + is->reboxInterval) {
		reboxParticles();
		is->lastReboxTime = getTime();
	}

	return TASK_OK;
}
static void integratorTaskStop(void *state)
{
	UNUSED(state);
	freeGrid();
}

Task makeIntegratorTask(IntegratorConf *conf)
{
	Task task;
	IntegratorConf *confCpy = malloc(sizeof(*confCpy));
	memcpy(confCpy, conf, sizeof(*confCpy));

	task.initialData = confCpy;
	task.start = &integratorTaskStart;
	task.tick  = &integratorTaskTick;
	task.stop  = &integratorTaskStop;
	return task;
}

char *integratorInfo(IntegratorConf *conf)
{
	char *integratorTypeAndSettings;

	switch (conf->integrator.type) {
	case LANGEVIN:
		integratorTypeAndSettings = asprintfOrDie(
				"# Integrator: Langevin\n"
				"# Gamma: %e\n", 
				conf->integrator.settings.langevin.gamma);
		break;
	case VERLET:
		integratorTypeAndSettings = asprintfOrDie(
				"# Integrator: Verlet\n"
				"# Tau: %e\n", 
				conf->integrator.settings.verlet.tau);
		break;
	default:
		die("Unknown integrator type!\n");
		return NULL; /* To avoid uninitialized usage warnings */
	}

	char *ret = asprintfOrDie("%s# Time step: %e\n # Rebox interval: %e\n",
			integratorTypeAndSettings,
			conf->timeStep, conf->reboxInterval);
	free(integratorTypeAndSettings);
	return ret;
}

