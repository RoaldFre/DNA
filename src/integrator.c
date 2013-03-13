#include <string.h>
#include <stdio.h>
#include "integrator.h"
#include "spgrid.h"

static double timeStep;

__inline__ double getTimeStep(void)
{
	return timeStep;
}
__inline__ void setTimeStep(double dt)
{
	timeStep = dt;
}

static double simulationTime = 0;

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


#ifdef ALTERNATIVE_LANGEVIN
/* Alternative Langevin integrator. This one should be accurate up to 
 * second order. It also uses two standard normally distributed random 
 * vectors per iteration instead of one and has to store three extra 
 * vectors per particle. Moreover, it needs two passes, with one 
 * calculateForces() pass in between middle.
 * Only use this for testing purposes, as it is slower (and 
 * less accurate[?]) than the regular integrator below. */
static void langevin2helper1(Particle *p, void *data)
{
	LangevinSettings *settings = (LangevinSettings*) data;

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos,   "start langevin2helper1");
		debugVectorSanity(p->vel,   "start langevin2helper1");
		debugVectorSanity(p->F,     "start langevin2helper1");
		debugVectorSanity(p->fPrev, "start langevin2helper1");
		debugVectorSanity(p->xi,    "start langevin2helper1");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
		assert(isSaneVector(p->fPrev));
		assert(isSaneVector(p->xi));
	}

	double dt  = getTimeStep();
	double sdt = sqrt(dt);
	double g   = settings->gamma;
	double T   = getHeatBathTemperature();
	double s   = sqrt(2 * BOLTZMANN_CONSTANT * T * g / p->m);

	Vec3 f = scale(p->F, 1/p->m); /* f(r(t)) */
	Vec3 gv = scale(p->vel, g);
	Vec3 xi = randNormVec(1);
	Vec3 theta = randNormVec(1);
	Vec3 R = add(scale(xi, 1/2.0), scale(theta, 1/(2.0 * sqrt(3.0))));
	Vec3 A = add(scale(sub(f, gv), SQUARE(dt)/2.0),  scale(R, s*dt*sdt));

	/* r(t) -> r(t + dt) */
	p->pos = add(add(p->pos, scale(p->vel, dt)), A);
	p->fPrev = f; /* f(r(t)) */
	p->xi = xi;   /* xi(t) */
	p->A = A;     /* A(t) */



	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos,   "end langevin2helper1");
		debugVectorSanity(p->vel,   "end langevin2helper1");
		debugVectorSanity(p->F,     "end langevin2helper1");
		debugVectorSanity(p->fPrev, "end langevin2helper1");
		debugVectorSanity(p->xi,    "end langevin2helper1");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
		assert(isSaneVector(p->fPrev));
		assert(isSaneVector(p->xi));
	}
}
static void langevin2helper2(Particle *p, void *data)
{
	LangevinSettings *settings = (LangevinSettings*) data;

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos,   "start langevin2helper2");
		debugVectorSanity(p->vel,   "start langevin2helper2");
		debugVectorSanity(p->F,     "start langevin2helper2");
		debugVectorSanity(p->fPrev, "start langevin2helper2");
		debugVectorSanity(p->xi,    "start langevin2helper2");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
		assert(isSaneVector(p->fPrev));
		assert(isSaneVector(p->xi));
	}

	double dt  = getTimeStep();
	double sdt = sqrt(dt);
	double g   = settings->gamma;
	double T   = getHeatBathTemperature();
	double s   = sqrt(2 * BOLTZMANN_CONSTANT * T * g / p->m);

	Vec3 f = scale(p->F, 1/p->m); /* f(r(t + dt)) */
	Vec3 gv = scale(p->vel, g);

	/* v(t) -> v(t + dt) */
	p->vel = add(add(add(add(
			p->vel,
			scale(add(p->fPrev, f), dt/2)),
			scale(gv, -dt)),
			scale(p->xi, s*sdt)),
			scale(p->A, -g));

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos,   "end langevin2helper2");
		debugVectorSanity(p->vel,   "end langevin2helper2");
		debugVectorSanity(p->F,     "end langevin2helper2");
		debugVectorSanity(p->fPrev, "end langevin2helper2");
		debugVectorSanity(p->xi,    "end langevin2helper2");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
		assert(isSaneVector(p->fPrev));
		assert(isSaneVector(p->xi));
	}
}

/* Alternative integrator for Langevin dynamics.
 * Based on the work of Vanden-Eijden and Ciccotti. Accurate up to dt^2. */
static void langevin2(LangevinSettings *settings)
{
	forEveryParticleD(&langevin2helper1, settings);
	calculateForces();
	forEveryParticleD(&langevin2helper2, settings);
}

#else //ALTERNATIVE_LANGEVIN

static void langevinBBKhelper(Particle *p, void *data)
{
	LangevinSettings *settings = (LangevinSettings*) data;

	double dt = getTimeStep();
	double g  = settings->gamma;
	double T  = getHeatBathTemperature();

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos,     "start langevinBBKhelper");
		debugVectorSanity(p->vel,     "start langevinBBKhelper");
		debugVectorSanity(p->F,       "start langevinBBKhelper");
		debugVectorSanity(p->prevPos, "start langevinBBKhelper");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
		assert(isSaneVector(p->prevPos));
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
		debugVectorSanity(p->pos,     "end langevinBBKhelper");
		debugVectorSanity(p->vel,     "end langevinBBKhelper");
		debugVectorSanity(p->F,       "end langevinBBKhelper");
		debugVectorSanity(p->prevPos, "end langevinBBKhelper");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
		assert(isSaneVector(p->prevPos));
	}
}

/* BBK integrator for Langevin dynamics.
 * See http://localscf.com/LangevinDynamics.aspx 
 * (relocated to http://localscf.com/localscf.com/LangevinDynamics.aspx.html atm)
 * Should be accurate up to dt^3 in the positions (and dt^2 in velocities) */
static void langevinBBK(LangevinSettings *settings)
{
	calculateForces();
	forEveryParticleD(&langevinBBKhelper, settings);
}
#endif //ALTERNATIVE_LANGEVIN




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
#ifdef ALTERNATIVE_LANGEVIN
		langevin2(&integrator.settings.langevin);
#else
		langevinBBK(&integrator.settings.langevin);
#endif
		break;
	default:
		fprintf(stderr, "ERROR: Unknown integrator!\n");
		assert(false);
		break;
	}

	simulationTime += timeStep;
}

static void initializeParticle(Particle *p)
{
	p->prevPos = sub(p->pos, scale(p->vel, timeStep));
#ifdef ALTERNATIVE_LANGEVIN
	p->xi = vec3(0, 0, 0);
	p->fPrev = vec3(0, 0, 0);
#endif
}

typedef struct
{
	Integrator integrator;
	double reboxInterval;
	double lastReboxTime;
} IntegratorState;
static void *integratorTaskStart(void *initialData)
{
	IntegratorConf *ic = (IntegratorConf*) initialData;

	forEveryParticle(&initializeParticle);

	IntegratorState *state = malloc(sizeof(*state));
	state->integrator = ic->integrator;
	state->reboxInterval = ic->reboxInterval;
	state->lastReboxTime = getTime();

	setTimeStep(ic->timeStep);

#ifdef ALTERNATIVE_LANGEVIN
	printf("Built for alternative Langevin integrator\n");
#endif

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

	char *ret = asprintfOrDie("%s# Time step: %e\n# Rebox interval: %e\n",
			integratorTypeAndSettings,
			conf->timeStep, conf->reboxInterval);
	free(integratorTypeAndSettings);
	return ret;
}

