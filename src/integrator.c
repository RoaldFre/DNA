#include <string.h>
#include <stdio.h>
#include "integrator.h"
#include "spgrid.h"

static double timeStep;

__inline__ double getIntegratorTimeStep(void)
{
	return timeStep;
}
__inline__ void setIntegratorTimeStep(double dt)
{
	timeStep = dt;
}


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
	double dt = getIntegratorTimeStep();
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
	double dt = getIntegratorTimeStep();

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
	double dt = getIntegratorTimeStep();

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
 * second order.
 * However, it also uses two standard normally distributed random 
 * vectors per iteration instead of one and has to store two extra 
 * vectors per particle instead of just one. */
/* Before calling this helper function, the values for position, velocity 
 * and force for each particle are evaluated at times:
 *   START:
 *     pos(t)
 *     vel(t - dt)
 *     F(t)
 *     fPrev = F(t - dt)/m
 *     rnd(t - dt)
 * This function propagates positions and velocities to
 *   END:
 *     pos(t + dt)
 *     vel(t)
 *     F(t)
 *     fPrev = F(t)/m
 *     rnd(t)
 * When calling calcForces() after this function, F(t + dt) gets calculated 
 * based on r(t + dt) and we are back at the starting position, but with t 
 * incremented by dt.
 */
static void langevin2helper(Particle *p, void *data)
{
	LangevinSettings *settings = (LangevinSettings*) data;

	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos,   "start langevin2helper");
		debugVectorSanity(p->vel,   "start langevin2helper");
		debugVectorSanity(p->F,     "start langevin2helper");
		debugVectorSanity(p->fPrev, "start langevin2helper");
		debugVectorSanity(p->rnd,   "start langevin2helper");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
		assert(isSaneVector(p->fPrev));
		assert(isSaneVector(p->rnd));
	}

	double dt  = getIntegratorTimeStep();
	double g   = settings->gamma;
	double T   = getHeatBathTemperature();
	double s   = sqrt(2 * BOLTZMANN_CONSTANT * T * g / p->m);
	double sdt = sqrt(dt);
	double gdt = g*dt;
	Vec3 f = scale(p->F, 1/p->m); /* f(t) */

	/* vel(t - dt) -> vel(t) */
	p->vel = add(add(scale(add(f,
	                           scale(p->fPrev, 1 - gdt)),
			       dt/2),
			p->rnd),
			scale(p->vel, (1 - gdt + gdt*gdt/2)));


	/* xi(t - dt) -> xi(t),   theta(t - dt) -> theta(t) */
	Vec3 xi = randNormVec(1);
	Vec3 theta = randNormVec(s*sdt*dt/(2*sqrt(3))); /* Common prefactor */

	/* Compute the random vector for the velocity propagation in the 
	 * next iteration.
	 * rnd(t - dt) -> rnd(t) */
	p->rnd = add(scale(xi, s*sdt*(1 - gdt/2)),
	             scale(theta, -g));
	/* The random vector for the position propagation below */
	Vec3 rndForPos = add(theta, scale(xi, sdt*sdt/2));


	/* pos(t) -> pos(t + dt) */
	p->pos = add(add(
			scale(add(scale(p->vel, (1 - gdt/2)),
			          scale(f, dt/2)),
			      dt),
			rndForPos),
			p->pos);


	/* Store f(t) */
	p->fPrev = f;


	if (DEBUG_VECTOR_SANITY) {
		debugVectorSanity(p->pos,   "end langevin2helper");
		debugVectorSanity(p->vel,   "end langevin2helper");
		debugVectorSanity(p->F,     "end langevin2helper");
		debugVectorSanity(p->fPrev, "end langevin2helper");
		debugVectorSanity(p->rnd,   "end langevin2helper");
	} else {
		assert(isSaneVector(p->pos));
		assert(isSaneVector(p->vel));
		assert(isSaneVector(p->F));
		assert(isSaneVector(p->fPrev));
		assert(isSaneVector(p->rnd));
	}
}

/* Alternative integrator for Langevin dynamics.
 * Based on the work of Vanden-Eijden and Ciccotti. Accurate up to dt^2. */
static void langevin2(LangevinSettings *settings)
{
	forEveryParticleD(&langevin2helper, settings);
	calculateForces();
}

#else //ALTERNATIVE_LANGEVIN

static void langevinBBKhelper(Particle *p, void *data)
{
	LangevinSettings *settings = (LangevinSettings*) data;

	double dt = getIntegratorTimeStep();
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
	tmp = scale(p->F, dt*dt / p->m);
	tmp = add(tmp, scale(p->prevPos, (g*dt/2 - 1)));
	tmp = add(tmp, scale(p->pos, 2));

	Vec3 newPos = scale(tmp, 1 / (1 + g*dt/2));

	p->vel = scale(sub(newPos, p->prevPos), 1/(2*dt));
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
 * Warning, only accurate up to first order! */
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

	advanceTimeBy(getIntegratorTimeStep());
}

#ifdef ALTERNATIVE_LANGEVIN
static void initializeParticle(Particle *p, void *data)
{
	LangevinSettings *settings = (LangevinSettings*) data;

	double dt  = getIntegratorTimeStep();
	double g   = settings->gamma;
	double T   = getHeatBathTemperature();
	double s   = sqrt(2 * BOLTZMANN_CONSTANT * T * g / p->m);
	double sdt = sqrt(dt);
	double gdt = g*dt;

	Vec3 xi = randNormVec(1);
	Vec3 theta = randNormVec(s*sdt*dt/(2*sqrt(3)));
	p->rnd = add(scale(xi, s*sdt*(1 - gdt/2)),
	             scale(theta, -g));

	p->fPrev = scale(p->F, 1/p->m);
}
#else
static void initializeParticle(Particle *p)
{
	p->prevPos = sub(p->pos, scale(p->vel, getIntegratorTimeStep()));
}
#endif

typedef struct
{
	Integrator integrator;
	double reboxInterval;
	double lastReboxTime;
} IntegratorState;
static void *integratorTaskStart(void *initialData)
{
	IntegratorConf *ic = (IntegratorConf*) initialData;

	setIntegratorTimeStep(ic->timeStep);

	calculateForces();
#ifdef ALTERNATIVE_LANGEVIN
	forEveryParticleD(&initializeParticle, &ic->integrator.settings.langevin);
#else
	forEveryParticle(&initializeParticle);
#endif

	IntegratorState *state = malloc(sizeof(*state));
	state->integrator = ic->integrator;
	state->reboxInterval = ic->reboxInterval;
	state->lastReboxTime = getTime();

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
				"# Gamma: %le\n", 
				conf->integrator.settings.langevin.gamma);
		break;
	case VERLET:
		integratorTypeAndSettings = asprintfOrDie(
				"# Integrator: Verlet\n"
				"# Tau: %le\n", 
				conf->integrator.settings.verlet.tau);
		break;
	default:
		die("Unknown integrator type!\n");
		return NULL; /* To avoid uninitialized usage warnings */
	}

	char *ret = asprintfOrDie("%s# Time step: %le\n# Rebox interval: %le\n",
			integratorTypeAndSettings,
			conf->timeStep, conf->reboxInterval);
	free(integratorTypeAndSettings);
	return ret;
}

