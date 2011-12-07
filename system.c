#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "vmath.h"
#include "system.h"

/* Masses (in kg) */
#define AU      1.660539e-27
#define MASS_A  (134.1 * AU)
#define MASS_P  (94.97 * AU)
#define MASS_S  (83.11 * AU)

/* Equilibrium distance of bonds (in m) */
#define D_S5P   3.899e-10
#define D_S3P   3.559e-10
#define D_SA    6.430e-10

static void verlet(void);
static void calculateAcceleration(void);

static double randNorm(void);
static void randNormVec(double stdev, Vec3 *vec);


/* GLOBALS */

World world;
Config config;

double time = 0;


/* UNDER THE BONNET STUFF */


/* Allocates the world.
 * Precondition: config MUST be valid, and allocWorld must not already have 
 * been called (unless followed by a freeWorld)
 * Returns true on success, false on failure. In the case of failure, 
 * nothing will be allocated */
bool allocWorld()
{
	assert(world.Ps == NULL);
	assert(world.As == NULL);
	assert(world.Ss == NULL);

	world.Ps = calloc(config.numMonomers, sizeof(*world.Ps));
	world.As = calloc(config.numMonomers, sizeof(*world.As));
	world.Ss = calloc(config.numMonomers, sizeof(*world.Ss));
	if (world.Ps == NULL || world.As == NULL || world.Ss == NULL) {
		free(world.Ps);
		free(world.As);
		free(world.Ss);
		return false;
	}

	return true;
}

/* Place monomers in a vertical column (in the x-y plane) in the center of 
 * the world. Distances between sugar, base and phospate are the 
 * equilibrium lenghts with some small gaussian jitter added.
 *
 * TODO: not yet, now just simple gaussian sampling.
 * Initial velocities are sampled from a normal distribution to get a 
 * temperature of config.temperature 'per monomer'.
 * Velocities are then shifted to set the total momentum to zero. */

/*
 * Indices work like this:
 *
 *        .  y
 *       /|\
 *        |      .
 *        |      .
 *        |      .
 *        |      Ps[1]
 *        |      |
 *        |      |
 *        |    5'|
 *        |      Ss[1]------As[1]     <-- i=1
 *        |    3'|
 *        |      |  . . . . . . . . . . . . . . . . . . . . 
 *        |      |                                       /|\
 *        |      Ps[0]                                    |
 *        |      |                                        |   one
 *        |      |                                        |  monomer
 *        |    5'|                                        |  
 *        |      Ss[0]------As[0]     <-- i=0            \|/
 *        |    3'    . . . . . . . . . . . . . . . . . . .'
 *        |  
 *        |  
 *        | 
 *        +-----------------------------------------------------> x
 *       / 
 *      /
 *     /
 *    /
 *  |/   z 
 *  ''' 
 * 
 */
void fillWorld()
{
	int n = config.numMonomers;
	double ws = config.worldSize;
	Vec3 totP = {0, 0, 0}, avgP; /* momentum */
	Vec3 tmp;
	Vec3 corrVelS, corrVelA, corrVelP; /* correction on velocity to get P=0 */

	double spacing = D_S5P + D_S3P; /* vertical spacing between monomers */
	double yoffset = (ws - n * spacing) / 2;
	double xoffset = (ws - D_SA) / 2;
	double posStdev = spacing / 1000;
	double velStdev = sqrt(config.temperature); // TODO factor 1/3 (3D) somethere?

	for (int i = 0; i < n; i++) {
		/* Positions */
		world.Ss[i].pos.z = world.Ps[i].pos.z = world.As[i].pos.z = ws / 2;

		world.Ss[i].pos.x = xoffset;
		world.As[i].pos.x = xoffset + D_SA;
		world.Ps[i].pos.x = xoffset;

		world.Ss[i].pos.y = yoffset + i*spacing;
		world.As[i].pos.y = yoffset + i*spacing;
		world.Ps[i].pos.y = yoffset + i*spacing + D_S5P;

		world.Ss[i].pos.x += posStdev * randNorm();
		world.Ss[i].pos.y += posStdev * randNorm();
		world.Ss[i].pos.z += posStdev * randNorm();
		world.As[i].pos.x += posStdev * randNorm();
		world.As[i].pos.y += posStdev * randNorm();
		world.As[i].pos.z += posStdev * randNorm();
		world.Ps[i].pos.x += posStdev * randNorm();
		world.Ps[i].pos.y += posStdev * randNorm();
		world.Ps[i].pos.z += posStdev * randNorm();

		/* Velocities */
		randNormVec(velStdev, &world.Ss[i].vel);
		randNormVec(velStdev, &world.As[i].vel);
		randNormVec(velStdev, &world.Ps[i].vel);

		scale(&world.Ss[i].vel, MASS_S, &tmp);
		add(&totP, &tmp, &totP);
		scale(&world.As[i].vel, MASS_S, &tmp);
		add(&totP, &tmp, &totP);
		scale(&world.Ps[i].vel, MASS_S, &tmp);
		add(&totP, &tmp, &totP);
	}
	/* Correct for zero momentum */
	scale(&totP, 1.0 / (3 * n), &avgP);
	scale(&avgP, 1.0 / MASS_S, &corrVelS);
	scale(&avgP, 1.0 / MASS_A, &corrVelA);
	scale(&avgP, 1.0 / MASS_P, &corrVelP);
	for (int i = 0; i < config.numMonomers; i++) {
		sub(&world.Ss[i].vel, &corrVelS, &world.Ss[i].vel);
		sub(&world.As[i].vel, &corrVelA, &world.As[i].vel);
		sub(&world.Ps[i].vel, &corrVelP, &world.Ps[i].vel);
	}
}

/* Returns a number sampled from a standard normal distribution. */
static double randNorm()
{
	/* Box-Muller transform */
	double u1 = ((double) rand()) / RAND_MAX;
	double u2 = ((double) rand()) / RAND_MAX;

	return sqrt(-2*log(u1)) * cos(2*M_PI*u2);
}

static void randNormVec(double stdev, Vec3 *vec)
{
	vec->x = stdev * randNorm();
	vec->y = stdev * randNorm();
	vec->z = stdev * randNorm();
}

void freeWorld()
{
	free(world.Ss);
	free(world.As);
	free(world.Ps);

	return;
}

void dumpWorld()
{
	//TODO
	return;
}



/* PHYSICS */

static void verlet()
{
	double dt = config.timeStep;
	// TODO
}

static void calculateAcceleration()
{
	/* Reset acceleration */
	for (int i = 0; i < config.numMonomers; i++) {
		world.Ss[i].acc.x = 0;
		world.Ss[i].acc.y = 0;
		world.Ss[i].acc.z = 0;
		world.As[i].acc.x = 0;
		world.As[i].acc.y = 0;
		world.As[i].acc.z = 0;
		world.Ps[i].acc.x = 0;
		world.Ps[i].acc.y = 0;
		world.Ps[i].acc.z = 0;
	}
	// TODO
}

void stepWorld(void)
{
	verlet();
	time += config.timeStep;
}

void dumpStats()
{
	return;
}
