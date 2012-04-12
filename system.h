#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <stdbool.h>
#include "vmath.h"

typedef struct config
{
	double timeStep;        /* The timestep (dt) in the simulation */
	double thermostatTemp;  /* Thermostat temperature */
	double thermostatTau;   /* Thermostat relaxation time */
	double langevinGamma;	/* Friction coefficient for Langevin dynamics */
	double worldSize;       /* World = periodic cube with edges of this length */
	double truncationLen;	/* Length at which potentials are truncated */
	double saltConcentration; /* Na+ concentration in the environment */
} Config;

typedef enum particleType
{
	PHOSPHATE,
	SUGAR,
	BASE_A,
	BASE_T,
	BASE_C,
	BASE_G 
} ParticleType;

typedef struct particle
{
	double m; /* Mass */
	Vec3 pos; /* Position */
	Vec3 vel; /* Velocity */
	Vec3 F;   /* Force */
	ParticleType type;
	struct particle *prev, *next; /* Previous/Next particle in box */
} Particle;

typedef struct strand
{
	int numMonomers;
	Particle *Ps; /* Phosphates */
	Particle *Ss; /* Sugars */
	Particle *Bs; /* Bases */
	Particle *all; /* List of *all* particles */
	/*
	 * Indices work like this:
	 *
	 *    .  i
	 *   /|\
	 *    |      .
	 *    |      .
	 *    |      .
	 *    |      Ps[1]
	 *    |      |
	 *    |      |
	 *    |    5'|
	 *    |      Ss[1]------Bs[1]     <-- i=1
	 *    |    3'|
	 *    |      |  . . . . . . . . . . . . . . . . . . . . 
	 *    |      |                                       /|\
	 *    |      Ps[0]                                    |
	 *    |      |                                        |   one
	 *    |      |                                        |  monomer
	 *    |    5'|                                        |  
	 *    |      Ss[0]------Bs[0]     <-- i=0            \|/
	 *    |    3'   . . . . . . . . . . . . . . . . . . . .
	 *    |  
	 * 
	 */
} Strand;

typedef struct world
{
	int numStrands;
	int compStrand;
	Strand *strands;
} World;

/* Tasks are used to do stuff during the simulation, such as (accumulating) 
 * measurements. */
typedef struct task
{
	/* Data pointer that gets passed to the start() function below. Can 
	 * be used to pass configuration data. Make sure that this gets 
	 * allocated on the heap, so it will stay resident when we start 
	 * the simulation! */
	void *initialData;

	/* Called at the start of the simulation run. Returns the state 
	 * pointer that gets passed to the functions below. */
	void *(*start)(void *initialData);

	/* Called at every iteration step. Gets passed the data pointer 
	 * that start() returned. Returns false if the simulation has to be 
	 * stopped, true if everything can continue.
	 * Time gets updated after the call to tick(), ie the simulation 
	 * time when tick() gets called first is always 0. */
	bool (*tick)(void *state);

	/* Called at the end of the simulation run.  Gets passed the data 
	 * pointer that start() returned. */
	void (*stop)(void *state);
} Task;



/* Globals */
extern World world;
extern Config config;



/* Run the given task in the simulation */
void run(Task *task);

/* Get the time in the simulation after the last completed iteration. */
double getTime(void);

/* Get the number of the last completed iteration. */
long getIteration(void);

#endif
