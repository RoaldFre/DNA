#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <stdbool.h>
#include "vmath.h"

typedef struct config
{
	double timeStep;        /* The timestep (dt) in the simulation */
	double measureInterval; /* Time between performing measurements */
	long   measureSamples;  /* Number of samples to measure */
	double measureWait;     /* Time to wait before starting measurement */
	int    numMonomers;     /* Number of monomers in the DNA strand TODO */
	double thermostatTemp;  /* Thermostat temperature */
	double thermostatTau;   /* Thermostat relaxation time */
	double langevinGamma;	/* Friction coefficient for Langevin dynamics */
	int    verbose;         /* Iterations between dumping info, <0 to disable */
	double worldSize;       /* World = periodic cube with edges of this length */
	double truncationLen;	/* Length at which potentials are truncated */
	int    numBoxes;	/* Number of boxes in each dimension */
	bool   render;          /* Whether or not to render the simulation */
	double framerate;       /* The desired framerate when rendering */
	double radius;          /* The radius of the particles to render */
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
	Strand *strands;
} World;

bool allocWorld(int numStrands, int numBoxes, double worldSize);
bool allocStrand(Strand *s, int numMonomers);
void fillWorld(void);
void freeWorld(void);
void freeStrand(Strand *strand);
void stepWorld(void);
void dumpStats(void);
void dumpEnergies(FILE *stream);
bool physicsCheck(void);
void forEveryParticle(void (*f)(Particle *p));
void forEveryParticleD(void (*f)(Particle *p, void *data), void *data);

extern World world;
extern Config config;

/* Current time in the simulation */
extern double sim_time;

#endif
