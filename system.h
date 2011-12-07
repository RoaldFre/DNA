#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <stdbool.h>
#include "vmath.h"

typedef struct config
{
	double timeStep;        /* The timestep (dt) in the simulation */
	int    numMonomers;     /* Number of monomers in the DNA strand */
	double temperature;     /* The temperature timestep (dt) in the simulation */
	int    verbose;         /* Iterations between dumping info, <0 to disable */
	double worldSize;       /* The world is a cube with edges of this length 
				   (Only used for rendering) */
	bool   render;          /* Whether or not to render the simulation */
	int    renderSteps;     /* Physics steps between rendering frames. */
	double radius;          /* The radius of the particles to render */
} Config;

typedef struct particle
{
	Vec3 pos; /* Position */
	Vec3 vel; /* Velocity */
	Vec3 acc; /* Acceleration */
} Particle;


typedef struct world
{
	Particle *Ps; /* Phosphates */
	Particle *Ss; /* Sugars */
	Particle *As; /* bases -- currently only Adenine */
} World;


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
 *    |      Ss[1]------As[1]     <-- i=1
 *    |    3'|
 *    |      |  . . . . . . . . . . . . . . . . . . . . 
 *    |      |                                       /|\
 *    |      Ps[0]                                    |
 *    |      |                                        |   one
 *    |      |                                        |  monomer
 *    |    5'|                                        |  
 *    |      Ss[0]------As[0]     <-- i=0            \|/
 *    |    3'   . . . . . . . . . . . . . . . . . . . .
 *    |  
 * 
 */

bool allocWorld(void);
void fillWorld(void);
void freeWorld(void);
void stepWorld(void);
void dumpWorld(void);
void dumpStats(void);
bool sanityCheck(void);
Vec3 nearestImageVector(Vec3 *v1, Vec3 *v2); 

extern World world;
extern Config config;

/* Current time in the simulation */
extern double time;

#endif
