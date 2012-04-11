#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <stdbool.h>
#include "system.h"
#include "vmath.h"
#include "task.h"

typedef enum
{
	LANGEVIN,
	VERLET
} Integrator;

typedef struct
{
	Integrator integrator;
	int numBoxes; /* For space partition grid */
} IntegratorConf;

bool allocWorld(int numStrands, double worldSize);
bool allocStrand(Strand *s, int numMonomers);
void fillStrand(Strand *s, const char *sequence);
void freeWorld(void);
void freeStrand(Strand *strand);
void dumpStats(void);
void dumpEnergies(FILE *stream);
bool physicsCheck(void);
void forEveryParticle(void (*f)(Particle *p));
void forEveryParticleD(void (*f)(Particle *p, void *data), void *data);
void forEveryParticleOf(Strand *s, void (*f)(Particle *p));
void forEveryParticleOfD(Strand *s,
			void (*f)(Particle *p, void *data), void *data);

double temperature(void);
/* Returns the position vector of the Center Of Mass. */
Vec3 getCOM(Particle *ps, int num);

Task makeIntegratorTask(IntegratorConf *conf);

#endif
