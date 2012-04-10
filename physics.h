#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <stdbool.h>
#include "system.h"
#include "vmath.h"
#include "task.h"

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
void forEveryParticleOf(Strand *s, void (*f)(Particle *p));
void forEveryParticleOfD(Strand *s,
			void (*f)(Particle *p, void *data), void *data);

double temperature(void);
/* Returns the position vector of the Center Of Mass. */
Vec3 getCOM(Particle *ps, int num);

extern Task integratorTask;

#endif
