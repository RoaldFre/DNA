#ifndef _WORLD_H_
#define _WORLD_H_
#include "system.h"

bool allocWorld(int numStrands, double worldSize);
bool allocStrand(Strand *s, int numMonomers);
void fillStrand(Strand *s, const char *sequence);
void fillComplementaryStrand(Strand *s, const char *baseSequence);
void freeWorld(void);
void freeStrand(Strand *strand);

int numParticles(void);
void forEveryParticle(void (*f)(Particle *p));
void forEveryParticleD(void (*f)(Particle *p, void *data), void *data);
void forEveryParticleOf(Strand *s, void (*f)(Particle *p));
void forEveryParticleOfD(Strand *s,
			void (*f)(Particle *p, void *data), void *data);

#endif
