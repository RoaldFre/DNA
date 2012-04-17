#ifndef _WORLD_H_
#define _WORLD_H_
#include "system.h"


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
	double worldSize; /* World = periodic cube with edges of this length */
			//TODO: worldsize only really relvant for spgrid, though
	int numStrands;
	Strand *strands;
} World;

extern World world;

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
