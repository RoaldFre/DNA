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
	BASE_G,
	/* WARNING, if for some reason you want to change the ordering of 
	 * this enum, watch out to avoid breaking getConnectedParticle() ! */
} ParticleType;

typedef struct particle
{
	double m; /* Mass */
	Vec3 pos; /* Position */
	Vec3 vel; /* Velocity */
	Vec3 F;   /* Force */
	ParticleType type;
	struct particle *prev, *next; /* Previous/Next particle in box */
	struct strand *strand; /* The strand I belong to */
	int strandIndex; /* My index in the appropriate array of my strand */
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

/* Returns the connected particle to the given one, as determined by the 
 * strand it's in. Returns NULL if no such particle exists.
 * Particles are ordered to avoid doubles: this has the pleasant 
 * consequence that there is at most one particle that can be returned.
 * Depending on the given particle, the returned particle is:
 *
 * given     | returned
 * ----------+-----------------------------------------------------
 * Phosphate | Sugar of same monomer
 * Base      | Sugar of same monomer
 * Sugar     | Phosphate of previous monomer or NULL if no such monomer
 */
Particle *getConnectedParticle(Particle *p);


/* Checks whether particles of strands are properly associated with their 
 * type, monomer and strand. */
bool strandSanityCheck(Strand *s);

/* Checks all the strands in the world */
bool worldSanityCheck(void);

#endif
