#ifndef _WORLD_H_
#define _WORLD_H_
#include "system.h"
#include "spgridBootstrap.h"

typedef enum particleType
{
	BASE_A = 0,
	BASE_T,
	BASE_C,
	BASE_G,
	BASE_X,
	BASE_Y,
	PHOSPHATE,
	SUGAR,

	/* WARNING:
	 * If for some reason you want to change the ordering of 
	 * this enum, watch out to avoid breaking getConnectedParticle(), 
	 * isBase() and feelExclusion()
	 * WARNING:
	 * If for some reason you change the base types to not be numbered 
	 * 0 to 3, then also change the dihedral angle cacheing code in 
	 * physics.c.
	 * WARNING:
	 * Finally, if for some reason you change the number of base types, 
	 * also change the define below, pretty much everything in the code 
	 * that involves particle types... */
} ParticleType;

#define NUM_BASE_TYPES 6 /* Number of base types in the enum above. */
#define NUM_PARTICLE_TYPES 8 /* Number of particle types in the enum above. */

static __inline__ bool isBase(ParticleType t)
{
	return t < PHOSPHATE;
}

typedef struct particle
{
	double m; /* Mass */
	Vec3 pos; /* Position */
	Vec3 prevPos; /* Position at previous time step. This gets managed 
			 by the integrator task.
			 NOTE: even with ALTERNATIVE_LANGEVIN, we still 
			 need this for the monte carlo dynamics!
			 NOTE2: translations in world.c also know of this, 
			 TODO nicer way? */
	Vec3 vel; /* Velocity */
	Vec3 F;   /* Force */
#ifdef ALTERNATIVE_LANGEVIN
	Vec3 fPrev; /* Force/mass at previous time step. Only needed for 
		       alternative Langevin integrator. */
	Vec3 rnd; /* Random component of force at previous time step. */
#endif
	ParticleType type;
	struct particle *prev, *next; /* Previous/Next particle in box */
	struct strand *strand; /* The strand I belong to */
	int strandIndex; /* My index in the appropriate array of my strand */
	struct box *myBox; /* The space patition box that I am in */
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
	double worldSize; /* World = periodic cube with edges of this 
			     length. Change with resizeWorld(). */
	int numStrands;
	Strand *strands;
} World;


/* Sets the heat bath temperature. This automatically calls syncPhysics() 
 * afterwards. */
void setHeatBathTemperature(double temperature);
double getHeatBathTemperature(void);

extern World world;

bool allocWorld(int numStrands, double worldSize);
bool allocStrand(Strand *s, int numMonomers);
/* Allocates the strand and fills it up with the given base sequence. The 
 * strand will be <3'-sequence-5'>, ie the 3' matches the first char of 
 * the string, and the last char of the sequence is the 5' end */
void fillStrand(Strand *s, const char *sequence);
void fillComplementaryStrand(Strand *s, const char *baseSequence);
/* Returns the base sequence of this strand. Don't forget to free the 
 * pointer afterwards! */
char *getSequence(Strand *s);
/* Returns a string containing information of the world. Each line is 
 * prefixed with '#'. You should free the returned pointer when you don't 
 * need it any more. */
char *getWorldInfo(void);
void freeWorld(void);
void freeStrand(Strand *strand);

/* Use this to resize the world: it automatically handles periodic boundary 
 * conditions correctly. */
void resizeWorld(double newWorldSize);

/* Write the current state of the world to the given file. Can be read with 
 * readWorld(). If filename == NULL, don't do anything. */
void writeWorld(const char *filename);
/* Read, allocate and fill the world with the state from the given file. 
 * The file must be produced by writeWorld(). */
void readWorld(const char *filename);

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

/* Translate the strand over the given vector.
 * WARNING: if calling this with space partitioning, you need to rebox the 
 * particles yourself! */
void translateStrand(Strand *s, Vec3 delta);

/* Center the given strand in the world so that its center of mass is the 
 * zero vector. Periodic boundary wrapping gets undone prior to shifting.
 * WARNING: if calling this with space partitioning, you need to rebox the 
 * particles yourself! */
void centerStrand(Strand *s);

/* If the particle is outside the world: translate it back. Also updates 
 * the prevPos correctly. */
void periodicPosition(Particle *p);

/* Undoes periodic boundary conditions. The sugar of the first monomer of 
 * the strand is used as an anchor.
 * The prevPos gets updated to keep the same difference (pos - prevPos) as 
 * before changing the positions */
void undoPeriodicBoundaryConditions(Strand *s);
/* This undoes undoPeriodicBoundaryConditions(). It won't change the boxes 
 * of the space partitioning grid, as it is assumed that these were correct 
 * before undoPeriodicBoundaryConditions() was called! */
void redoPeriodicBoundaryConditions(Strand *s);

/* Checks whether particles of strands are properly associated with their 
 * type, monomer and strand. */
bool strandSanityCheck(Strand *s);

/* Checks all the strands in the world */
bool worldSanityCheck(void);

#endif
