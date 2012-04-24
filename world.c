#include "main.h" //for die TODO clean
#include "world.h"
#include "spgrid.h"
#include "physics.h" /* for spacings/angles/... */
#include <string.h>


/* Structure of the helix. See table I in Knotts */

#define HELIX_DELTA_Z   (3.38*A) 	/* vertical spacing between layers */
#define HELIX_DELTA_PHI (36*RAD)	/* twist at each consecutive layer */

/* Angles */
#define P_PHI (94.038*RAD)
#define S_PHI (70.197*RAD)
#define A_PHI (41.905*RAD)
#define T_PHI (86.119*RAD)
#define C_PHI (85.027*RAD)
#define G_PHI (40.691*RAD)
/* Radial distance */
#define P_R (8.916*A)
#define S_R (6.981*A)
#define A_R (0.773*A)
#define T_R (2.349*A)
#define C_R (2.296*A)
#define G_R (0.828*A)
/* Heights */
#define P_Z (2.186*A)
#define S_Z (1.280*A)
#define A_Z (0.051*A)
#define T_Z (0.191*A)
#define C_Z (0.187*A)
#define G_Z (0.053*A)


World world;

/* Allocates the world to hold the given number of strands, and sets this 
 * value in world.
 * Precondition: allocWorld must not already have been called (unless 
 * followed by a freeWorld). */
bool allocWorld(int numStrands, double worldSize)
{
	assert(world.strands == NULL);
	world.strands = calloc(numStrands, sizeof(*world.strands));
	if (world.strands == NULL)
		return false;
	world.worldSize = worldSize;
	world.numStrands = numStrands;
	return true;
}

/* Returns true on success, false on failure. In the case of failure, 
 * nothing will be allocated */
bool allocStrand(Strand *s, int numMonomers) {
	/* Allocate one big continuous list */
	s->all = calloc(3 * numMonomers, sizeof(*s->all));
	if (s->all == NULL)
		return false;
	
	/* Split the list in three sublists */
	s->Ss = &s->all[0];
	s->Bs = &s->all[1 * numMonomers];
	s->Ps = &s->all[2 * numMonomers];

	s->numMonomers = numMonomers;
	return true;
}


/* 
 * Builds a single DNA strand with the given base sequence.
 * The strand is a single helix (in the y direction) in the center
 * of the world. Distances between sugar, base and phospate are the 
 * equilibrium lenghts with some small gaussian jitter added.
 *
 * The given strand may not have been allocated before. It is also your 
 * responsibility to free the strand afterwards (either via freeStrand() or 
 * freeWorld()).
 *
 * Flags:
 * complementarySequence: Use the complement of the given base sequence.
 * complementaryHelix:    Build the helix turning the other way around.
 *
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
 *        |      Ss[1]------Bs[1]     <-- i=1
 *        |    3'|
 *        |      |  . . . . . . . . . . . . . . . . . . . . 
 *        |      |                                       /|\
 *        |      Ps[0]                                    |
 *        |      |                                        |   one
 *        |      |                                        |  monomer
 *        |    5'|                                        |  
 *        |      Ss[0]------Bs[0]     <-- i=0            \|/
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
static void fillStrandHelper(Strand *s, const char *baseSequence,
		bool complementarySequence, bool complementaryHelix)
{
	//TODO complementaryHelix;

	int n = strlen(baseSequence);
	allocStrand(s, n);

	double ws = world.worldSize;
	double velStdev = sqrt(config.thermostatTemp);

	Vec3 offset;
	offset.x = ws / 2;
	offset.z = ws / 2;
	offset.y = (ws - n * HELIX_DELTA_Z) / 2;

	double phi = 0;
	double z = 0;
	for (int i = 0; i < n; i++) {
		/* Default to Adenine */
		ParticleType b_t = BASE_A;
		double b_m = MASS_A;
		double b_r = A_R;
		double b_z = A_Z;
		double b_phi = A_PHI;
		switch (baseSequence[i]) {
		case 'A':
			b_t = complementarySequence ? BASE_T : BASE_A;
			b_m = MASS_A;
			b_r = A_R;
			b_z = A_Z;
			b_phi = A_PHI;
			break;
		case 'T':
			b_t = complementarySequence ? BASE_A : BASE_T;
			b_m = MASS_T;
			b_r = T_R;
			b_z = T_Z;
			b_phi = T_PHI;
			break;
		case 'C':
			b_t = complementarySequence ? BASE_G : BASE_C;
			b_m = MASS_C;
			b_r = C_R;
			b_z = C_Z;
			b_phi = C_PHI;
			break;
		case 'G':
			b_t = complementarySequence ? BASE_C : BASE_G;
			b_m = MASS_G;
			b_r = G_R;
			b_z = G_Z;
			b_phi = G_PHI;
			break;
		default:
			fprintf(stderr, "Unknown base type '%c' at "
					"position %d in base sequence '%s'! "
					"Defaulting to Adenine!\n",
					baseSequence[i], i, baseSequence);
		}

		/* Type */
		s->Bs[i].type = b_t;
		s->Ss[i].type = SUGAR;
		s->Ps[i].type = PHOSPHATE;

		/* Mass */
		s->Bs[i].m = b_m;
		s->Ss[i].m = MASS_S;
		s->Ps[i].m = MASS_P;

		/* Positions */
		s->Bs[i].pos = fromCilindrical(b_r, phi + b_phi, z + b_z);
		s->Ss[i].pos = fromCilindrical(S_R, phi + S_PHI, z + S_Z);
		s->Ps[i].pos = fromCilindrical(P_R, phi + P_PHI, z + P_Z);
		add(&s->Bs[i].pos, &offset, &s->Bs[i].pos);
		add(&s->Ss[i].pos, &offset, &s->Ss[i].pos);
		add(&s->Ps[i].pos, &offset, &s->Ps[i].pos);

		/* Velocity */
		s->Bs[i].vel = randNormVec(velStdev);
		s->Ss[i].vel = randNormVec(velStdev);
		s->Ps[i].vel = randNormVec(velStdev);

		/* Particle's strand */
		s->Bs[i].strand = s; s->Bs[i].strandIndex = i;
		s->Ss[i].strand = s; s->Ss[i].strandIndex = i;
		s->Ps[i].strand = s; s->Ps[i].strandIndex = i;

		z += HELIX_DELTA_Z;
		phi += HELIX_DELTA_PHI;
	}
}
void fillStrand(Strand *s, const char *baseSequence)
{
	fillStrandHelper(s, baseSequence, false, false);
}
void fillComplementaryStrand(Strand *s, const char *baseSequence)
{
	fillStrandHelper(s, baseSequence, true, true);
}

void freeStrand(Strand *strand)
{
	free(strand->all);
	strand->numMonomers = 0;
}

void freeWorld(void)
{
	assert(world.strands != NULL);
	for (int s = 0; s < world.numStrands; s++)
		freeStrand(&world.strands[s]);
	free(world.strands);
	freeGrid();
	return;
}



int numParticles(void)
{
	int num = 0;
	for (int s = 0; s < world.numStrands; s++)
		num += world.strands[s].numMonomers;
	return 3*num; /* three particles per monomer. */
}



/* ===== ITERATION FUNCTIONS ===== */

/* loop over every particle in the world */
void forEveryParticle(void (*f)(Particle *p))
{
	for (int s = 0; s < world.numStrands; s++)
		forEveryParticleOf(&world.strands[s], f);
}
/* loop over every particle in the world, pass [D]ata to the function */
void forEveryParticleD(void (*f)(Particle *p, void *data), void *data)
{
	for (int s = 0; s < world.numStrands; s++)
		forEveryParticleOfD(&world.strands[s], f, data);
}
/* loop over every particle in the strand */
void forEveryParticleOf(Strand *s, void (*f)(Particle *p))
{
	for (int i = 0; i < 3 * s->numMonomers; i++)
		f(&s->all[i]);
}
/* loop over every particle in the strand, pass [D]ata to the function */
void forEveryParticleOfD(Strand *s,
			void (*f)(Particle *p, void *data), void *data)
{
	for (int i = 0; i < 3 * s->numMonomers; i++)
		f(&s->all[i], data);
}

Particle *getConnectedParticle(Particle *p)
{
	Strand *s = p->strand;
	assert(s != NULL);
	int i = p->strandIndex;

	switch (p->type) {
	case PHOSPHATE:
		return &s->Ss[i];
	case SUGAR:
		if (i-1 < 0)
			return NULL;
		return &s->Ps[i-1];
	default: /* Base */
		return &s->Ss[i];
	}
}

