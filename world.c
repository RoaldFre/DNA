#include "world.h"
#include "spgrid.h"
#include "physics.h" /* for spacings/angles/... */
#include <string.h>

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
	int n = strlen(baseSequence);
	allocStrand(s, n);

	double ws = world.worldSize;
	double spacing = D_S5P + D_S3P; /* vertical spacing between monomers */
	double xoffset = -D_SA / 2;
	double zoffset = -D_SA / 2;
	double yoffset = (ws - n * spacing) / 2;
	double posStdev = spacing / 100;

	for (int i = 0; i < n; i++) {
		/* screw factor */
		double screwFactorCos = cos(i * SCREW_SYM_PHI);
		double screwFactorSin = sin(i * SCREW_SYM_PHI);

		if (complementaryHelix) {
			screwFactorCos *= -1;
			screwFactorSin *= -1;
		}
		
		/* Positions */
		s->Ss[i].pos.x = ws/2 + (xoffset - D_SA) * screwFactorCos;
		s->Bs[i].pos.x = ws/2 + (xoffset       ) * screwFactorCos;
		s->Ps[i].pos.x = ws/2 + (xoffset - D_SA) * screwFactorCos;

		s->Ss[i].pos.y = yoffset + i*spacing;
		s->Bs[i].pos.y = yoffset + i*spacing;
		s->Ps[i].pos.y = yoffset + i*spacing + D_S5P;
		
		
		s->Ss[i].pos.z = ws/2 + (zoffset - D_SA) * screwFactorSin;
		s->Bs[i].pos.z = ws/2 + (zoffset       ) * screwFactorSin;
		s->Ps[i].pos.z = ws/2 + (zoffset - D_SA) * screwFactorSin;

		s->Ss[i].pos.x += posStdev * randNorm();
		s->Ss[i].pos.y += posStdev * randNorm();
		s->Ss[i].pos.z += posStdev * randNorm();
		s->Bs[i].pos.x += posStdev * randNorm();
		s->Bs[i].pos.y += posStdev * randNorm();
		s->Bs[i].pos.z += posStdev * randNorm();
		s->Ps[i].pos.x += posStdev * randNorm();
		s->Ps[i].pos.y += posStdev * randNorm();
		s->Ps[i].pos.z += posStdev * randNorm();

		/* Velocity */
		s->Ss[i].vel = (Vec3) {0, 0, 0};
		s->Bs[i].vel = (Vec3) {0, 0, 0};
		s->Ps[i].vel = (Vec3) {0, 0, 0};

		/* Mass */
		s->Ss[i].m = MASS_S;
		s->Bs[i].m = MASS_A;
		s->Ps[i].m = MASS_P;

		/* Type */
		ParticleType b;
		switch (baseSequence[i]) {
		case 'A': b = complementarySequence ? BASE_T : BASE_A; break;
		case 'T': b = complementarySequence ? BASE_A : BASE_T; break;
		case 'C': b = complementarySequence ? BASE_G : BASE_C; break;
		case 'G': b = complementarySequence ? BASE_C : BASE_G; break;
		default: 
			fprintf(stderr, "Unknown base type '%c' at "
				"position %d in base sequence %s. "
				"Defaulting to Adenine!\n",
				baseSequence[i], i, baseSequence);
			b = BASE_A;
		}

		s->Bs[i].type = b;
		s->Ss[i].type = SUGAR;
		s->Ps[i].type = PHOSPHATE;

		/* Particle's strand */
		s->Ss[i].strand = s; s->Ss[i].strandIndex = i;
		s->Bs[i].strand = s; s->Bs[i].strandIndex = i;
		s->Ps[i].strand = s; s->Ps[i].strandIndex = i;
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
	int n = s->numMonomers;

	switch (p->type) {
	case PHOSPHATE:
		return &s->Ss[i];
	case SUGAR:
		if (i+1 >= n)
			return NULL;
		return &s->Ps[i+1];
	default: /* Base */
		return &s->Ss[i];
	}
}

