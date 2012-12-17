#include "world.h"
#include "spgrid.h"
#include "physics.h"
#include "integrator.h"
#include <string.h>
#include <stdio.h>


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
 * complementaryHelix:    Build the helix 'upside down' and turning the other way around.
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

/* ********************
 * BIG GIANT MAJOR TODO:
 * ********************
 * COMPLEMENTARY STRAND DOES NOT LIKE THE WAY WE NUMBER MONOMERS, OR 
 * SOMETHING LIKE THAT. THIS IS APPARANT IN, FOR INSTANCE THE EXCLUSION 
 * INTERACTION, WHERE PHOSPHATES AND MONOMERS OF THE COMPLEMENTARY STRAND 
 * REPULSE EACH OTHER! */
static void fillStrandHelper(Strand *s, const char *baseSequence,
		bool complementarySequence, bool complementaryHelix)
{
	/* right order when we are building the complementary strand */
	double order;
	if (complementarySequence)
		order = -1;
	else
		order = 1;


	int n = strlen(baseSequence);
	allocStrand(s, n);

	Vec3 offset;
	offset.x = 0;
	offset.z = 0;
	offset.y = -n * HELIX_DELTA_Z / 2.0;

	/* <v^2> = 3 T k_B / m
	 * -> per dimension: gaussian with variance T k_B / m */
	double T = getHeatBathTemperature();
	double velVarPerInvMass = T * BOLTZMANN_CONSTANT;
	double dt = getTimeStep();

	double phi = 0;
	double z = 0;
	int j;
	
	if (complementarySequence) {
		z = (n - 1) * HELIX_DELTA_Z;
		phi = (n - 1) * HELIX_DELTA_PHI;
	}
	
	for (int i = 0; i < n; i++) {
		/* Default to Adenine */
		ParticleType b_t = BASE_A;
		double b_m = A_M;
		double b_r = A_R;
		double b_z = A_Z;
		double b_phi = A_PHI;
		j = i;
		if (complementarySequence)
			j = n - 1 - i;
		
		switch (baseSequence[j]) {
		case 'A':
			b_t = complementarySequence ? BASE_T : BASE_A;
			b_m = complementarySequence ? T_M : A_M;
			b_r = complementarySequence ? T_R : A_R;
			b_z = complementarySequence ? T_Z : A_Z;
			b_phi = complementarySequence ? T_PHI : A_PHI;
			break;
		case 'T':
			b_t = complementarySequence ? BASE_A : BASE_T;
			b_m = complementarySequence ? A_M : T_M;
			b_r = complementarySequence ? A_R : T_R;
			b_z = complementarySequence ? A_R : T_Z;
			b_phi = complementarySequence ? A_PHI : T_PHI;
			break;
		case 'C':
			b_t = complementarySequence ? BASE_G : BASE_C;
			b_m = complementarySequence ? G_M : C_M;
			b_r = complementarySequence ? G_R : C_R;
			b_z = complementarySequence ? G_Z : C_Z;
			b_phi = complementarySequence ? G_PHI : C_PHI;
			break;
		case 'G':
			b_t = complementarySequence ? BASE_C : BASE_G;
			b_m = complementarySequence ? C_M : G_M;
			b_r = complementarySequence ? C_R : G_R;
			b_z = complementarySequence ? C_Z : G_Z;
			b_phi = complementarySequence ? C_PHI : G_PHI;
			break;
		case 'X':
			b_t = complementarySequence ? BASE_C : BASE_X;
			b_m = complementarySequence ? C_M : X_M;
			b_r = complementarySequence ? C_R : X_R;
			b_z = complementarySequence ? C_Z : X_Z;
			b_phi = complementarySequence ? C_PHI : X_PHI;
			break;
		case 'Y':
			b_t = complementarySequence ? BASE_C : BASE_Y;
			b_m = complementarySequence ? C_M : Y_M;
			b_r = complementarySequence ? C_R : Y_R;
			b_z = complementarySequence ? C_Z : Y_Z;
			b_phi = complementarySequence ? C_PHI : Y_PHI;
			break;
		default:
			fprintf(stderr, "Unknown base type '%c' at "
					"position %d in base sequence '%s'! "
					"Defaulting to Adenine!\n",
					baseSequence[j], j, baseSequence);
		}

		/* Type */
		s->Bs[i].type = b_t;
		s->Ss[i].type = SUGAR;
		s->Ps[i].type = PHOSPHATE;

		/* Mass */
		s->Bs[i].m = b_m;
		s->Ss[i].m = S_M;
		s->Ps[i].m = P_M;

		/* Positions */
		if (!complementaryHelix) {
			s->Bs[i].pos = fromCilindrical(b_r, phi + b_phi, z + b_z);
			s->Ss[i].pos = fromCilindrical(S_R, phi + S_PHI, z + S_Z);
			s->Ps[i].pos = fromCilindrical(P_R, phi + P_PHI, z + P_Z);
		} else {
		/* Screw assymetry */	
			s->Bs[i].pos = fromCilindrical(b_r, (phi - b_phi), z - b_z);
			s->Ss[i].pos = fromCilindrical(S_R, (phi - S_PHI), z - S_Z);
			s->Ps[i].pos = fromCilindrical(P_R, (phi - P_PHI), z - P_Z);
		}
		
		
		s->Bs[i].pos = add(s->Bs[i].pos, offset);
		s->Ss[i].pos = add(s->Ss[i].pos, offset);
		s->Ps[i].pos = add(s->Ps[i].pos, offset);

		/* Velocity */
		s->Bs[i].vel = randNormVec(sqrt(velVarPerInvMass / s->Bs[i].m));
		s->Ss[i].vel = randNormVec(sqrt(velVarPerInvMass / s->Ss[i].m));
		s->Ps[i].vel = randNormVec(sqrt(velVarPerInvMass / s->Ps[i].m));

		/* Previous position */
		s->Bs[i].prevPos = sub(s->Bs[i].pos, scale(s->Bs[i].vel, dt));
		s->Ss[i].prevPos = sub(s->Ss[i].pos, scale(s->Ss[i].vel, dt));
		s->Ps[i].prevPos = sub(s->Ps[i].pos, scale(s->Ps[i].vel, dt));

		/* Particle's strand */
		s->Bs[i].strand = s; s->Bs[i].strandIndex = j;
		s->Ss[i].strand = s; s->Ss[i].strandIndex = j;
		s->Ps[i].strand = s; s->Ps[i].strandIndex = j;

		z += order * HELIX_DELTA_Z;
		phi += order * HELIX_DELTA_PHI;
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

/* Returns the base sequence string of the given strand. You need to free 
 * the returned pointer yourself! */
char *getSequence(Strand *s)
{
	int n = s->numMonomers;
	char *seq = calloc(n + 1, sizeof(*seq));

	for (int i = 0; i < n; i++) {
		switch (s->Bs[i].type) {
		case BASE_A: seq[i] = 'A'; break;
		case BASE_T: seq[i] = 'T'; break;
		case BASE_C: seq[i] = 'C'; break;
		case BASE_G: seq[i] = 'G'; break;
		case BASE_X: seq[i] = 'X'; break;
		case BASE_Y: seq[i] = 'Y'; break;
		default:
			assert(false);
			die("getSequence: got invalid strand!\n");
		}
	}

	seq[n] = '\0';
	return seq;
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
		return ((i-1 < 0) ? NULL : &s->Ps[i-1]);
	default: /* Base */
		return &s->Ss[i];
	}
}




/* ===== MISC FUNCTIONS ===== */

static void translateParticle(Particle *p, void *data)
{
	Vec3 *delta = (Vec3*) data;
	p->pos = add(p->pos, *delta);
}
void translateStrand(Strand *s, Vec3 delta)
{
	forEveryParticleOfD(s, translateParticle, (void*) &delta);
}
char *getWorldInfo(void)
{
	int n = world.numStrands;
	char *ret = asprintfOrDie("# Number of strands in the world: %d\n", n);

	for (int i = 0; i < n; i++) {
		/* This isn't very malloc friendly and quadratic in the 
		 * string length, but meh. */
		char *tmp = ret;
		char *seq = getSequence(&world.strands[i]);
		ret = asprintfOrDie("%s# Strand %d: %s\n", tmp, i + 1, seq);
		free(tmp);
		free(seq);
	}

	return ret;
}





/* ===== CHECK FUNCTIONS ===== */
bool strandSanityCheck(Strand *s)
{
	assert(s != NULL);
	bool OK = true;
	for (int i = 0; i < s->numMonomers; i++) {
		if (!isBase(s->Bs[i].type)) {
			fprintf(stderr, "Particle %d of strand %p should be "
					"a base but isn't!\n", i, (void*)s);
			OK = false;
		}
		if (s->Ps[i].type != PHOSPHATE) {
			fprintf(stderr, "Particle %d of strand %p should be a "
					"phosphate but isn't!\n", i, (void*)s);
			OK = false;
		}
		if (s->Ss[i].type != SUGAR) {
			fprintf(stderr, "Particle %d of strand %p should be "
					"a sugar but isn't!\n", i, (void*)s);
			OK = false;
		}
		if (s->Bs[i].strand != s || s->Ps[i].strand != s
				|| s->Ss[i].strand != s) {
			fprintf(stderr, "Monomer %d of strand %p has "
					"particles associated with the wrong "
					"strand!\n", i, (void*)s);
			OK = false;
		}
		//XXX disable this because we break this above when filling 
		//the strand (for now!...)
#if 0
		if (s->Bs[i].strandIndex != i || s->Ps[i].strandIndex != i
				|| s->Ss[i].strandIndex != i) {
			fprintf(stderr, "Monomer %d of strand %p has "
					"particles with a wrong "
					"strand index!\n", i, (void*)s);
			OK = false;
		}
#endif
	}

	return OK;
}

bool worldSanityCheck(void)
{
	bool OK = true;
	for (int i = 0; i < world.numStrands; i++) {
		if (strandSanityCheck(&world.strands[i]))
			continue;
		fprintf(stderr, "-> Strand %d had errors!\n\n", i);
		OK = false;
	}

	return OK;
}
