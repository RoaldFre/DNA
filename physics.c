#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "physics.h"
#include "vmath.h"
#include "spgrid.h"

/* Masses (in kg) */
#define AU      1.660539e-27
#define MASS_A  (134.1 * AU)
#define MASS_P  (94.97 * AU)
#define MASS_S  (83.11 * AU)

/* Equilibrium distance of bonds (in m) */
#define D_S5P   3.899e-10
#define D_S3P   3.559e-10
#define D_SA    6.430e-10

/* Equilibrium distance of stacking potential (in m) */
#define STACK_SIGMA  (3.414e-10)

/* Energy unit */
#define EPSILON 1.81e-21 /* 0.26kcal/mol == 1.81 * 10^-21 J (per particle) */


/* Bond stretch */
#define BOND_K1      10*(EPSILON * FROM_ANGSTROM_SQUARED)
#define BOND_K2      (100 * EPSILON * FROM_ANGSTROM_SQUARED)
/* Bond bend */
#define BOND_Ktheta  (400 * EPSILON) /* per radian^2 */
/* Bond twist */
#define BOND_Kphi    (4 * EPSILON)
/* Bond stack */
#define BOND_STACK   EPSILON

/* Screw symmetry constants */
#define SCREW_SYM_PHI	(36 * TO_RADIANS)
// #define SCREW_SYM_Z		(3.38e-10) /* 3.38 angstrom */

/* Bond angle */
#define ANGLE_S5_P_3S	( 94.49 * TO_RADIANS)
#define ANGLE_P_5S3_P	(120.15 * TO_RADIANS)
#define ANGLE_P_5S_A	(113.13 * TO_RADIANS)
#define ANGLE_P_3S_A	(108.38 * TO_RADIANS)

/* Dihedral angle */
#define DIHEDRAL_P_5S3_P_5S	(-154.80 * TO_RADIANS)
#define DIHEDRAL_S3_P_5S3_P	(-179.17 * TO_RADIANS)
#define DIHEDRAL_A_S3_P_5S	( -22.60 * TO_RADIANS)
#define DIHEDRAL_S3_P_5S_A	(  50.69 * TO_RADIANS)

/* Base-Pair couplings */
#define COUPLING_BP_AT 	19.28e-11  /* 2.77 kcal/mol (per particle) */
#define COUPLING_BP_GC	28.96e-11  /* 4.16 kcal/mol (per particle) */
#define DISTANCE_r0_AT	2.9002e-10 /* Knotts et al 2007, table III, 2.9002 A */
#define DISTANCE_r0_GC	2.8694e-10 /* Knotts et al 2007, table III, 2.8694 A */

/* Coulomb interaction between phosphates */
#define CHARGE_ELECTRON		1.602e-19  /* 1.602 Coulomb */
#define VACUUM_PERMITTIVITY	8.8541e-12 /* 8.854e-12 Farads/m */
#define COUPLING_EPS_H2O	(78*VACUUM_PERMITTIVITY) /* 78 epsilon_0 */
#define DEBYE_LENGTH 		13.603e-10 /* 13.603 Angstrom for 50mM = [Na+]*/

/* Thermodynamics */
#define ENERGY_FACTOR	      (1/1.602177e-19) /* Energy in electronvolt */
#define BOLTZMANN_CONSTANT    1.38065e-23
#define FROM_ANGSTROM_SQUARED 1e20 /* Bond constants are given for Angstrom */
#define TO_RADIANS	      (M_PI / 180)

#define DIELECTRIC_CST_H20    80
#define AVOGADRO              6.023e23 /* particles per mol */




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
	config.worldSize = worldSize; //TODO better place for this?
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


/* Returns a number sampled from a standard normal distribution. */
static double randNorm(void)
{
	/* Box-Muller transform */
	double u1 = ((double) rand()) / RAND_MAX;
	double u2 = ((double) rand()) / RAND_MAX;

	return sqrt(-2*log(u1)) * cos(2*M_PI*u2);
}

/* Returns a vector with components sampled from a standard normal 
 * distribution. */
static Vec3 randNormVec(void)
{
	Vec3 res;
	res.x = randNorm();
	res.y = randNorm();
	res.z = randNorm();
	return res;
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

	double ws = config.worldSize;
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

static int numParticles(void)
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



/* ===== FORCES AND POTENTIALS ===== */

/* V = k1 * (dr - d0)^2  +  k2 * (d - d0)^4
 * where dr is the distance between the particles */
static double Vbond(Particle *p1, Particle *p2, double d0)
{
	double k1 = BOND_K1;
	double k2 = BOND_K2;
	double d = nearestImageDistance(&p1->pos, &p2->pos) - d0;
	double d2 = d * d;
	double d4 = d2 * d2;
	return k1 * d2  +  k2 * d4;
}
static void Fbond(Particle *p1, Particle *p2, double d0)
{
	double k1 = BOND_K1;
	double k2 = BOND_K2;
	Vec3 drVec, drVecNormalized, F;
	drVec = nearestImageVector(&p1->pos, &p2->pos);
	double dr = length(&drVec);
	double d  = dr - d0;
	double d3 = d * d * d;

	scale(&drVec, 1/dr, &drVecNormalized);
	scale(&drVecNormalized, 2*k1*d + 4*k2*d3, &F);

	add(&p1->F, &F, &p1->F);
	sub(&p2->F, &F, &p2->F);
}

/* V = ktheta * (theta - theta0) 
 *
 * p1 \       /p3
 *     \theta/
 *      \   /
 *       \ /
 *        p2
 */
static double Vangle(Particle *p1, Particle *p2, Particle *p3, double theta0)
{
	Vec3 a, b;
	double ktheta = BOND_Ktheta;
	
	a = nearestImageVector(&p2->pos, &p1->pos);
	b = nearestImageVector(&p2->pos, &p3->pos);

	double dtheta = angle(&a, &b) - theta0;
	return ktheta/2 * dtheta*dtheta;
}
static void Fangle(Particle *p1, Particle *p2, Particle *p3, double theta0)
{
	Vec3 a, b;
	double ktheta = BOND_Ktheta;
	
	a = nearestImageVector(&p2->pos, &p1->pos);
	b = nearestImageVector(&p2->pos, &p3->pos);
	
	double lal = length(&a);
	double lbl = length(&b);
	double adotb = dot(&a, &b);
	double costheta = adotb / (lal * lbl);
	double theta = acos(costheta);
	double sintheta = sqrt(1 - costheta*costheta);

	// TODO correct cut off?
	if (fabs(sintheta) < 1e-30)
		/* "No" force (unstable equilibrium), numerical instability 
		 * otherwise */
		return;

	Vec3 tmp1, tmp2, F1, F2, F3;

	scale(&b, 1/(lal * lbl), &tmp1);
	scale(&a, adotb / (lal*lal*lal * lbl), &tmp2);
	sub(&tmp1, &tmp2, &F1);
	scale(&F1, ktheta * (theta - theta0) / sintheta, &F1);
	add(&p1->F, &F1, &p1->F);	

	scale(&a, 1/(lal * lbl), &tmp1);
	scale(&b, adotb / (lbl*lbl*lbl * lal), &tmp2);
	sub(&tmp1, &tmp2, &F3);
	scale(&F3, ktheta * (theta - theta0) / sintheta, &F3);
	add(&p3->F, &F3, &p3->F);	

	add(&F1, &F3, &F2);
	sub(&p2->F, &F2, &p2->F);	

	assert(ktheta == 0 
		|| fabs(dot(&a, &F1) / length(&a) / length(&F1)) < 1e-5);
	assert(ktheta == 0
		|| fabs(dot(&b, &F3) / length(&b) / length(&F3)) < 1e-5);
}

static double Vdihedral(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
								double phi0)
{
	double kphi = BOND_Kphi;
	Vec3 r1, r2, r3;
	sub(&p2->pos, &p1->pos, &r1);
	sub(&p3->pos, &p2->pos, &r2);
	sub(&p4->pos, &p3->pos, &r3);
	
	double phi = dihedral(&r1, &r2, &r3);
	//printf("phi = %f\n",phi / TO_RADIANS);
	return kphi * (1 - cos(phi - phi0));
}
static void FdihedralParticle(Particle *target, 
		Particle *p1, Particle *p2, Particle *p3, Particle *p4, 
		double Vorig, double phi0)
{
	double hfactor = 1e-8; /* roughly sqrt(epsilon) for a double */
	double h;
	Vec3 F;

	h = target->pos.x * hfactor;
	target->pos.x += h;
	F.x = (Vorig - Vdihedral(p1, p2, p3, p4, phi0)) / h;
	target->pos.x -= h;

	h = target->pos.y * hfactor;
	target->pos.y += h;
	F.y = (Vorig - Vdihedral(p1, p2, p3, p4, phi0)) / h;
	target->pos.y -= h;

	h = target->pos.z * hfactor;
	target->pos.z += h;
	F.z = (Vorig - Vdihedral(p1, p2, p3, p4, phi0)) / h;
	target->pos.z -= h;

	add(&target->F, &F, &target->F);
}
static void Fdihedral(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
								double phi0)
{
	/* This is a *mess* to do analytically, so we do a numerical 
	 * differentiation instead. */
	double Vorig = Vdihedral(p1, p2, p3, p4, phi0);
	FdihedralParticle(p1, p1, p2, p3, p4, Vorig, phi0);
	FdihedralParticle(p2, p1, p2, p3, p4, Vorig, phi0);
	FdihedralParticle(p3, p1, p2, p3, p4, Vorig, phi0);
	FdihedralParticle(p4, p1, p2, p3, p4, Vorig, phi0);
}

static double Vstack(Particle *p1, Particle *p2)
{
	double kStack = BOND_STACK;
	double sigma = STACK_SIGMA;
	double sigma2 = sigma * sigma;
	double sigma6 = sigma2 * sigma2 * sigma2;
	double sigma12 = sigma6 * sigma6;
	double r2 = distance2(&p1->pos, &p2->pos);
	double r6 = r2 * r2 * r2;
	double r12 = r6 * r6;

	return kStack * (sigma12/r12 - 2*sigma6/r6 + 1);
	//return kStack * (sigma12/r12 - 2*sigma6/r6);
}
static void Fstack(Particle *p1, Particle *p2)
{
	double kStack = BOND_STACK;
	double sigma = STACK_SIGMA;
	double sigma2 = sigma * sigma;
	double sigma6 = sigma2 * sigma2 * sigma2;
	double sigma12 = sigma6 * sigma6;

	Vec3 Fi;
	Vec3 drVec;
	sub(&p2->pos, &p1->pos, &drVec);
	double dr = length(&drVec);

	assert(dr != 0);

	double dr2 = dr*dr;
	double dr3 = dr*dr*dr;
	double dr6 = dr3*dr3;
	double dr8 = dr6*dr2;
	double dr12 = dr6*dr6;
	double dr14 = dr12*dr2;

	scale(&drVec, -12 * kStack * (sigma12/dr14 - sigma6/dr8), &Fi);
	add(&p1->F, &Fi, &p1->F);
	sub(&p2->F, &Fi, &p2->F);
}


static double calcVBasePai(double coupling, double rij0, double rijVar2)
{
	double rfrac2 = rij0 * rij0 / rijVar2;
	double rfrac4 = rfrac2 * rfrac2;
	double rfrac8 = rfrac4 * rfrac4;
	double rfrac10 = rfrac8 * rfrac2;
	double rfrac12 = rfrac8 * rfrac4;
	
	return coupling*(5*rfrac12 - 6*rfrac10 + 1);
}
static double VbasePair(Particle *p1, Particle *p2)
{
	double rij = nearestImageDistance(&p1->pos, &p2->pos);
	if (rij > config.truncationLen)
		return 0; /* Too far away */

	double rij2 = rij*rij;
	
	double bpCoupling;
	double bpForceDist;
	double bpPotential;
	double truncCorrection = 0;
	
	/* Apply right potential constant for AT-bonding and GC-bonding, 
	 * if not AT or GC then zero */
	
	/* Lennard-Jones potential:
	 * potential = bpCoupling * 5*(r0 / r)^12 - 6*(r0/r)^10 + 1. */
	
	if ((p1->type==BASE_A && p2->type==BASE_T) 
			|| (p1->type==BASE_T && p2->type==BASE_A)) {
		bpCoupling = COUPLING_BP_AT;
		bpForceDist = DISTANCE_r0_AT;
		/* calculate the correction by which the force should be lifted */
		double truncLen2 = config.truncationLen * config.truncationLen;
		truncCorrection = calcVBasePai(bpCoupling, bpForceDist, truncLen2);
		bpPotential = calcVBasePai(bpCoupling, bpForceDist, rij2) - truncCorrection;
		
	} else if ((p1->type==BASE_G && p2->type==BASE_C) 
			|| (p1->type==BASE_C && p2->type==BASE_G)) {
		bpCoupling = COUPLING_BP_GC;
		bpForceDist = DISTANCE_r0_GC;
		/* calculate the correction by which the force should be lifted */
		double truncLen2 = config.truncationLen * config.truncationLen;
		truncCorrection = calcVBasePai(bpCoupling, bpForceDist, truncLen2);
		bpPotential = calcVBasePai(bpCoupling, bpForceDist, rij2) - truncCorrection;
		
	/* Else, no L-J potential and return potential = 0 */
	} else {
		return 0; 
	}

	return bpPotential;
}

static double calcFbasePair(double coupling, double rij0, double rijVar)
{
	double rfrac = rij0 / rijVar;
	double rfrac2 = rfrac * rfrac;
	double rfrac4 = rfrac2 * rfrac2;
	double rfrac8 = rfrac4 * rfrac4;
	double rfrac10 = rfrac8 * rfrac2;
	double rfrac12 = rfrac10 * rfrac2;
				
	return coupling*60*( rfrac12 / rijVar - rfrac10 / rijVar );
}
static void FbasePair(Particle *p1, Particle *p2)
{
	double bpCoupling;
	double bpForceDist;
	double truncLen = config.truncationLen;
	Vec3 forceVec;
	double force;
	
	double rij = nearestImageDistance(&p1->pos, &p2->pos);
	if (rij > truncLen)
		return; /* Too far away */

	Vec3 direction = nearestImageUnitVector(&p1->pos, &p2->pos);
	
	/* Apply right force constant for AT-bonding and GC-bonding, 
	 * if not AT or GC then zero */
	
	/* Lennard-Jones potential:
	 * potential = bpCoupling * 5*(r0 / r)^12 - 6*(r0/r)^10 + 1.
	 * 
	 * For the force we differentiate with respect to r, so we get
	 * bpCoupling * 60*[ r0^10 / r^11 - r0^12/r^13 ] */

	if ((p1->type==BASE_A && p2->type==BASE_T) 
			|| (p1->type==BASE_T && p2->type==BASE_A)) {
		bpCoupling = COUPLING_BP_AT;
		bpForceDist = DISTANCE_r0_AT;
		force = calcFbasePair(bpCoupling, bpForceDist, rij);
		
	} else if ((p1->type==BASE_G && p2->type==BASE_C) 
			|| (p1->type==BASE_C && p2->type==BASE_G)) {
		bpCoupling = COUPLING_BP_GC;
		bpForceDist = DISTANCE_r0_GC;
		force = calcFbasePair(bpCoupling, bpForceDist, rij);
	
	/* Else, no force and return without modifying particles */
	} else {
		return; 
	}
	/* Scale the direction with the calculated force */
	scale(&direction, force, &forceVec);

	/* Add force to particle objects */
	add(&p1->F, &forceVec, &p1->F);
	sub(&p2->F, &forceVec, &p2->F);
}


static double calcDebyeLength(void)
{
	double T = config.thermostatTemp;
	double saltCon = (double)(config.saltConcentration);
	double elCharge = CHARGE_ELECTRON;
	double lambdaBDenom, lambdaB, kInvSq, kInv, kDebye;
	
	lambdaBDenom = 4 * M_PI * VACUUM_PERMITTIVITY * DIELECTRIC_CST_H20 *
			BOLTZMANN_CONSTANT * T;
	lambdaB = elCharge*elCharge / lambdaBDenom;
	kInvSq = 8*M_PI*lambdaB*AVOGADRO*saltCon;
	kInv = sqrt(kInvSq);
	kDebye = 1/kInv;
	
	return kDebye;
}
static double calcVCoulomb(double distanceLength)
{
	double couplingConstant = CHARGE_ELECTRON*CHARGE_ELECTRON/
				(4*M_PI*COUPLING_EPS_H2O);
	double debyeLength = calcDebyeLength();
	
	double expArgument = - distanceLength/debyeLength;
	double exponentialPart = exp(expArgument);
	
	double potentialQQ = couplingConstant * exponentialPart / distanceLength;
	
	return potentialQQ*exponentialPart;	
}
static double VCoulomb(Particle *p1, Particle *p2)
{
	double rij = nearestImageDistance(&p1->pos, &p2->pos);
	if (rij > config.truncationLen)
		return 0; /* Too far away */
	
	double truncLength = config.truncationLen;	
	double phPotential;
	double truncCorrection = 0;
	
	if (p1->type==PHOSPHATE && p2->type==PHOSPHATE) {
		/* calculate the correction by which the force should be lifted */
		truncCorrection = calcVCoulomb(truncLength);
		phPotential = calcVCoulomb(rij) - truncCorrection;
	
	/* Else, no Coulomb potential and return potential = 0 */
	} else {
		return 0; 
	}
	
	return phPotential;
}


static double calcFCoulomb(double distanceLength)
{
	double couplingConstant = CHARGE_ELECTRON*CHARGE_ELECTRON/
					(4*M_PI*COUPLING_EPS_H2O);
	double debyeLength = calcDebyeLength();
		
	double expArgument = - distanceLength/debyeLength;
	double exponentialPart = exp(expArgument);
	
	double forceQQ = couplingConstant* exponentialPart*
				( 1/(distanceLength*distanceLength) 
				+ 1/(debyeLength*distanceLength) );
	
	return forceQQ;
	
}
static void FCoulomb(Particle *p1, Particle *p2)
{	
	double truncLen = config.truncationLen;
	double rij = nearestImageDistance(&p2->pos, &p1->pos);
	if (rij > truncLen)
		return; /* Too far away */
	
	/* If phosphates, apply Coulomb interaction between them */

	double force;

	if (p1->type==PHOSPHATE && p2->type==PHOSPHATE) {
		force = calcFCoulomb(rij);
	} else {
		return;
	}
	
	Vec3 direction = nearestImageUnitVector(&p1->pos, &p2->pos);
	Vec3 forceVec;
	
	/* Scale the direction with the calculated force */
	scale(&direction, force, &forceVec);

	/* Add force to particle objects */
	add(&p1->F, &forceVec, &p1->F);
	sub(&p2->F, &forceVec, &p2->F);
}





/* ===== FORCE FUNCTIONS ===== */

static void resetForce(Particle *p)
{
	p->F.x = p->F.y = p->F.z = 0;
}

/* Calculate the forces on the particles of the strand that are attributed 
 * to the structure of the strand itself. These are:
 * Fbond, Fangle, Fdihedral, Fstack. */
static void strandForces(Strand *s) {
	if (s->numMonomers < 0)
		return;
	assert(s->Ss != NULL && s->Ps != NULL && s->Bs != NULL);

	/* Bottom monomer */
	Fbond(&s->Ss[0], &s->Bs[0], D_SA);
	Fbond(&s->Ss[0], &s->Ps[0], D_S5P);
	Fangle(&s->Ps[0], &s->Ss[0], &s->Bs[0], ANGLE_P_5S_A);
	/* Rest of the monomers */
	for (int i = 1; i < s->numMonomers; i++) {
		Fbond(&s->Ss[i], &s->Bs[i],   D_SA);
		Fbond(&s->Ss[i], &s->Ps[i],   D_S5P);
		Fbond(&s->Ss[i], &s->Ps[i-1], D_S3P);

		Fstack(&s->Bs[i], &s->Bs[i-1]);

		Fangle(&s->Ps[ i ], &s->Ss[ i ], &s->Bs[ i ], ANGLE_P_5S_A);
		Fangle(&s->Ps[ i ], &s->Ss[ i ], &s->Ps[i-1], ANGLE_P_5S3_P);
		Fangle(&s->Ps[i-1], &s->Ss[ i ], &s->Bs[ i ], ANGLE_P_3S_A);
		Fangle(&s->Ss[i-1], &s->Ps[i-1], &s->Ss[ i ], ANGLE_S5_P_3S);

		Fdihedral(&s->Ps[i], &s->Ss[ i ], &s->Ps[i-1], &s->Ss[i-1],
							DIHEDRAL_P_5S3_P_5S);
		Fdihedral(&s->Bs[i], &s->Ss[ i ], &s->Ps[i-1], &s->Ss[i-1],
							DIHEDRAL_A_S3_P_5S);
		Fdihedral(&s->Ss[i], &s->Ps[i-1], &s->Ss[i-1], &s->Bs[i-1],
							DIHEDRAL_S3_P_5S_A);
		if (i >= 2)
		Fdihedral(&s->Ss[i], &s->Ps[i-1], &s->Ss[i-1], &s->Ps[i-2],
							DIHEDRAL_S3_P_5S3_P);
	}
}

static void pairForces(Particle *p1, Particle *p2)
{
	FbasePair(p1, p2);
	FCoulomb(p1, p2);
}

static void calculateForces(void)
{
	/* Reset forces */
	forEveryParticle(&resetForce);

	/* Strand-based forces */
	for (int s = 0; s < world.numStrands; s++)
		strandForces(&world.strands[s]);

	/* Particle-based forces */
	forEveryPair(&pairForces);
}



/* ===== ENERGY FUNCTIONS ===== */

static void kineticHelper(Particle *p, void *data)
{
	double *twiceK = (double*) data;
	*twiceK += p->m * length2(&p->vel);
}
static double kineticEnergy(void)
{
	double twiceK = 0;
	forEveryParticleD(&kineticHelper, (void*) &twiceK);
	return twiceK/2;
}
double temperature(void)
{
	return 2.0 / (3.0 * BOLTZMANN_CONSTANT)
			* kineticEnergy() / numParticles();
}


typedef struct PotentialEnergies {
	double bond, angle, dihedral, stack, basePair, Coulomb;
} PotentialEnergies;
static void pairPotentials(Particle *p1, Particle *p2, void *data)
{
	PotentialEnergies *pe = (PotentialEnergies*) data;
	
	pe->basePair += VbasePair(p1, p2);
	pe->Coulomb  += VCoulomb(p1, p2);
}

/* Add energy stats of given strand, in electronvolts. */
static void addPotentialEnergies(Strand *s, PotentialEnergies *pe)
{
	double Vb = 0;
	double Va = 0;
	double Vd = 0;
	double Vs = 0;
	Vb += Vbond(&s->Ss[0], &s->Bs[0],   D_SA);
	Vb += Vbond(&s->Ss[0], &s->Ps[0],   D_S5P);

	Va += Vangle(&s->Bs[0], &s->Ss[0], &s->Ps[0], ANGLE_P_5S_A);
	for (int i = 1; i < s->numMonomers; i++) {
		Vb += Vbond(&s->Ss[i], &s->Bs[i],   D_SA);
		Vb += Vbond(&s->Ss[i], &s->Ps[i],   D_S5P);
		Vb += Vbond(&s->Ss[i], &s->Ps[i-1], D_S3P);

		Vs += Vstack(&s->Bs[i], &s->Bs[i-1]);

		Va += Vangle(&s->Ps[ i ], &s->Ss[ i ], &s->Bs[ i ], ANGLE_P_5S_A);
		Va += Vangle(&s->Ps[ i ], &s->Ss[ i ], &s->Ps[i-1], ANGLE_P_5S3_P);
		Va += Vangle(&s->Ps[i-1], &s->Ss[ i ], &s->Bs[ i ], ANGLE_P_3S_A);
		Va += Vangle(&s->Ss[i-1], &s->Ps[i-1], &s->Ss[ i ], ANGLE_S5_P_3S);

		Vd += Vdihedral(&s->Ps[i], &s->Ss[ i ], &s->Ps[i-1], &s->Ss[i-1],
							DIHEDRAL_P_5S3_P_5S);
		Vd += Vdihedral(&s->Bs[i], &s->Ss[ i ], &s->Ps[i-1], &s->Ss[i-1],
							DIHEDRAL_A_S3_P_5S);
		Vd += Vdihedral(&s->Ss[i], &s->Ps[i-1], &s->Ss[i-1], &s->Bs[i-1],
							DIHEDRAL_S3_P_5S_A);
		if (i >= 2)
		Vd += Vdihedral(&s->Ss[i], &s->Ps[i-1], &s->Ss[i-1], &s->Ps[i-2],
							DIHEDRAL_S3_P_5S3_P);
	}

	pe->bond     = Vb;
	pe->angle    = Va;
	pe->dihedral = Vd;
	pe->stack    = Vs;
}

/* Return energy stats of world, in electronvolts. */
static PotentialEnergies calcPotentialEnergies(void) {
	
	PotentialEnergies pe = {0, 0, 0, 0, 0, 0};
	for (int s = 0; s < world.numStrands; s++)
		addPotentialEnergies(&world.strands[s], &pe);
	
	forEveryPairD(&pairPotentials, &pe);
	
	/* Convert to eV */
	pe.bond     *= ENERGY_FACTOR;
	pe.angle    *= ENERGY_FACTOR;
	pe.dihedral *= ENERGY_FACTOR;
	pe.stack    *= ENERGY_FACTOR;
	pe.basePair *= ENERGY_FACTOR;
	pe.Coulomb  *= ENERGY_FACTOR;
	
	return pe;
}




/* ===== MOMENTUM FUNCTIONS ===== */

static void momentumHelper(Particle *p, void *data)
{
	Vec3 *Ptot = (Vec3*) data;
	Vec3 P;
	scale(&p->vel, p->m, &P);
	add(&P, Ptot, Ptot);
}
static Vec3 momentum(void)
{
	Vec3 Ptot = {0, 0, 0};
	forEveryParticleD(&momentumHelper, (void*) &Ptot);
	return Ptot;
}

bool physicsCheck(void)
{
	Vec3 P = momentum();
	double PPM = length(&P) / numParticles();
	if (PPM > 1e-20) {
		fprintf(stderr, "\nMOMENTUM CONSERVATION VIOLATED! "
				"Momentum per monomer: |P| = %e\n", PPM);
		return false;
	}
	return true;
}




/* ===== INTEGRATOR ===== */

static void thermostatHelper(Particle *p, void *data)
{
	double lambda = *(double*) data;
	scale(&p->vel, lambda, &p->vel);
}
static void thermostat(void)
{
	if (config.thermostatTau <= 0)
		return;

	/* Mass and Boltzmann constant are 1 */ 
	double Tk  = temperature();
	double T0  = config.thermostatTemp;
	double dt  = config.timeStep;
	double tau = config.thermostatTau;
	double lambda = sqrt(1 + dt/tau * (T0/Tk - 1));

	forEveryParticleD(&thermostatHelper, (void*) &lambda);
}

static void verletHelper1(Particle *p)
{
	double dt = config.timeStep;
	Vec3 tmp;

	/* vel(t + dt/2) = vel(t) + acc(t)*dt/2 */
	scale(&p->F, dt / (2 * p->m), &tmp);
	add(&p->vel, &tmp, &p->vel);

	assert(!isnan(p->vel.x) && !isnan(p->vel.y) && !isnan(p->vel.z));

	/* pos(t + dt) = pos(t) + vel(t + dt/2)*dt */
	scale(&p->vel, dt, &tmp);
	add(&p->pos, &tmp, &p->pos);
}
static void verletHelper2(Particle *p)
{
	double dt = config.timeStep;
	Vec3 tmp;

	/* vel(t + dt) = vel(t + dt/2) + acc(t + dt)*dt/2 */
	scale(&p->F, dt / (2 * p->m), &tmp);
	add(&p->vel, &tmp, &p->vel);
}
static void verlet(void)
{
	// The compiler better inlines all of this. TODO if not: force it.
	forEveryParticle(&verletHelper1);
	calculateForces(); /* acc(t + dt) */
	forEveryParticle(&verletHelper2);
	reboxParticles(); //TODO only once every N iterations...
}

static void langevinBBKhelper1(Particle *p)
{
	double dt    = config.timeStep;
	double gamma = config.langevinGamma;

	Vec3 tmp1, tmp2;

	/* from v(t) to v(t + dt/2) */
	scale(&p->vel, 1 - gamma*dt/2, &tmp1);

	/* p->F is regular force + random collision force */
	scale(&p->F, dt / (2 * p->m), &tmp2);

	add(&tmp1, &tmp2, &p->vel);

	/* from r(t) to r(t + dt) */
	scale(&p->vel, dt, &tmp1);
	add(&p->pos, &tmp1, &p->pos);
}
static void langevinBBKhelper2(Particle *p)
{
	double dt    = config.timeStep;
	double gamma = config.langevinGamma;
	double T     = config.thermostatTemp;

	/* Regular forces have been calculated. Add the random force due 
	 * to collisions to the total force. The result is:
	 * p->F = F(t + dt) + R(t + dt) */
	/* TODO check that compiler inlines this and precalculates the 
	 * prefactor before p->m when looping over all particles. */
	double Rstddev = sqrt(2 * BOLTZMANN_CONSTANT * T * gamma * p->m / dt);
	Vec3 R = randNormVec();
	scale(&R, Rstddev, &R);
	add(&p->F, &R, &p->F);

	/* from v(t + dt/2) to v(t + dt) */
	Vec3 tmp;
	scale(&p->F, dt / (2 * p->m), &tmp);
	add(&p->vel, &tmp, &tmp);
	scale(&tmp, 1 / (1 + gamma*dt/2), &p->vel);
}

/* BBK integrator for Langevin dynamics. Uses the one based on 
 * velocity-verlet to include calculation of the velocities. 
 * See http://localscf.com/LangevinDynamics.aspx */
static void langevinBBK(void)
{
	forEveryParticle(&langevinBBKhelper1); /* updates positions */
	reboxParticles(); //TODO only once every N iterations(?)
	calculateForces();
	forEveryParticle(&langevinBBKhelper2);
}



static void stepWorld(Integrator integrator)
{
	switch(integrator) {
	case VERLET:
		verlet();
		assert(physicsCheck());
		thermostat();
		assert(physicsCheck());
		break;
	case LANGEVIN:
		langevinBBK();
		break;
	default:
		fprintf(stderr, "ERROR: Unknown integrator!\n");
		assert(false);
		break;
	}
}

typedef struct
{
	Integrator integrator;
} IntegratorState;
/* The integrator task is responsible for handeling the space partition 
 * grid */
static void *integratorTaskStart(void *initialData)
{
	IntegratorConf *ic = (IntegratorConf*) initialData;

	if(!allocGrid(ic->numBoxes, config.worldSize))
		return NULL;

	forEveryParticle(&addToGrid);

	IntegratorState *state = malloc(sizeof(*state));
	state->integrator = ic->integrator;
	free(initialData);
	return state;
}
static bool integratorTaskTick(void *state)
{
	IntegratorState *is = (IntegratorState*) state;
	stepWorld(is->integrator);
	return true;
}
static void integratorTaskStop(void *state)
{
	UNUSED(state);
	freeGrid();
}

Task makeIntegratorTask(IntegratorConf *conf)
{
	Task task;
	IntegratorConf *confCpy = malloc(sizeof(*confCpy));
	memcpy(confCpy, conf, sizeof(*confCpy));

	task.initialData = confCpy;
	task.start = &integratorTaskStart;
	task.tick  = &integratorTaskTick;
	task.stop  = &integratorTaskStop;
	return task;
}



/* ===== INFORMATION FUNCTIONS ===== */

void dumpStats()
{
	PotentialEnergies pe = calcPotentialEnergies();
	double K = kineticEnergy() * ENERGY_FACTOR;
	double T = temperature();
	double E = K + pe.bond + pe.angle + pe.dihedral + pe.stack + pe.basePair + pe.Coulomb;

	printf("E = %e, K = %e, Vb = %e, Va = %e, Vd = %e, Vs = %e, Vbp = %e, Vpp = %e, T = %f\n",
			E, K, pe.bond, pe.angle, pe.dihedral, pe.stack, pe.basePair, pe.Coulomb, T);
}

void dumpEnergies(FILE *stream)
{
#if 1
	assert(stream != NULL);
	PotentialEnergies pe = calcPotentialEnergies();
	double K = kineticEnergy() * ENERGY_FACTOR;
	double E = K + pe.bond + pe.angle + pe.dihedral + pe.stack + pe.basePair + pe.Coulomb;
	fprintf(stream, "%e %e %e %e %e %e %e %e %e \n",
			getTime(), E, K, pe.bond, pe.angle, pe.dihedral, pe.stack, pe.basePair, pe.Coulomb);
#else
	/* DEBUG equipartition theorem */
	dumpEquipartitionStats();
#endif
}


/* ===== MISC FUNCTIONS ===== */

Vec3 getCOM(Particle *ps, int num)
{
	Vec3 COM = {0, 0, 0};
	double M = 0; /* total mass */
	for (int i = 0; i < num; i++) {
		Vec3 tmp;
		scale(&ps[i].pos, ps[i].m, &tmp);
		add(&tmp, &COM, &COM);
		M += ps[i].m;
	}
	scale(&COM, 1/M, &COM);
	return COM;
}


