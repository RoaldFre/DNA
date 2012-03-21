#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "system.h"
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
#define BOND_K1      (EPSILON * FROM_ANGSTROM_SQUARED)
#define BOND_K2      (100 * EPSILON * FROM_ANGSTROM_SQUARED)
/* Bond bend */
#define BOND_Ktheta  (400 * EPSILON) /* per radian^2 */
/* Bond twist */
#define BOND_Kphi    (4 * EPSILON)
/* Bond stack */
#define BOND_STACK   EPSILON


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
#define COUPLING_BP_AT 	19.28e-21 /* 2.77 kcal/mol == 19.28e-21 J (per particle) */
#define COUPLING_BP_GC	28.96e-21 /* 4.16 kcal/mol == 28.96e-21 J (per particle) */
#define DISTANCE_r0_AT	2.9002e-10	/* Knotts et al 2007, table III, 2.9002 Angstrom */
#define DISTANCE_r0_GC	2.8694e-10	/* Knotts et al 2007, table III, 2.8694 Angstrom */

#define ENERGY_FACTOR	(1/1.602177e-19) /* Energy in electronvolt */
#define BOLTZMANN_CONSTANT    1.38065e-23
#define FROM_ANGSTROM_SQUARED 1e20 /* Bond constants are given for angstrom */
#define TO_RADIANS	      (M_PI / 180)

static void verlet(void);
static void calculateForces(void);

static double randNorm(void);

static double kineticEnergy(void);
static double Vbond(Particle *p1, Particle *p2, double d0);
static void   Fbond(Particle *p1, Particle *p2, double d0);
static double Vstack(Particle *p1, Particle *p2);
static void   Fstack(Particle *p1, Particle *p2);
static double Vangle(Particle *p1, Particle *p2, Particle *p3, double theta0);
static void   Fangle(Particle *p1, Particle *p2, Particle *p3, double theta0);
static double Vdihedral(Particle*, Particle*, Particle*, Particle*, double);
static void   Fdihedral(Particle*, Particle*, Particle*, Particle*, double);
static void   FdihedralParticle(Particle *target, Particle *p1, Particle *p2,
			Particle *p3, Particle *p4, double Vorig, double phi0);

static void basePairForce(Particle *p1, Particle *p2);

/* GLOBALS */

World world;
Config config;

double sim_time = 0;


/* Allocates the world to hold the given number of strands, and sets this 
 * value in world.
 * Precondition: allocWorld must not already have been called (unless 
 * followed by a freeWorld). */
bool allocWorld(int numStrands, int numBoxes, double worldSize)
{
	assert(world.strands == NULL);
	world.strands = calloc(numStrands, sizeof(*world.strands));
	if (world.strands == NULL)
		return false;
	world.numStrands = numStrands;


	if(!allocGrid(numBoxes, worldSize)) {
		freeWorld();
		return false;
	}

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
	s->Bs = &s->all[1 * config.numMonomers];
	s->Ps = &s->all[2 * config.numMonomers];

	s->numMonomers = numMonomers;
	return true;
}


/* Place monomers in a vertical column (in the x-y plane) in the center of 
 * the world. Distances between sugar, base and phospate are the 
 * equilibrium lenghts with some small gaussian jitter added.
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
void fillWorld()
{
	int n = config.numMonomers;

	double spacing = D_S5P + D_S3P; /* vertical spacing between monomers */
	double yoffset = -n * spacing / 2 + config.worldSize/2;
	double xoffset = -D_SA / 2 + config.worldSize/2;
	double posStdev = spacing / 100;

	for (int s = 0; s < world.numStrands; s++) {
		Strand strand = world.strands[s];
		for (int i = 0; i < n; i++) {
			/* Positions */
			strand.Ss[i].pos.z = strand.Ps[i].pos.z = strand.Bs[i].pos.z 
					= config.worldSize / 2;

			strand.Ss[i].pos.x = xoffset;
			strand.Bs[i].pos.x = xoffset + D_SA;
			strand.Ps[i].pos.x = xoffset;

			strand.Ss[i].pos.y = yoffset + i*spacing;
			strand.Bs[i].pos.y = yoffset + i*spacing;
			strand.Ps[i].pos.y = yoffset + i*spacing + D_S5P;

			strand.Ss[i].pos.x += posStdev * randNorm();
			strand.Ss[i].pos.y += posStdev * randNorm();
			strand.Ss[i].pos.z += posStdev * randNorm();
			strand.Bs[i].pos.x += posStdev * randNorm();
			strand.Bs[i].pos.y += posStdev * randNorm();
			strand.Bs[i].pos.z += posStdev * randNorm();
			strand.Ps[i].pos.x += posStdev * randNorm();
			strand.Ps[i].pos.y += posStdev * randNorm();
			strand.Ps[i].pos.z += posStdev * randNorm();

			/* Velocity */
			strand.Ss[i].vel.x = strand.Ss[i].vel.y = strand.Ss[i].vel.z = 0;
			strand.Bs[i].vel.x = strand.Bs[i].vel.y = strand.Bs[i].vel.z = 0;
			strand.Ps[i].vel.x = strand.Ps[i].vel.y = strand.Ps[i].vel.z = 0;

			/* Mass */
			strand.Ss[i].m = MASS_S;
			strand.Bs[i].m = MASS_A;
			strand.Ps[i].m = MASS_P;

			/* Type */
			strand.Ss[i].type = SUGAR;
			strand.Ps[i].type = PHOSPHATE;
			
			int basetype_number = rand() % 4;
			
			switch (basetype_number) {
				case 0: strand.Bs[i].type = BASE_A; break;
				case 1: strand.Bs[i].type = BASE_T; break;
				case 2: strand.Bs[i].type = BASE_C; break;
				case 3: strand.Bs[i].type = BASE_G; break;
			}
		}
	}

	/* Add everything to the space partition grid */
	forEveryParticle(&addToGrid);
}

/* Returns a number sampled from a standard normal distribution. */
static double randNorm()
{
	/* Box-Muller transform */
	double u1 = ((double) rand()) / RAND_MAX;
	double u2 = ((double) rand()) / RAND_MAX;

	return sqrt(-2*log(u1)) * cos(2*M_PI*u2);
}

void freeWorld()
{
	assert(world.strands != NULL);
	for (int s = 0; s < world.numStrands; s++)
		freeStrand(&world.strands[s]);
	free(world.strands);
	freeGrid();
	return;
}
void freeStrand(Strand *strand)
{
	free(strand->all);
	strand->numMonomers = 0;
}

/* loop over every particle in the world */
void forEveryParticle(void (*f)(Particle *p))
{
	for (int s = 0; s < world.numStrands; s++) {
		Strand *strand = &world.strands[s];
		for (int i = 0; i < 3*strand->numMonomers; i++)
			f(&strand->all[i]);
	}
}
/* loop over every particle in the world, pass [D]ata to the function */
void forEveryParticleD(void (*f)(Particle *p, void *data), void *data)
{
	for (int s = 0; s < world.numStrands; s++) {
		Strand *strand = &world.strands[s];
		for (int i = 0; i < 3*strand->numMonomers; i++)
			f(&strand->all[i], data);
	}
}




/* PHYSICS */

//TODO Verlet implementation is wrong, iirc; but needs to be changed by Langevin anyway.
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
static void verletHelper2(Particle *p) {
	double dt = config.timeStep;
	Vec3 tmp;

	/* vel(t + dt) = vel(t + dt/2) + acc(t + dt)*dt/2 */
	scale(&p->F, dt / (2 * p->m), &tmp);
	add(&p->vel, &tmp, &p->vel);
}
static void verlet()
{
	// The compiler better inlines all of this. TODO if not: force it.
	forEveryParticle(&verletHelper1);
	calculateForces(); /* acc(t + dt) */
	forEveryParticle(&verletHelper2);
}



static double temperature(void)
{
	return 2.0 / (3.0 * BOLTZMANN_CONSTANT)
			* kineticEnergy() / (config.numMonomers * 3.0);
}

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

static void basePairForce(Particle *p1, Particle *p2)
{
	double rij = distance(&p1->pos, &p2->pos);
	double bp_coupling;
	double bp_force_distance;
	double force;
	Vec3 force_vec;
	
	/* Apply right force constant for AT-bonding and GC-bonding, 
	 * if not AT or GC then zero */
	
	/* Leonard-Jones potential: force = bp_coupling 
	 *  * 5*(r^0 / r)^12 - 6*(r^0/r)^10 + 1 ].
	 * 
	 * For the force we differentiate with respect to r,so we get
	 * bp_coupling * 60*[ r^0^10 / r^11 - r^0^12/r^13 ] */	
	
	if ((p1->type==BASE_A && p2->type==BASE_T) 
			|| (p1->type==BASE_T && p2->type==BASE_A)) {
		bp_coupling = COUPLING_BP_AT;
		bp_force_distance = DISTANCE_r0_AT;
		
	} else if ((p1->type==BASE_G && p2->type==BASE_C) 
			|| (p1->type==BASE_C && p2->type==BASE_G)) {
		bp_coupling = COUPLING_BP_GC;
		bp_force_distance = DISTANCE_r0_GC;
	} else {
		return; /* no force */
	}

	double rfrac = bp_force_distance / rij;
	double rfrac2 = rfrac * rfrac;
	double rfrac4 = rfrac2 * rfrac2;
	double rfrac8 = rfrac4 * rfrac4;
	double rfrac10 = rfrac8 * rfrac2;
	double rfrac12 = rfrac10 * rfrac2;
				
	force = bp_coupling*60*( rfrac12 / rij - rfrac10 / rij );
	
	/* calculate the direction of the force between the basepairs */
	Vec3 direction;
	
	/* subtract and normalize */
	sub(&p1->pos, &p2->pos, &direction);
	normalize(&direction, &direction);
	
	/* scale the direction with the calculated force */
	scale(&direction, force, &force_vec);

	/* add force to particle objects */
	add(&p1->F, &force_vec, &p1->F);
	sub(&p2->F, &force_vec, &p2->F);
}


static void calculateForces()
{
	/* Reset forces */
	forEveryParticle(&resetForce);

	/* Strand-based forces */
	for (int s = 0; s < world.numStrands; s++)
		strandForces(&world.strands[s]);

	/* Particle-based forces */
	//forEveryPair(&basePairForce);
}





/* V = k1 * (dr - d0)^2  +  k2 * (d - d0)^4
 * where dr is the distance between the particles */
static double Vbond(Particle *p1, Particle *p2, double d0)
{
	double k1 = BOND_K1;
	double k2 = BOND_K2;
	double d = distance(&p1->pos, &p2->pos) - d0;
	double d2 = d * d;
	double d4 = d2 * d2;
	return k1 * d2  +  k2 * d4;
}
static void Fbond(Particle *p1, Particle *p2, double d0)
{
	double k1 = BOND_K1;
	double k2 = BOND_K2;
	Vec3 drVec, drVecNormalized, F;
	sub(&p2->pos, &p1->pos, &drVec);
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
	sub(&p1->pos, &p2->pos, &a);
	sub(&p3->pos, &p2->pos, &b);
	double dtheta = angle(&a, &b) - theta0;
	return ktheta/2 * dtheta*dtheta;
}
static void Fangle(Particle *p1, Particle *p2, Particle *p3, double theta0)
{
	Vec3 a, b;
	double ktheta = BOND_Ktheta;
	sub(&p1->pos, &p2->pos, &a);
	sub(&p3->pos, &p2->pos, &b);
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


void stepWorld(void)
{
	verlet();
	assert(physicsCheck());
	thermostat();
	assert(physicsCheck());
	sim_time += config.timeStep;
}

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
	double PPM = length(&P) / config.numMonomers;
	if (PPM > 1e-20) {
		fprintf(stderr, "\nMOMENTUM CONSERVATION VIOLATED! "
				"Momentum per monomer: |P| = %e\n", PPM);
		return false;
	}
	return true;
}

typedef struct potentialEnergies {
	double bond, angle, dihedral, stack;
} PotentialEnergies;

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
	for (int i = 1; i < config.numMonomers; i++) {
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

	pe->bond     = Vb * ENERGY_FACTOR;
	pe->angle    = Va * ENERGY_FACTOR;
	pe->dihedral = Vd * ENERGY_FACTOR;
	pe->stack    = Vs * ENERGY_FACTOR;
}

/* Return energy stats of world, in electronvolts. */
static PotentialEnergies calcPotentialEnergies(void) {
	PotentialEnergies pe = {0, 0, 0, 0};
	for (int s = 0; s < world.numStrands; s++)
		addPotentialEnergies(&world.strands[s], &pe);
	return pe;
}

void dumpStats()
{
	PotentialEnergies pe = calcPotentialEnergies();
	double K = kineticEnergy() * ENERGY_FACTOR;
	double T = temperature();
	double E = K + pe.bond + pe.angle + pe.dihedral + pe.stack;

	printf("E = %e, K = %e, Vb = %e, Va = %e, Vd = %e, Vs = %e, T = %f\n",
			E, K, pe.bond, pe.angle, pe.dihedral, pe.stack, T);
}

void dumpEnergies(FILE *stream)
{
#if 1
	assert(stream != NULL);
	PotentialEnergies pe = calcPotentialEnergies();
	double K = kineticEnergy() * ENERGY_FACTOR;
	double E = K + pe.bond + pe.angle + pe.dihedral + pe.stack;
	fprintf(stream, "%e %e %e %e %e %e %e\n",
			sim_time, E, K, pe.bond, pe.angle, pe.dihedral, pe.stack);
#else
	/* DEBUG equipartition theorem */
	dumpEquipartitionStats();
#endif
}

