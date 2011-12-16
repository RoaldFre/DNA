#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "vmath.h"
#include "system.h"

/* Masses (in kg) */
#define AU      1.660539e-27
#define MASS_A  (134.1 * AU)
#define MASS_P  (94.97 * AU)
#define MASS_S  (83.11 * AU)

/* Equilibrium distance of bonds (in m) */
#define D_S5P   3.899e-10
#define D_S3P   3.559e-10
#define D_SA    6.430e-10

/* Bond stretch */
#define TO_ANGSTROM_SQUARED 1e-20 /* Bond constants are given for angstrom */
#define BOND_K1      (EPSILON / TO_ANGSTROM_SQUARED)
#define BOND_K2      (100 * EPSILON / TO_ANGSTROM_SQUARED)
/* Bond bend */
#define BOND_Ktheta  (400 * EPSILON) /* per radian^2 */
/* Bond twist */
#define BOND_Kphi    (4 * EPSILON)

#define EPSILON 1.81e-21 /* 0.26kcal/mol == 1.81 * 10^-21 J (per particle) */

#define TO_RADIANS	(M_PI / 180)
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

#define ENERGY_FACTOR	(1/1.602177e-19) /* Energy in electronvolt */
#define MOMENTUM_FACTOR	(1/1.602177e-19) /* Momentum factor */

#define BOLTZMANN_CONSTANT 1.38065e-23

static void verlet(void);
static void calculateForces(void);

static double randNorm(void);
static void   randNormVec(double stdev, Vec3 *vec);

static double kineticEnergy(void);
static double Vbond(Particle *p1, Particle *p2, double d0);
static void   Fbond(Particle *p1, Particle *p2, double d0);
static double Vangle(Particle *p1, Particle *p2, Particle *p3, double theta0);
static void   Fangle(Particle *p1, Particle *p2, Particle *p3, double theta0);
static double Vdihedral(Particle*, Particle*, Particle*, Particle*, double);
static void   Fdihedral(Particle*, Particle*, Particle*, Particle*, double);
static void FdihedralParticle(Particle *target, Particle *p1, Particle *p2,
			Particle *p3, Particle *p4, double Vorig, double phi0);


/* GLOBALS */

World world;
Config config;

double time = 0;



/* Allocates the world.
 * Precondition: config MUST be valid, and allocWorld must not already have 
 * been called (unless followed by a freeWorld)
 * Returns true on success, false on failure. In the case of failure, 
 * nothing will be allocated */
bool allocWorld()
{
	assert(world.all == NULL);

	/* Allocate one big continuous list */
	world.all = calloc(3 * config.numMonomers, sizeof(*world.Ps));
	if (world.all == NULL)
		return false;
	
	/* Split the list in three sublists */
	world.Ss = &world.all[0];
	world.As = &world.all[1 * config.numMonomers];
	world.Ps = &world.all[2 * config.numMonomers];
	return true;
}

/* Place monomers in a vertical column (in the x-y plane) in the center of 
 * the world. Distances between sugar, base and phospate are the 
 * equilibrium lenghts with some small gaussian jitter added.
 *
 * TODO: not yet, now just simple gaussian sampling.
 * Initial velocities are sampled from a normal distribution to get a 
 * temperature of config.temperature 'per monomer'.
 * Velocities are then shifted to set the total momentum to zero. */

/*
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
 *        |      Ss[1]------As[1]     <-- i=1
 *        |    3'|
 *        |      |  . . . . . . . . . . . . . . . . . . . . 
 *        |      |                                       /|\
 *        |      Ps[0]                                    |
 *        |      |                                        |   one
 *        |      |                                        |  monomer
 *        |    5'|                                        |  
 *        |      Ss[0]------As[0]     <-- i=0            \|/
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
	double ws = config.worldSize;
	Vec3 totP = {0, 0, 0}, avgP; /* momentum */
	Vec3 tmp;
	Vec3 corrVelS, corrVelA, corrVelP; /* correction on velocity to get P=0 */

	double spacing = D_S5P + D_S3P; /* vertical spacing between monomers */
	double yoffset = (ws - n * spacing) / 2;
	double xoffset = (ws - D_SA) / 2;
	double posStdev = spacing / 100;
	//XXX XXX
	//TODO: boltzmann cte etc etc
	double velStdev = 0 * sqrt(config.thermostatTemp); // TODO factor 1/3 (3D) somethere?

	for (int i = 0; i < n; i++) {
		/* Positions */
		world.Ss[i].pos.z = world.Ps[i].pos.z = world.As[i].pos.z = ws / 2;

		world.Ss[i].pos.x = xoffset;
		world.As[i].pos.x = xoffset + D_SA;
		world.Ps[i].pos.x = xoffset;

		world.Ss[i].pos.y = yoffset + i*spacing;
		world.As[i].pos.y = yoffset + i*spacing;
		world.Ps[i].pos.y = yoffset + i*spacing + D_S5P;

		world.Ss[i].pos.x += posStdev * randNorm();
		world.Ss[i].pos.y += posStdev * randNorm();
		world.Ss[i].pos.z += posStdev * randNorm();
		world.As[i].pos.x += posStdev * randNorm();
		world.As[i].pos.y += posStdev * randNorm();
		world.As[i].pos.z += posStdev * randNorm();
		world.Ps[i].pos.x += posStdev * randNorm();
		world.Ps[i].pos.y += posStdev * randNorm();
		world.Ps[i].pos.z += posStdev * randNorm();

		/* Velocities */
		randNormVec(velStdev, &world.Ss[i].vel);
		randNormVec(velStdev, &world.As[i].vel);
		randNormVec(velStdev, &world.Ps[i].vel);

		scale(&world.Ss[i].vel, MASS_S, &tmp);
		add(&totP, &tmp, &totP);
		scale(&world.As[i].vel, MASS_A, &tmp);
		add(&totP, &tmp, &totP);
		scale(&world.Ps[i].vel, MASS_P, &tmp);
		add(&totP, &tmp, &totP);

		/* Mass */
		world.Ss[i].m = MASS_S;
		world.As[i].m = MASS_A;
		world.Ps[i].m = MASS_P;
	}
	/* Correct for zero momentum */
	scale(&totP, 1.0 / (3 * n), &avgP);
	scale(&avgP, 1.0 / MASS_S, &corrVelS);
	scale(&avgP, 1.0 / MASS_A, &corrVelA);
	scale(&avgP, 1.0 / MASS_P, &corrVelP);
	for (int i = 0; i < config.numMonomers; i++) {
		sub(&world.Ss[i].vel, &corrVelS, &world.Ss[i].vel);
		sub(&world.As[i].vel, &corrVelA, &world.As[i].vel);
		sub(&world.Ps[i].vel, &corrVelP, &world.Ps[i].vel);
	}
}

/* Returns a number sampled from a standard normal distribution. */
static double randNorm()
{
	/* Box-Muller transform */
	double u1 = ((double) rand()) / RAND_MAX;
	double u2 = ((double) rand()) / RAND_MAX;

	return sqrt(-2*log(u1)) * cos(2*M_PI*u2);
}

static void randNormVec(double stdev, Vec3 *vec)
{
	vec->x = stdev * randNorm();
	vec->y = stdev * randNorm();
	vec->z = stdev * randNorm();
}

void freeWorld()
{
	free(world.Ss);
	free(world.As);
	free(world.Ps);

	return;
}

void dumpWorld()
{
	//TODO
	return;
}







/* PHYSICS */

static void verlet()
{
	double dt = config.timeStep;

	/* Velocity Verlet */
	for (int i = 0; i < 3 * config.numMonomers; i++) {
		Particle *p = &world.all[i];
		Vec3 tmp;

		/* vel(t + dt/2) = vel(t) + acc(t)*dt/2 */
		scale(&p->F, dt / (2 * p->m), &tmp);
		add(&p->vel, &tmp, &p->vel);

		assert(!isnan(p->vel.x) && !isnan(p->vel.y) && !isnan(p->vel.z));

		/* pos(t + dt) = pos(t) + vel(t + dt/2)*dt */
		scale(&p->vel, dt, &tmp);
		add(&p->pos, &tmp, &p->pos);
	}
	calculateForces(); /* acc(t + dt) */
	for (int i = 0; i < 3 * config.numMonomers; i++) {
		Particle *p = &world.all[i];
		Vec3 tmp;

		/* vel(t + dt) = vel(t + dt/2) + acc(t + dt)*dt/2 */
		scale(&p->F, dt / (2 * p->m), &tmp);
		add(&p->vel, &tmp, &p->vel);
	}
}

static double temperature()
{
	return 2.0 / (3.0 * BOLTZMANN_CONSTANT)
			* kineticEnergy() / (config.numMonomers * 3.0);
}

static void thermostat()
{
	if (config.thermostatTau <= 0)
		return;

	/* Mass and Boltzmann constant are 1 */ 
	double Tk  = temperature();
	double T0  = config.thermostatTemp;
	double dt  = config.timeStep;
	double tau = config.thermostatTau;
	double lambda = sqrt(1 + dt/tau * (T0/Tk - 1));

	for (int i = 0; i < config.numMonomers * 3; i++) {
		Particle *p = &world.all[i];
		scale(&p->vel, lambda, &p->vel);
	}
}

static void calculateForces()
{
	World *w = &world;

	/* Reset forces */
	for (int i = 0; i < 3 * config.numMonomers; i++) {
		w->all[i].F.x = 0;
		w->all[i].F.y = 0;
		w->all[i].F.z = 0;
	}

	/* Bottom monomer */
	Fbond(&w->Ss[0], &w->As[0], D_SA);
	Fbond(&w->Ss[0], &w->Ps[0], D_S5P);
	Fangle(&w->Ps[0], &w->Ss[0], &w->As[0], ANGLE_P_5S_A);
	/* Rest of the monomers */
	for (int i = 1; i < config.numMonomers; i++) {
		Fbond(&w->Ss[i], &w->As[i],   D_SA);
		Fbond(&w->Ss[i], &w->Ps[i],   D_S5P);
		Fbond(&w->Ss[i], &w->Ps[i-1], D_S3P);

		Fangle(&w->Ps[ i ], &w->Ss[ i ], &w->As[ i ], ANGLE_P_5S_A);
		Fangle(&w->Ps[ i ], &w->Ss[ i ], &w->Ps[i-1], ANGLE_P_5S3_P);
		Fangle(&w->Ps[i-1], &w->Ss[ i ], &w->As[ i ], ANGLE_P_3S_A);
		Fangle(&w->Ss[i-1], &w->Ps[i-1], &w->Ss[ i ], ANGLE_S5_P_3S);

		Fdihedral(&w->Ps[ i ], &w->Ss[ i ], &w->Ps[i-1], &w->Ss[i-1],
							DIHEDRAL_P_5S3_P_5S);
		Fdihedral(&w->As[ i ], &w->Ss[ i ], &w->Ps[i-1], &w->Ss[i-1],
							DIHEDRAL_A_S3_P_5S);
		Fdihedral(&w->Ss[ i ], &w->Ps[i-1], &w->Ss[i-1], &w->As[i-1],
							DIHEDRAL_S3_P_5S_A);
		if (i >= 2)
		Fdihedral(&w->Ss[ i ], &w->Ps[i-1], &w->Ss[i-1], &w->Ps[i-2],
							DIHEDRAL_S3_P_5S3_P);
	}
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

	Vec3 tmp1, tmp2, F1, F2, F3, F4;

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

	assert(fabs(dot(&a, &F1) / length(&a) / length(&F1)) < 1e-5);
	assert(fabs(dot(&b, &F3) / length(&b) / length(&F3)) < 1e-5);
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
	F.x = (Vdihedral(p1, p2, p3, p4, phi0) - Vorig) / h;
	target->pos.x -= h;

	h = target->pos.y * hfactor;
	target->pos.y += h;
	F.y = (Vdihedral(p1, p2, p3, p4, phi0) - Vorig) / h;
	target->pos.y -= h;

	h = target->pos.z * hfactor;
	target->pos.z += h;
	F.z = (Vdihedral(p1, p2, p3, p4, phi0) - Vorig) / h;
	target->pos.z -= h;

	add(&target->F, &F, &target->F);
}

void stepWorld(void)
{
	verlet();
	assert(physicsCheck());
	thermostat();
	assert(physicsCheck());
	time += config.timeStep;
}

static double kineticEnergy(void)
{
	double twiceK = 0;
	for (int i = 0; i < 3 * config.numMonomers; i++)
		twiceK += world.all[i].m * length2(&world.all[i].vel);
	return twiceK/2;
}
		
static Vec3 momentum(void)
{
	Vec3 Ptot = {0, 0, 0};
	for (int i = 0; i < 3 * config.numMonomers; i++) {
		Vec3 P;
		scale(&world.all[i].vel, world.all[i].m, &P);
		add(&P, &Ptot, &Ptot);
	}
	return Ptot;
}

bool physicsCheck(void)
{
	Vec3 P = momentum();
	return length(&P) < 1e-10;
}

void dumpStats()
{
	World *w = &world;

	double K = kineticEnergy();

	double Vb = 0;
	double Va = 0;
	double Vd = 0;
	Vb += Vbond(&w->Ss[0], &w->As[0],   D_SA);
	Vb += Vbond(&w->Ss[0], &w->Ps[0],   D_S5P);

	Va += Vangle(&w->As[0], &w->Ss[0], &w->Ps[0], ANGLE_P_5S_A);
	for (int i = 1; i < config.numMonomers; i++) {
		Vb += Vbond(&w->Ss[i], &w->As[i],   D_SA);
		Vb += Vbond(&w->Ss[i], &w->Ps[i],   D_S5P);
		Vb += Vbond(&w->Ss[i], &w->Ps[i-1], D_S3P);

		Va += Vangle(&w->Ps[ i ], &w->Ss[ i ], &w->As[ i ], ANGLE_P_5S_A);
		Va += Vangle(&w->Ps[ i ], &w->Ss[ i ], &w->Ps[i-1], ANGLE_P_5S3_P);
		Va += Vangle(&w->Ps[i-1], &w->Ss[ i ], &w->As[ i ], ANGLE_P_3S_A);
		Va += Vangle(&w->Ss[i-1], &w->Ps[i-1], &w->Ss[ i ], ANGLE_S5_P_3S);

		Vd += Vdihedral(&w->Ps[ i ], &w->Ss[ i ], &w->Ps[i-1], &w->Ss[i-1],
							DIHEDRAL_P_5S3_P_5S);
		Vd += Vdihedral(&w->As[ i ], &w->Ss[ i ], &w->Ps[i-1], &w->Ss[i-1],
							DIHEDRAL_A_S3_P_5S);
		Vd += Vdihedral(&w->Ss[ i ], &w->Ps[i-1], &w->Ss[i-1], &w->As[i-1],
							DIHEDRAL_S3_P_5S_A);
		if (i >= 2)
		Vd += Vdihedral(&w->Ss[ i ], &w->Ps[i-1], &w->Ss[i-1], &w->Ps[i-2],
							DIHEDRAL_S3_P_5S3_P);
	}

	K  *= ENERGY_FACTOR;
	Vb *= ENERGY_FACTOR;
	Va *= ENERGY_FACTOR;
	Vd *= ENERGY_FACTOR;

	double T = temperature();

	Vec3 P = momentum();
	double normP = length(&P) * MOMENTUM_FACTOR;
	printf("E = %e, K = %e, Vb = %e, Va = %e, Vd = %e, |P| = %e, T = %f\n",
			K+Vb+Va+Vd, K, Vb, Va, Vd, normP, T);
}
