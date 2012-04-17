#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "physics.h"
#include "world.h"
#include "spgrid.h"



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


static double calcVBasePair(double coupling, double rij0, double rijVar2)
{
	double rfrac2 = rij0 * rij0 / rijVar2;
	double rfrac4 = rfrac2 * rfrac2;
	double rfrac8 = rfrac4 * rfrac4;
	double rfrac10 = rfrac8 * rfrac2;
	double rfrac12 = rfrac8 * rfrac4;
	
	return coupling*(5*rfrac12 - 6*rfrac10 + 1);
}
double VbasePair(Particle *p1, Particle *p2)
{
	double rij = nearestImageDistance(&p1->pos, &p2->pos);
	if (rij > config.truncationLen)
		return 0; /* Too far away */

	double rij2 = rij*rij;
	
	double bpCoupling;
	double bpForceDist;
	
	/* Apply right potential constant for AT-bonding and GC-bonding, 
	 * if not AT or GC then zero */
	
	/* Lennard-Jones potential:
	 * potential = bpCoupling * 5*(r0 / r)^12 - 6*(r0/r)^10 + 1. */
	if ((p1->type==BASE_A && p2->type==BASE_T) 
			|| (p1->type==BASE_T && p2->type==BASE_A)) {
		bpCoupling = COUPLING_BP_AT;
		bpForceDist = DISTANCE_r0_AT;
		
	} else if ((p1->type==BASE_G && p2->type==BASE_C) 
			|| (p1->type==BASE_C && p2->type==BASE_G)) {
		bpCoupling = COUPLING_BP_GC;
		bpForceDist = DISTANCE_r0_GC;
	} else {
		return 0; /* no L-J potential */
	}
	double truncLen2 = config.truncationLen * config.truncationLen;
	double truncCorrection = calcVBasePair(bpCoupling, bpForceDist, truncLen2);
	double bpPotential = calcVBasePair(bpCoupling, bpForceDist, rij2) - truncCorrection;

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



static void stepPhysics(Integrator integrator)
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

	if(!allocGrid(ic->numBoxes, world.worldSize))
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
	stepPhysics(is->integrator);
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


