#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "physics.h"
#include "world.h"
#include "spgrid.h"

/* Disable interactions by commenting these defines */
#define ENABLE_BOND		true
#define ENABLE_ANGLE		true
#define ENABLE_DIHEDRAL		true
#define ENABLE_STACK		true //TODO check whether distances are correct
#define ENABLE_EXCLUSION	false //TODO SERIOUSLY FUCKED UP
#define ENABLE_BASE_PAIR	true
#define ENABLE_COULOMB		true



/* ===== FORCES AND POTENTIALS ===== */

/* BOND */
/* V = k1 * (dr - d0)^2  +  k2 * (d - d0)^4
 * where dr is the distance between the particles */
static double Vbond(Particle *p1, Particle *p2, double d0)
{
	if (!ENABLE_BOND)
		return 0;
	double k1 = BOND_K1;
	double k2 = BOND_K2;
	double d = nearestImageDistance(p1->pos, p2->pos) - d0;
	double d2 = d * d;
	double d4 = d2 * d2;
	return k1 * d2  +  k2 * d4;
}
static void Fbond(Particle *p1, Particle *p2, double d0)
{
	if (!ENABLE_BOND)
		return;
	double k1 = BOND_K1;
	double k2 = BOND_K2;
	Vec3 drVec = nearestImageVector(p1->pos, p2->pos);
	double dr = length(drVec);
	double d  = dr - d0;
	double d3 = d * d * d;

	Vec3 drVecNormalized = scale(drVec, 1/dr);
	Vec3 F = scale(drVecNormalized, 2*k1*d + 4*k2*d3);

	p1->F = add(p1->F, F);
	p2->F = sub(p2->F, F);
}

/* ANGLE */
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
	if (!ENABLE_ANGLE)
		return 0;
	Vec3 a, b;
	double ktheta = BOND_Ktheta;
	
	a = nearestImageVector(p2->pos, p1->pos);
	b = nearestImageVector(p2->pos, p3->pos);

	double dtheta = angle(a, b) - theta0;
	return ktheta/2 * dtheta*dtheta;
}
static void Fangle(Particle *p1, Particle *p2, Particle *p3, double theta0)
{
	if (!ENABLE_ANGLE)
		return;
	Vec3 a, b;
	double ktheta = BOND_Ktheta;
	
	a = nearestImageVector(p2->pos, p1->pos);
	b = nearestImageVector(p2->pos, p3->pos);
	
	double lal = length(a);
	double lbl = length(b);
	double adotb = dot(a, b);
	double costheta = adotb / (lal * lbl);
	double theta = acos(costheta);
	double sintheta = sqrt(1 - costheta*costheta);

	// TODO correct cut off?
	if (fabs(sintheta) < 1e-30)
		/* "No" force (unstable equilibrium), numerical instability 
		 * otherwise */
		return;

	Vec3 tmp1, tmp2, F1, F2, F3;

	tmp1 = scale(b, 1/(lal * lbl));
	tmp2 = scale(a, adotb / (lal*lal*lal * lbl));
	F1 = scale(sub(tmp1, tmp2), ktheta * (theta - theta0) / sintheta);
	p1->F = add(p1->F, F1);	

	tmp1 = scale(a, 1/(lal * lbl));
	tmp2 = scale(b, adotb / (lbl*lbl*lbl * lal));
	F3 = scale(sub(tmp1, tmp2), ktheta * (theta - theta0) / sintheta);
	p3->F = add(p3->F, F3);	

	F2 = add(F1, F3); /* actually -F2 */
	p2->F = sub(p2->F, F2);	

	assert(ktheta == 0 
		|| fabs(dot(a, F1) / length(a) / length(F1)) < 1e-5);
	assert(ktheta == 0
		|| fabs(dot(b, F3) / length(b) / length(F3)) < 1e-5);
}

/* DIHEDRAL */
static double Vdihedral(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
								double phi0)
{
	if (!ENABLE_DIHEDRAL)
		return 0;

	Vec3 r1 = nearestImageVector(p1->pos, p2->pos);
	Vec3 r2 = nearestImageVector(p2->pos, p3->pos);
	Vec3 r3 = nearestImageVector(p3->pos, p4->pos);
	
	double phi = dihedral(r1, r2, r3);
	return BOND_Kphi * (1 - cos(phi - phi0));
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

	target->F = add(target->F, F);
}
static void Fdihedral(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
								double phi0)
{
	if (!ENABLE_DIHEDRAL)
		return;

	/* This is a *mess* to do analytically, so we do a numerical 
	 * differentiation instead. */
	double Vorig = Vdihedral(p1, p2, p3, p4, phi0);
	FdihedralParticle(p1, p1, p2, p3, p4, Vorig, phi0);
	FdihedralParticle(p2, p1, p2, p3, p4, Vorig, phi0);
	FdihedralParticle(p3, p1, p2, p3, p4, Vorig, phi0);
	FdihedralParticle(p4, p1, p2, p3, p4, Vorig, phi0);
}

/* STACKING */
#define INV_CUBE_ROOT_OF_TWO 0.793700525984099737375852819636 /* 2^(-1/3) */
typedef struct {
	double r, phi, z;
} HelixInfo;
static HelixInfo getHelixInfo(ParticleType t)
{
	HelixInfo info = (HelixInfo) {0, 0, 0};
	switch (t) {
	case PHOSPHATE:	info.r = P_R;  info.phi = P_PHI;  info.z = P_Z;  break;
	case SUGAR:	info.r = S_R;  info.phi = S_PHI;  info.z = S_Z;  break;
	case BASE_A:	info.r = A_R;  info.phi = A_PHI;  info.z = A_Z;  break;
	case BASE_T:	info.r = T_R;  info.phi = T_PHI;  info.z = T_Z;  break;
	case BASE_C:	info.r = C_R;  info.phi = C_PHI;  info.z = C_Z;  break;
	case BASE_G:	info.r = G_R;  info.phi = G_PHI;  info.z = G_Z;  break;
	default: 	assert(false);
	}
	return info;
}
static double helixDistance2(HelixInfo bot, HelixInfo top,
		double delta_z, double delta_phi)
{
	double dz = top.z - bot.z + delta_z;
	double dphi = top.phi - bot.phi + delta_phi;
	Vec3 low  = (Vec3) {bot.r, 0, 0};
	Vec3 high = (Vec3) {top.r * cos(dphi), top.r * sin(dphi), dz};
	return distance2(low, high);
}

/* Distance between next neighbouring particles in the beta helix, ie 
 * particle i and i+2.
 * Particle 1 is the 'bottom' particle in the helix, particle 2 is the 
 * 'top' particle of the helix.
 * monomerDistance:
 *   1 for immediate neighbours i and i+1
 *   2 for   next    neighbours i and i+2
 *   n for           neighbours i and i+n */
static double neighbourStackDistance2(ParticleType t1, ParticleType t2, int monomerDistance)
{
	HelixInfo hi1 = getHelixInfo(t1);
	HelixInfo hi2 = getHelixInfo(t2);
	return helixDistance2(hi1, hi2,
			monomerDistance * HELIX_DELTA_Z,
			monomerDistance * HELIX_DELTA_PHI);
}
/* monomerDistance:
 *   1 for immediate neighbours (i and i+1)
 *   2 for   next    neighbours (i and i+2) */
static double Vstack(Particle *p1, Particle *p2, int monomerDistance)
{
	assert(isBase(p1->type) && isBase(p2->type));

	if (!ENABLE_STACK)
		return 0;

	double kStack = BOND_STACK;
	/* sigma = 2^(-1/6) * r_min  =>  sigma^2 = 2^(-1/3) * r_min^2 */
	double sigma2 = INV_CUBE_ROOT_OF_TWO * neighbourStackDistance2(
			p1->type, p2->type, monomerDistance);
	double sigma6 = sigma2 * sigma2 * sigma2;
	double sigma12 = sigma6 * sigma6;
	double r2 = distance2(p1->pos, p2->pos);
	double r6 = r2 * r2 * r2;
	double r12 = r6 * r6;

	return kStack * (sigma12/r12 - 2*sigma6/r6 + 1);
	//return kStack * (sigma12/r12 - 2*sigma6/r6);
}
static void Fstack(Particle *p1, Particle *p2, int monomerDistance)
{
	assert(isBase(p1->type) && isBase(p2->type));

	if (!ENABLE_STACK)
		return;

	double kStack = BOND_STACK;
	double sigma2 = INV_CUBE_ROOT_OF_TWO * neighbourStackDistance2(
			p1->type, p2->type, monomerDistance);
	double sigma6 = sigma2 * sigma2 * sigma2;
	double sigma12 = sigma6 * sigma6;

	Vec3 drVec = nearestImageVector(p1->pos, p2->pos);
	double dr2 = length2(drVec);

	assert(dr2 != 0);
	double dr4 = dr2 * dr2;
	double dr8 = dr4 * dr4;
	double dr14 = dr8 * dr4 * dr2;

	Vec3 F = scale(drVec, -12 * kStack * (sigma12/dr14 - sigma6/dr8));
	p1->F = add(p1->F, F);
	p2->F = sub(p2->F, F);
}


/* BASE PAIRING */
typedef struct {
	double coupling; /* Coupling strength */
	double distance2; /* Coupling distance squared */
} BasePairInfo;
/* Returns the BasePairInfo for the given types. If the given types do not 
 * constitute a valid base pair, then coupling and distance are set to -1. */
static BasePairInfo getBasePairInfo(ParticleType t1, ParticleType t2)
{
	BasePairInfo bpi;
	if ((t1 == BASE_A  &&  t2 == BASE_T) 
			||  (t1 == BASE_T  &&  t2 == BASE_A)) {
		bpi.coupling = COUPLING_BP_AT;
		bpi.distance2 = SQUARE(DISTANCE_r0_AT);
	} else if ((t1 == BASE_C  &&  t2 == BASE_G)
			||  (t1 == BASE_G  &&  t2 == BASE_C)) {
		bpi.coupling = COUPLING_BP_GC;
		bpi.distance2 = SQUARE(DISTANCE_r0_GC);
	} else {
		bpi.coupling = -1;
		bpi.distance2 = -1;
	}
	return bpi;
}
static bool isBondedBasePair(ParticleType t1, ParticleType t2)
{
	BasePairInfo bpi = getBasePairInfo(t1, t2);
	return bpi.coupling > 0;
}
static double calcVbasePair(BasePairInfo bpi, double rsquared)
{
	double rfrac2 = bpi.distance2 / rsquared;
	double rfrac4 = rfrac2 * rfrac2;
	double rfrac8 = rfrac4 * rfrac4;
	double rfrac10 = rfrac8 * rfrac2;
	double rfrac12 = rfrac8 * rfrac4;
	
	return bpi.coupling * (5*rfrac12 - 6*rfrac10 + 1);
}
double VbasePair(Particle *p1, Particle *p2)
{
	if (!ENABLE_BASE_PAIR)
		return 0;

	BasePairInfo bpi = getBasePairInfo(p1->type, p2->type);
	if (bpi.coupling < 0)
		return 0; /* Wrong pair */
	
	double rsq = nearestImageDistance2(p1->pos, p2->pos);
	double truncSq = SQUARE(config.truncationLen);
	if (rsq > truncSq)
		return 0; /* Too far away */

	return calcVbasePair(bpi, rsq) - calcVbasePair(bpi, truncSq);
}

static double calcFbasePair(BasePairInfo bpi, double r)
{
	double rfrac2 = bpi.distance2 / (r * r);
	double rfrac4 = rfrac2 * rfrac2;
	double rfrac8 = rfrac4 * rfrac4;
	double rfrac10 = rfrac8 * rfrac2;
	double rfrac12 = rfrac10 * rfrac2;
				
	return bpi.coupling * 60 * (rfrac12 - rfrac10) / r;
}
static void FbasePair(Particle *p1, Particle *p2)
{
	if (!ENABLE_BASE_PAIR)
		return;

	BasePairInfo bpi = getBasePairInfo(p1->type, p2->type);
	if (bpi.coupling < 0)
		return; /* Wrong pair */
	
	Vec3 rVec = nearestImageVector(p1->pos, p2->pos);
	double r = length(rVec);
	if (r > config.truncationLen)
		return; /* Too far away */

	Vec3 direction = scale(rVec, 1/r);

	double force = calcFbasePair(bpi, r);
	Vec3 forceVec = scale(direction, force);

	/* Add force to particle objects */
	p1->F = sub(p1->F, forceVec);
	p2->F = add(p2->F, forceVec);
}

static void Fexclusion(Particle *p1, Particle *p2)
{
	if (!ENABLE_EXCLUSION)
		return;
	double sig, rInv;
	Vec3 forceVec;
	double force;
	double cutOff;
	
	double rij = nearestImageDistance(p1->pos, p2->pos);
	if (rij > D_CUT)
		return; /* Too far away */

	Vec3 direction = nearestImageUnitVector(p1->pos, p2->pos);
	
	/* Apply right parameters for types of molecules, 
	 * if bases (and mismatched): SIGMA_0_CST*1.0
	 * otherwise: SIGMA_0_CST*D_CUT */
	
	/* Lennard-Jones potential:
	 * potential = exCoupling ((sigma0 / r)^12 - (sigma0/r)^6) + epsilon.
	 * 
	 * For the force we differentiate with respect to r, so we get
	 * 4*EPSILON*[ - 12* sigma0^12 / r^13 + 6* sigma0^6/r^7 ] */

	cutOff = getExclusionCutOff(p1->type, p2->type);
	
	
	sig = SIGMA_0_CST*cutOff;
	rInv = 1.0/rij;
	
	double sig12 = powerCalc(sig, 12);
	double sig6 = powerCalc(sig, 6);
	double rInv7 = powerCalc(rInv, 7);
	double rInv13 = powerCalc(rInv, 13);
	

	force = 4*EPSILON*(12*sig12*rInv13 - 6*sig6*rInv7); 
	
#if 0
	if (force > 1e-8){
	
	printf("Powers: distance: %e\t fractions: %e\t%e\n", rij, sig12*rInv13, sig6*rInv7);
	printf("Cutoff: \t%e\n", cutOff);
	printf("Force exclusion: %e\n", force);	
	}
#endif
		
	/* Scale the direction with the calculated force */
	forceVec = scale(direction, force);

	/* Add force to particle objects */
	p1->F = add(p1->F, forceVec);
	p2->F = sub(p2->F, forceVec);
}

double getExclusionCutOff(ParticleType t1, ParticleType t2){
	
	if (isBase(t1) && isBase(t2)){
		/* mis-matched basepairs */	
		return D_CUT_BASE;			
	} else {
		/* otherwise */
		return D_CUT;
	}
}

static double Vexclusion(Particle *p1, Particle *p2)
{
	if (!ENABLE_EXCLUSION)
		return 0;
	double sig, rInv;
	double potential;
	double cutOff;
	
	double rij = nearestImageDistance(p1->pos, p2->pos);
	if (rij > D_CUT)
		return 0; /* Too far away */

	
	/* Apply right parameters for types of molecules, 
	 * if bases and mismatched: SIGMA_0_CST*1.0
	 * if bases and matched: 0
	 * otherwise: SIGMA_0_CST*D_CUT */
	
	/* Lennard-Jones potential:
	 * potential = exCoupling ((sigma0 / r)^12 - (sigma0/r)^6) + epsilon.*/
	
	cutOff = getExclusionCutOff(p1->type, p2->type);

	sig = SIGMA_0_CST*cutOff;
	rInv = 1.0/rij;
	
	double sig12 = powerCalc(sig, 12);
	double sig6 = powerCalc(sig, 6);
	double rInv6 = powerCalc(rInv, 6);
	double rInv12 = powerCalc(rInv, 12);

	potential = 4*EPSILON*(sig12*rInv12 - sig6*rInv6) + EPSILON;
	
	/* correct so that at D_CUT, V zero */
	double dInv = 1.0/D_CUT;
	double dInv6 = powerCalc(dInv, 6);
	double dInv12 = powerCalc(dInv, 12);
	double correction = 4*EPSILON*(sig12*dInv12 - sig6*dInv6) + EPSILON;
	
	potential -= correction;
	
	return potential;
}


/* COULOMB */
static double calcInvDebyeLength(void)
{
	double T = config.thermostatTemp;
	double saltCon = config.saltConcentration;
	double lambdaBDenom, lambdaB;
	
	lambdaBDenom = 4 * M_PI * H2O_PERMETTIVITY
			* BOLTZMANN_CONSTANT * T;
	lambdaB = CHARGE_ELECTRON * CHARGE_ELECTRON / lambdaBDenom;
	return sqrt(8 * M_PI * lambdaB * AVOGADRO * saltCon);
}
static double calcVCoulomb(double r)
{
	double couplingConstant = CHARGE_ELECTRON * CHARGE_ELECTRON
				/ (4 * M_PI * H2O_PERMETTIVITY);
	double exponentialPart = exp(-r * calcInvDebyeLength());
	
	return couplingConstant * exponentialPart / r;
}
static bool isChargedPair(ParticleType t1, ParticleType t2)
{
	return t1 == PHOSPHATE  &&  t2 == PHOSPHATE;
}
static double VCoulomb(Particle *p1, Particle *p2)
{
	if (!ENABLE_COULOMB)
		return 0;

	if (!isChargedPair(p1->type, p2->type))
		return 0;

	double truncLength = config.truncationLen;
	double rij = nearestImageDistance(p1->pos, p2->pos);

	if (rij > truncLength)
		return 0; /* Too far away */
	
	return calcVCoulomb(rij) - calcVCoulomb(truncLength);
}

static double calcFCoulomb(double r)
{
	double couplingConstant = CHARGE_ELECTRON * CHARGE_ELECTRON
				/ (4 * M_PI * H2O_PERMETTIVITY);
	double k0 = calcInvDebyeLength();
	double exponentialPart = exp(-r * k0);
	
	return couplingConstant * exponentialPart * (k0 + 1/r) / r;
}
static void FCoulomb(Particle *p1, Particle *p2)
{	
	if (!ENABLE_COULOMB)
		return;

	if (!isChargedPair(p1->type, p2->type))
		return;

	double truncLen = config.truncationLen;
	double rij = nearestImageDistance(p2->pos, p1->pos);
	if (rij > truncLen)
		return; /* Too far away */
	
	Vec3 direction = nearestImageUnitVector(p1->pos, p2->pos);
	Vec3 forceVec;
	
	forceVec = scale(direction, calcFCoulomb(rij));
	p1->F = sub(p1->F, forceVec);
	p2->F = add(p2->F, forceVec);
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

		Fstack(&s->Bs[i], &s->Bs[i-1], 1);

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

		if (i < 2) continue;
		Fdihedral(&s->Ss[i], &s->Ps[i-1], &s->Ss[i-1], &s->Ps[i-2],
							DIHEDRAL_S3_P_5S3_P);
		Fstack(&s->Bs[i], &s->Bs[i-2], 2);
	}
}

static void pairForces(Particle *p1, Particle *p2)
{
	/* Nonbonded pair interactions are mutually exclusive. See Knotts. */
	if (isBondedBasePair(p1->type, p2->type))
		FbasePair(p1, p2);
	else if (isChargedPair(p1->type, p2->type))
		FCoulomb(p1, p2);
	else
		Fexclusion(p1, p2);
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
	assert(isSaneVector(p->vel));

	double *twiceK = (double*) data;
	*twiceK += p->m * length2(p->vel);
}
static double kineticEnergy(void)
{
	double twiceK = 0;
	forEveryParticleD(&kineticHelper, (void*) &twiceK);
	assert(isSaneNumber(twiceK));
	return twiceK/2;
}
double temperature(void)
{
	return 2.0 / (3.0 * BOLTZMANN_CONSTANT)
			* kineticEnergy() / numParticles();
}


typedef struct PotentialEnergies {
	double bond, angle, dihedral, stack, basePair, Coulomb, exclusion;
} PotentialEnergies;
static void pairPotentials(Particle *p1, Particle *p2, void *data)
{
	PotentialEnergies *pe = (PotentialEnergies*) data;
	
	pe->basePair += VbasePair(p1, p2);
	pe->Coulomb  += VCoulomb(p1, p2);
	pe->exclusion += Vexclusion(p1, p2);
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

		Vs += Vstack(&s->Bs[i], &s->Bs[i-1], 1);

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

		if (i < 2) continue;
		Vd += Vdihedral(&s->Ss[i], &s->Ps[i-1], &s->Ss[i-1], &s->Ps[i-2],
							DIHEDRAL_S3_P_5S3_P);
		Vs += Vstack(&s->Bs[i], &s->Bs[i-2], 2);
	}

	pe->bond     = Vb;
	pe->angle    = Va;
	pe->dihedral = Vd;
	pe->stack    = Vs;
}

/* Return energy stats of world, in electronvolts. */
static PotentialEnergies calcPotentialEnergies(void) {
	
	PotentialEnergies pe = {0, 0, 0, 0, 0, 0, 0};
	for (int s = 0; s < world.numStrands; s++)
		addPotentialEnergies(&world.strands[s], &pe);
	
	forEveryPairD(&pairPotentials, &pe);

	
	/* Convert to eV */
	pe.bond      *= ENERGY_FACTOR;
	pe.angle     *= ENERGY_FACTOR;
	pe.dihedral  *= ENERGY_FACTOR;
	pe.stack     *= ENERGY_FACTOR;
	pe.basePair  *= ENERGY_FACTOR;
	pe.Coulomb   *= ENERGY_FACTOR;
	pe.exclusion *= ENERGY_FACTOR;
	
	return pe;
}




/* ===== MOMENTUM FUNCTIONS ===== */

static void momentumHelper(Particle *p, void *data)
{
	Vec3 *Ptot = (Vec3*) data;
	*Ptot = add(*Ptot, scale(p->vel, p->m));
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
	double PPM = length(P) / numParticles();
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
	p->vel = scale(p->vel, lambda);
}
static void thermostat(void)
{
	if (config.thermostatTau <= 0)
		return;

	/* Mass and Boltzmann constant are 1 */ 
	double Tk  = temperature();
	assert(isSaneNumber(Tk));
	double T0  = config.thermostatTemp;
	double dt  = config.timeStep;
	double tau = config.thermostatTau;
	double lambda2 = 1 + dt/tau * (T0/Tk - 1);
	double lambda;
	if (lambda2 >= 0)
		lambda = sqrt(lambda2);
	else
		lambda = 1e-20; //TODO sane?

	forEveryParticleD(&thermostatHelper, (void*) &lambda);
}

static void verletHelper1(Particle *p)
{
	double dt = config.timeStep;

	assert(isSaneVector(p->pos));
	assert(isSaneVector(p->vel));
	assert(isSaneVector(p->F));

	/* vel(t + dt/2) = vel(t) + acc(t)*dt/2 */
	p->vel = add(p->vel, scale(p->F, dt / (2 * p->m)));

	/* pos(t + dt) = pos(t) + vel(t + dt/2)*dt */
	p->pos = add(p->pos, scale(p->vel, dt));

	assert(isSaneVector(p->pos));
	assert(isSaneVector(p->vel));
}
static void verletHelper2(Particle *p)
{
	double dt = config.timeStep;

	assert(isSaneVector(p->pos));
	assert(isSaneVector(p->vel));
	assert(isSaneVector(p->F));

	/* vel(t + dt) = vel(t + dt/2) + acc(t + dt)*dt/2 */
	p->vel = add(p->vel, scale(p->F, dt / (2*p->m)));

	assert(isSaneVector(p->vel));
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

	assert(isSaneVector(p->pos));
	assert(isSaneVector(p->vel));
	assert(isSaneVector(p->F));

	/* from v(t) to v(t + dt/2) */
	tmp1 = scale(p->vel, 1 - gamma*dt/2);
	tmp2 = scale(p->F, dt / (2*p->m)); /* p->F = regular + random force */
	p->vel = add(tmp1, tmp2);

	/* from r(t) to r(t + dt) */
	p->pos = add(p->pos, scale(p->vel, dt));

	assert(isSaneVector(p->pos));
	assert(isSaneVector(p->vel));
}
static void langevinBBKhelper2(Particle *p)
{
	double dt    = config.timeStep;
	double gamma = config.langevinGamma;
	double T     = config.thermostatTemp;

	assert(isSaneVector(p->pos));
	assert(isSaneVector(p->vel));
	assert(isSaneVector(p->F));

	/* Regular forces have been calculated. Add the random force due 
	 * to collisions to the total force. The result is:
	 * p->F = F(t + dt) + R(t + dt) */
	/* TODO check that compiler inlines this and precalculates the 
	 * prefactor before p->m when looping over all particles. */
	double Rstddev = sqrt(2 * BOLTZMANN_CONSTANT * T * gamma * p->m / dt);
	Vec3 R = randNormVec(Rstddev);
	p->F = add(p->F, R);

	/* from v(t + dt/2) to v(t + dt) */
	Vec3 tmp = add(p->vel, scale(p->F, dt / (2*p->m)));
	p->vel = scale(tmp, 1 / (1 + gamma*dt/2));

	assert(isSaneVector(p->pos));
	assert(isSaneVector(p->F));
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
	assert(worldSanityCheck());

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
	double E = K + pe.bond + pe.angle + pe.dihedral + pe.stack + pe.basePair + pe.Coulomb + pe.exclusion;

	printf("E = %e, K = %e, Vb = %e, Va = %e, Vd = %e, Vs = %e, Vbp = %e, Vpp = %e, Ve = %e, T = %f\n",
			E, K, pe.bond, pe.angle, pe.dihedral, pe.stack, pe.basePair, pe.Coulomb, pe.exclusion, T);
}



/* ===== MISC FUNCTIONS ===== */

Vec3 getCOM(Particle *ps, int num)
{
	Vec3 COM = {0, 0, 0};
	double M = 0; /* total mass */
	for (int i = 0; i < num; i++) {
		COM = add(COM, scale(ps[i].pos, ps[i].m));
		M += ps[i].m;
	}
	return scale(COM, 1/M);
}


