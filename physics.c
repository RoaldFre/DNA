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
	double d = nearestImageDistance(&p1->pos, &p2->pos) - d0;
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
	
	a = nearestImageVector(&p2->pos, &p1->pos);
	b = nearestImageVector(&p2->pos, &p3->pos);

	double dtheta = angle(&a, &b) - theta0;
	return ktheta/2 * dtheta*dtheta;
}
static void Fangle(Particle *p1, Particle *p2, Particle *p3, double theta0)
{
	if (!ENABLE_ANGLE)
		return;
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

/* DIHEDRAL */
static double Vdihedral(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
								double phi0)
{
	if (!ENABLE_DIHEDRAL)
		return 0;
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
	return distance2(&low, &high);
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
	double r2 = distance2(&p1->pos, &p2->pos);
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

	Vec3 drVec;
	sub(&p2->pos, &p1->pos, &drVec);
	double dr2 = length2(&drVec);

	assert(dr2 != 0);
	double dr4 = dr2 * dr2;
	double dr8 = dr4 * dr4;
	double dr14 = dr8 * dr4 * dr2;

	Vec3 F;
	scale(&drVec, -12 * kStack * (sigma12/dr14 - sigma6/dr8), &F);
	add(&p1->F, &F, &p1->F);
	sub(&p2->F, &F, &p2->F);
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
	
	double rsq = nearestImageDistance2(&p1->pos, &p2->pos);
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
	
	Vec3 rVec = nearestImageVector(&p1->pos, &p2->pos);
	double r = length(&rVec);
	if (r > config.truncationLen)
		return; /* Too far away */

	Vec3 direction;
	scale(&rVec, 1/r, &direction);

	double force = calcFbasePair(bpi, r);
	Vec3 forceVec;
	scale(&direction, force, &forceVec);

	/* Add force to particle objects */
	sub(&p1->F, &forceVec, &p1->F);
	add(&p2->F, &forceVec, &p2->F);
}

static void Fexclusion(Particle *p1, Particle *p2)
{
	if (!ENABLE_EXCLUSION)
		return;
	double sig;
	Vec3 forceVec;
	double force;
	
	double rij = nearestImageDistance(&p1->pos, &p2->pos);
	if (rij > D_CUT)
		return; /* Too far away */

	Vec3 direction = nearestImageUnitVector(&p1->pos, &p2->pos);
	
	/* Apply right parameters for types of molecules, 
	 * if bases and mismatched: SIGMA_0_CST*1.0
	 * if bases and matched: 0
	 * otherwise: SIGMA_0_CST*D_CUT */
	
	/* Lennard-Jones potential:
	 * potential = exCoupling ((sigma0 / r)^12 - (sigma0/r)^6) + epsilon.
	 * 
	 * For the force we differentiate with respect to r, so we get
	 * 4*EPSILON*[ - 12* sigma0^12 / r^13 + 6* sigma0^6/r^7 ] */

	if ((p1->type!=SUGAR) && (p1->type!=PHOSPHATE) && (p2->type!=SUGAR) 
						&& (p2->type!=PHOSPHATE)){
		
		if ( ((p1->type==BASE_A && p2->type==BASE_T) 
			|| (p1->type==BASE_T && p2->type==BASE_A))
			|| ((p1->type==BASE_G && p2->type==BASE_C) 
			|| (p1->type==BASE_C && p2->type==BASE_G)) )
			/* no force */
			return; 
		sig = SIGMA_0_CST*D_CUT_BASE;			
		
	} else {
		
		sig = SIGMA_0_CST*D_CUT;		
	}
	
	double sig2, sig4, sig6, sig12;
	double rInv, rInv2, rInv4, rInv7, rInv13;
	
	sig2 = sig*sig;
	sig4 = sig2*sig2;
	sig6 = sig2*sig4;
	sig12 = sig6*sig6;
	
	rInv = 1.0/rij;
	rInv2 = rInv*rInv;
	rInv4 = rInv2*rInv2;
	rInv7 = rInv4*rInv2*rInv;
	rInv13 = rInv7*rInv4*rInv4;

	force = 4*EPSILON*(-12*sig12*rInv13 + 6*sig6*rInv7); 
		
	/* Scale the direction with the calculated force */
	scale(&direction, force, &forceVec);

	/* Add force to particle objects */
	add(&p1->F, &forceVec, &p1->F);
	sub(&p2->F, &forceVec, &p2->F);
}

static double Vexclusion(Particle *p1, Particle *p2)
{
	if (!ENABLE_EXCLUSION)
		return 0;
	double sig;
	double potential;
	
	double rij = nearestImageDistance(&p1->pos, &p2->pos);
	if (rij > D_CUT)
		return 0; /* Too far away */

	
	/* Apply right parameters for types of molecules, 
	 * if bases and mismatched: SIGMA_0_CST*1.0
	 * if bases and matched: 0
	 * otherwise: SIGMA_0_CST*D_CUT */
	
	/* Lennard-Jones potential:
	 * potential = exCoupling ((sigma0 / r)^12 - (sigma0/r)^6) + epsilon.*/

	if ((p1->type!=SUGAR) && (p1->type!=PHOSPHATE) && (p2->type!=SUGAR) 
						&& (p2->type!=PHOSPHATE)){
		
		if ( ((p1->type==BASE_A && p2->type==BASE_T) 
			|| (p1->type==BASE_T && p2->type==BASE_A))
			|| ((p1->type==BASE_G && p2->type==BASE_C) 
			|| (p1->type==BASE_C && p2->type==BASE_G)) )
			/* no potential */
			return 0; 
		sig = SIGMA_0_CST*1.0;			
		
	} else {
		
		sig = SIGMA_0_CST*D_CUT;		
	}
	
	double sig2, sig4, sig6, sig12;
	double rInv, rInv2, rInv4, rInv6, rInv12;
	
	sig2 = sig*sig;
	sig4 = sig2*sig2;
	sig6 = sig2*sig4;
	sig12 = sig6*sig6;
	
	rInv = 1.0/rij;
	rInv2 = rInv*rInv;
	rInv4 = rInv2*rInv2;
	rInv6 = rInv4*rInv2;
	rInv12 = rInv6*rInv6;

	potential = 4*EPSILON*(sig12*rInv12 - sig6*rInv6) + EPSILON;
	
	/* correct so that at D_CUT, V zero */
	double dInv, dInv2, dInv4, dInv6, dInv12;
	dInv = 1.0/D_CUT;
	dInv2 = dInv*dInv;
	dInv4 = dInv2*dInv2;
	dInv6 = dInv4*dInv2;
	dInv12 = dInv6*dInv6;
	
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
static double VCoulomb(Particle *p1, Particle *p2)
{
	if (!ENABLE_COULOMB)
		return 0;

	if (p1->type != PHOSPHATE  ||  p2->type != PHOSPHATE)
		return 0; /* Only phosphates carry a charge */

	double truncLength = config.truncationLen;
	double rij = nearestImageDistance(&p1->pos, &p2->pos);

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

	if (p1->type != PHOSPHATE || p2->type != PHOSPHATE)
		return; /* Only phosphates carry a charge */

	double truncLen = config.truncationLen;
	double rij = nearestImageDistance(&p2->pos, &p1->pos);
	if (rij > truncLen)
		return; /* Too far away */
	
	Vec3 direction = nearestImageUnitVector(&p1->pos, &p2->pos);
	Vec3 forceVec;
	
	scale(&direction, calcFCoulomb(rij), &forceVec);
	sub(&p1->F, &forceVec, &p1->F);
	add(&p2->F, &forceVec, &p2->F);
}

#if 0
/* Fuzzy rope dynamics */
static void Frope(Particle *p1, Particle *p2, Particle *p3, Particle *p4)
{
	/* If the connections are on the same strand and on the same 
	 * particle: ignore. For instance something like this:
	 *
	 *   p1 (P[i])
	 *   |
	 *   |
	 * p2=p3 (S[i])
	 *   |
	 *   |
	 *   p4 (P[i+1])
	 *
	 * or:
	 *
	 *   p1 (P[i])
	 *   |
	 *   |
	 * p2=p4 (S[i]) --- p3 (B[i])
	 *
	 * etc...
	 *
	 * This gives a zero-distance because both connections meet at the 
	 * shared particle!
	 *
	 * TODO Currently brute forcing combinations. Smarter would be to 
	 * only test combinations that are possible (given the way how 
	 * forEveryConnection loops over particles). [Does that simplify 
	 * things?]
	 * */
	if (p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4)
		return;

	Vec3 *pos1 = &p1->pos;
	Vec3 *pos2 = &p3->pos;
	double truncLen = config.truncationLen;
	
	/* break if all points are further apart than 2*truncation length */
	double posDiff1 = nearestImageDistance(&p1->pos, &p3->pos);
	double posDiff2 = nearestImageDistance(&p1->pos, &p4->pos);
	double posDiff3 = nearestImageDistance(&p2->pos, &p3->pos);
	double posDiff4 = nearestImageDistance(&p2->pos, &p4->pos);
	
	if ( (posDiff1> (2*truncLen)) && (posDiff2 > (2*truncLen))
	&& (posDiff3 > (2*truncLen)) && (posDiff4 > (2*truncLen)) )
		return;
	
	Vec3 dir1 = nearestImageVector(&p1->pos, &p2->pos);
	Vec3 dir2 = nearestImageVector(&p3->pos, &p4->pos);
	
	/* calculate shortest distance between the two lines */
	double dist = nearestLineDistance(pos1, pos2, &dir1, &dir2);
	
	/* break if further apart than truncation length */
	if (dist > truncLen)
		return;

	//DEBUG XXX TODO
	if (dist < 1e-30)
		return;
		
	/* calculate the fuzzy rope force size */
	double coupling = ROPE_COUPLING;
	Vec3 forceVec;
	Vec3 direction;
	double rInv = ROPE_DIST / dist;
	double rInv2 = rInv * rInv;
	//double rInv4 = rInv2 * rInv2;
	//double rInv6 = rInv4 * rInv2;
	double force = coupling * rInv2;


	/* calculate and normalize the direction in which the force works */
	direction = cross(&dir1, &dir2);
	normalize(&direction, &direction);
	
	/* scale direction with calculated force */
	scale(&direction, force, &forceVec);

#if 0
	if (force > 1e-12)
		printf("%e\t%e\n",dist, force);
#endif

	assert(!(isnan(forceVec.x) || isnan(forceVec.y) || isnan(forceVec.z)));

	force = -force; //TODO this is correct sign?

	/* add force to particle objects (add for 1,2; sub for 3,4) */
	//TODO is this properly distributed?
	add(&p1->F, &forceVec, &p1->F);
	add(&p2->F, &forceVec, &p2->F);
	sub(&p3->F, &forceVec, &p3->F);
	sub(&p4->F, &forceVec, &p4->F);

	assert(isSaneVector(&p1->F));
	assert(isSaneVector(&p2->F));
	assert(isSaneVector(&p3->F));
	assert(isSaneVector(&p4->F));
}

static double Vrope(Particle *p1, Particle *p2, Particle *p3, Particle *p4)
{
	double posDiff1 = nearestImageDistance(&p1->pos, &p3->pos);
	double posDiff2 = nearestImageDistance(&p1->pos, &p4->pos);
	double posDiff3 = nearestImageDistance(&p2->pos, &p3->pos);
	double posDiff4 = nearestImageDistance(&p2->pos, &p4->pos);
	double truncLen = config.truncationLen;
	
	if ( (posDiff1> (2*truncLen)) && (posDiff2 > (2*truncLen))
	&& (posDiff3 > (2*truncLen)) && (posDiff4 > (2*truncLen)) )
		return 0;
		
	Vec3 *pos1 = &p1->pos;
	Vec3 *pos2 = &p3->pos;
	
	Vec3 dir1 = nearestImageVector(&p1->pos, &p2->pos);
	Vec3 dir2 = nearestImageVector(&p3->pos, &p4->pos);
	
	/* calculate shortest distance between the two (infinite) lines */
	double dist = nearestLineDistance(pos1, pos2, &dir1, &dir2);
	
	/* break if further apart than truncation length */
	if (dist > truncLen)
		return 0;
		
	/* calculate the fuzzy rope force size */
	double coupling = ROPE_COUPLING;
	double rInv = ROPE_DIST / dist;
	double rInv2 = rInv * rInv;
	double rInv3 = rInv2 * rInv;
	// double rInv4 = rInv2 * rInv2;
	// double rInv7 = rInv4 * rInv2 * rInv;

	double potential = 2 * coupling * rInv3;
	
	/* correct the potential to be zero at truncation length */
	double rInvTr = 1/truncLen;
	double rInvTr2 = rInvTr * rInvTr;
	double rInvTr3 = rInvTr2 * rInvTr;
	//double rInvTr4 = rInvTr2 * rInvTr2;
	//double rInvTr7 = rInvTr4 * rInvTr2 * rInvTr;
	double potentialCorr = 2 * coupling * rInvTr3;
	potential -= potentialCorr;
	
	return potential;
}
#endif

double nearestLineDistance(Vec3 *pos1, Vec3 *pos2, Vec3 *dist1, Vec3 *dist2)
{
	Vec3 rVec;
	double a11, a12,a22;
	double b1,b2;
	double D0;
	
	/* ugly hack */
	double eps = 0.00001e-40; /* Angstrom^4 */
	
	rVec = nearestImageVector(pos1, pos2);
		
	a11 = dot(dist1, dist1);
	a22 = dot(dist2, dist2);
	a12 = dot(dist1, dist2);
	b1 = dot(dist1, &rVec);
	b2 = dot(dist2, &rVec);
	
	D0 = a11 * a22 - a12 * a12;

#if 0
	//DEBUG
	printVector(dist1);
	printVector(dist2);
	printf("%e\n",D0);
#endif

	double sc, sN, sD = D0;
	double tc, tN, tD = D0;
	
	/* check endpoints */
	
	if (D0 < eps) {
		sN = 0.0;
		sD = 1.0;
		tN = b2;
		tD = a22;
	} else {
		sN = (a12 * b2 - a22 * b1);
		tN = (a11 * b2 - a12 * b1);
		
		if (sN < 0.0) {
			sN = 0.0;
			tN = b2;
			tD = a22;
		} else if (sN > sD) {
			sN = sD;
			tN = b2 + a12;
			tD = a22;
		}
	}
	
	if (tN < 0.0) {
		tN = 0.0;
 
		if (-b1 < 0.0)
			sN = 0.0;
		else if (-b1 > a11)
			sN = sD;
		else {
			sN = -b1;
			sD = a11;
		}
	}
	
	else if (tN > tD) {
		tN = tD;
		
		if ((-b1 + a12) < 0.0)
			sN = 0;
		else if ((-b1 + a12) > a11)
			sN = sD;
		else {
			sN = (-b1 + a12);
			sD = a11;
		}
	}
	
	if (fabs(sN) < eps)
		sc = 0.0;
	else
		sc = sN / sD;
	
	if (fabs(tN) < eps)
		tc = 0.0;
	else
		tc = tN / tD;
		
	Vec3 distVec1, distVec2, distVec3, distVecFinal;
	
	/* calculate distVecFinal = rVec + sc*dist1 - tc*dist2 */
	
	scale(dist1, sc, &distVec1);
	scale(dist2, tc, &distVec2);
	sub(&distVec1, &distVec2, &distVec3);
	add(&distVec3, &rVec, &distVecFinal);
	
	 
	return length(&distVecFinal);
	
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
	FbasePair(p1, p2);
	FCoulomb(p1, p2);
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

#if 0
	/* Connection-based forces */
	forEveryConnectionPair(&Frope);
#endif
}



/* ===== ENERGY FUNCTIONS ===== */

static void kineticHelper(Particle *p, void *data)
{
	assert(isSaneVector(&p->vel));

	double *twiceK = (double*) data;
	*twiceK += p->m * length2(&p->vel);
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
	double bond, angle, dihedral, stack, basePair, Coulomb, rope, exclusion;
} PotentialEnergies;
static void pairPotentials(Particle *p1, Particle *p2, void *data)
{
	PotentialEnergies *pe = (PotentialEnergies*) data;
	
	pe->basePair += VbasePair(p1, p2);
	pe->Coulomb  += VCoulomb(p1, p2);
	pe->exclusion += Vexclusion(p1, p2);
}

#if 0
static void addRopePotentialEnergy(Particle *p1, Particle *p2,
		Particle *p3, Particle *p4, void *data)
{
	double *V = (double*)data;
	*V += Vrope(p1, p2, p3, p4);
}
#endif

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
	
	PotentialEnergies pe = {0, 0, 0, 0, 0, 0, 0, 0};
	for (int s = 0; s < world.numStrands; s++)
		addPotentialEnergies(&world.strands[s], &pe);
	
	forEveryPairD(&pairPotentials, &pe);

#if 0
	forEveryConnectionPairD(&addRopePotentialEnergy, &pe.rope);
#endif
	
	/* Convert to eV */
	pe.bond      *= ENERGY_FACTOR;
	pe.angle     *= ENERGY_FACTOR;
	pe.dihedral  *= ENERGY_FACTOR;
	pe.stack     *= ENERGY_FACTOR;
	pe.basePair  *= ENERGY_FACTOR;
	pe.Coulomb   *= ENERGY_FACTOR;
#if 0
	pe.rope      *= ENERGY_FACTOR;
#endif
	pe.exclusion *= ENERGY_FACTOR;
	
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
	Vec3 tmp;

	assert(isSaneVector(&p->pos));
	assert(isSaneVector(&p->vel));
	assert(isSaneVector(&p->F));

	/* vel(t + dt/2) = vel(t) + acc(t)*dt/2 */
	scale(&p->F, dt / (2 * p->m), &tmp);
	add(&p->vel, &tmp, &p->vel);

	/* pos(t + dt) = pos(t) + vel(t + dt/2)*dt */
	scale(&p->vel, dt, &tmp);
	add(&p->pos, &tmp, &p->pos);

	assert(isSaneVector(&p->pos));
	assert(isSaneVector(&p->vel));
}
static void verletHelper2(Particle *p)
{
	double dt = config.timeStep;
	Vec3 tmp;

	assert(isSaneVector(&p->pos));
	assert(isSaneVector(&p->vel));
	assert(isSaneVector(&p->F));

	/* vel(t + dt) = vel(t + dt/2) + acc(t + dt)*dt/2 */
	scale(&p->F, dt / (2 * p->m), &tmp);
	add(&p->vel, &tmp, &p->vel);

	assert(isSaneVector(&p->vel));
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

	assert(isSaneVector(&p->pos));
	assert(isSaneVector(&p->vel));
	assert(isSaneVector(&p->F));

	/* from v(t) to v(t + dt/2) */
	scale(&p->vel, 1 - gamma*dt/2, &tmp1);

	/* p->F is regular force + random collision force */
	scale(&p->F, dt / (2 * p->m), &tmp2);

	add(&tmp1, &tmp2, &p->vel);

	/* from r(t) to r(t + dt) */
	scale(&p->vel, dt, &tmp1);
	add(&p->pos, &tmp1, &p->pos);

	assert(isSaneVector(&p->pos));
	assert(isSaneVector(&p->vel));
}
static void langevinBBKhelper2(Particle *p)
{
	double dt    = config.timeStep;
	double gamma = config.langevinGamma;
	double T     = config.thermostatTemp;

	assert(isSaneVector(&p->pos));
	assert(isSaneVector(&p->vel));
	assert(isSaneVector(&p->F));

	/* Regular forces have been calculated. Add the random force due 
	 * to collisions to the total force. The result is:
	 * p->F = F(t + dt) + R(t + dt) */
	/* TODO check that compiler inlines this and precalculates the 
	 * prefactor before p->m when looping over all particles. */
	double Rstddev = sqrt(2 * BOLTZMANN_CONSTANT * T * gamma * p->m / dt);
	Vec3 R = randNormVec(Rstddev);
	add(&p->F, &R, &p->F);

	/* from v(t + dt/2) to v(t + dt) */
	Vec3 tmp;
	scale(&p->F, dt / (2 * p->m), &tmp);
	add(&p->vel, &tmp, &tmp);
	scale(&tmp, 1 / (1 + gamma*dt/2), &p->vel);

	assert(isSaneVector(&p->pos));
	assert(isSaneVector(&p->F));
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
	double E = K + pe.bond + pe.angle + pe.dihedral + pe.stack + pe.basePair + pe.Coulomb + pe.rope + pe.exclusion;

	printf("E = %e, K = %e, Vb = %e, Va = %e, Vd = %e, Vs = %e, Vbp = %e, Vpp = %e, Vr = %e, Ve = %e, T = %f\n",
			E, K, pe.bond, pe.angle, pe.dihedral, pe.stack, pe.basePair, pe.Coulomb, pe.rope, pe.exclusion, T);
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


