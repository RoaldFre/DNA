#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "physics.h"
#include "integrator.h"
#include "world.h"
#include "spgrid.h"

static InteractionSettings interactions;
static double truncationLenSq; /* Cached: interactions.truncationLen^2 */

void registerInteractionSettings(InteractionSettings interactionSettings)
{
	interactions = interactionSettings;
	truncationLenSq = SQUARE(interactionSettings.truncationLen);
}



/* ===== FORCES AND POTENTIALS ===== */

/* Generic Lennard-Jones with potential depth epsilon.
 * args:
 *   rEq2 = (equilibrium distance r_eq)^2
 *   r2   = (actual distance r)^2
 * result:
 *   V = epsilon * [ (r_eq / r)^12  -  2*(r_eq / r)^6] */
static double calcVLJ(double epsilon, double rEq2, double r2)
{
	double rfrac2  = rEq2 / r2;
	double rfrac4  = SQUARE(rfrac2);
	double rfrac6  = rfrac4 * rfrac2;
	double rfrac12 = SQUARE(rfrac6);
	return epsilon * (rfrac12 - 2*rfrac6);
}
/* Generic Lennard-Jones with potential depth epsilon. Returns the Force 
 * *dividid by displacement distance r*!! So use this to scale the 
 * displacement vector to get the correct force!
 * args:
 *   rEq2 = (equilibrium distance r_eq)^2
 *   r2   = (actual distance r)^2
 * result:
 *   F/r = 12 * epsilon * [ r_eq^6 / r^8  -  r_eq^12 / r^14 ] */
static double calcFLJperDistance(double epsilon, double rEq2, double r2)
{
	double rfrac2  = rEq2 / r2;
	double rfrac4  = SQUARE(rfrac2);
	double rfrac6  = rfrac4 * rfrac2;
	double rfrac12 = SQUARE(rfrac6);
	return 12 * epsilon * (rfrac6 - rfrac12) / r2;
}


/* BOND */
/* V = k1 * (dr - d0)^2  +  k2 * (d - d0)^4
 * where dr is the distance between the particles */
static double Vbond(Particle *p1, Particle *p2, double d0)
{
	if (!interactions.enableBond)
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
	if (!interactions.enableBond)
		return;
	double k1 = BOND_K1;
	double k2 = BOND_K2;
	Vec3 drVec = nearestImageVector(p1->pos, p2->pos);
	double dr = length(drVec);
	double d  = dr - d0;
	double d3 = d * d * d;

	Vec3 drVecNormalized = scale(drVec, 1/dr);
	Vec3 F = scale(drVecNormalized, 2*k1*d + 4*k2*d3);

	debugVectorSanity(F, "Fbond");

	p1->F = add(p1->F, F);
	p2->F = sub(p2->F, F);
}
static double getSugarBaseBondDistance(ParticleType base)
{
	switch (base) {
		case BASE_A: return BOND_S_A;
		case BASE_T: return BOND_S_T;
		case BASE_C: return BOND_S_C;
		case BASE_G: return BOND_S_G;
		default: assert(false); return 0;
	}
}
/* Bond between a sugar and a base */
static double VbondSB(Particle *sugar, Particle *base)
{
	assert(sugar->type == SUGAR);
	return Vbond(sugar, base, getSugarBaseBondDistance(base->type));
}
static void FbondSB(Particle *sugar, Particle *base)
{
	assert(sugar->type == SUGAR);
	Fbond(sugar, base, getSugarBaseBondDistance(base->type));
}

/* ANGLE */
/* V = ktheta/2 * (theta - theta0)^2
 *
 * p1 \       /p3
 *     \theta/
 *      \   /
 *       \ /
 *        p2
 */
typedef struct {
	double P5SB;
	double P3SB;
} AngleBaseInfo;
static AngleBaseInfo getAngleBaseInfo(ParticleType base)
{
	AngleBaseInfo info;
	switch (base) {
	case BASE_A:
		info.P5SB = ANGLE_P_5S_A;
		info.P3SB = ANGLE_P_3S_A;
		break;
	case BASE_T:
		info.P5SB = ANGLE_P_5S_T;
		info.P3SB = ANGLE_P_3S_T;
		break;
	case BASE_C:
		info.P5SB = ANGLE_P_5S_C;
		info.P3SB = ANGLE_P_3S_C;
		break;
	case BASE_G:
		info.P5SB = ANGLE_P_5S_G;
		info.P3SB = ANGLE_P_3S_G;
		break;
	default:
		fprintf(stderr, "Unknown base type in getAngleBaseInfo!\n");
		assert(false);
	}
	return info;
}

static double Vangle(Particle *p1, Particle *p2, Particle *p3, double theta0)
{
	if (!interactions.enableAngle)
		return 0;
	Vec3 a, b;
	double ktheta = ANGLE_COUPLING;
	
	a = nearestImageVector(p2->pos, p1->pos);
	b = nearestImageVector(p2->pos, p3->pos);

	double dtheta = angle(a, b) - theta0;
	return ktheta/2 * dtheta*dtheta;
}
static double VangleP5SB(Particle *p, Particle *s, Particle *b)
{
	AngleBaseInfo info = getAngleBaseInfo(b->type);
	return Vangle(p, s, b, info.P5SB);
}
static double VangleP3SB(Particle *p, Particle *s, Particle *b)
{
	AngleBaseInfo info = getAngleBaseInfo(b->type);
	return Vangle(p, s, b, info.P3SB);
}
static void Fangle(Particle *p1, Particle *p2, Particle *p3, double theta0)
{
	if (!interactions.enableAngle)
		return;
	Vec3 a, b;
	double ktheta = ANGLE_COUPLING;
	
	a = nearestImageVector(p2->pos, p1->pos);
	b = nearestImageVector(p2->pos, p3->pos);
	
	double lal = length(a);
	double lbl = length(b);
	double adotb = dot(a, b);
	double costheta = adotb / (lal * lbl);

	if (UNLIKELY(fabs(costheta) >= 1 - 1e-5))
		/* Note, sometimes |costheta| > 1, due to numerical errors!
		 * Either way, if |costheta| almost equal to 1: theta is 
		 * almost equal to 0 or pi and we are at an (unstable) 
		 * equilibrium. We just bail out without applying a force, 
		 * because otherwise we would end up dividing by zero 
		 * below. */
		return;

	double theta = acos(costheta);
	double sintheta = sqrt(1 - costheta*costheta);
	/* Note: the positive sign for the square root is always correct, 
	 * because the angle is always less than 180 degrees! */

	Vec3 tmp1, tmp2, F1, F2, F3;

	tmp1 = scale(b, 1/(lal * lbl));
	tmp2 = scale(a, adotb / (lal*lal*lal * lbl));
	F1 = scale(sub(tmp1, tmp2), ktheta * (theta - theta0) / sintheta);
	debugVectorSanity(F1, "Fangle F1");
	p1->F = add(p1->F, F1);	

	tmp1 = scale(a, 1/(lal * lbl));
	tmp2 = scale(b, adotb / (lbl*lbl*lbl * lal));
	F3 = scale(sub(tmp1, tmp2), ktheta * (theta - theta0) / sintheta);
	debugVectorSanity(F3, "Fangle F3");
	p3->F = add(p3->F, F3);	

	F2 = add(F1, F3); /* actually -F2 */
	debugVectorSanity(F2, "Fangle F2");
	p2->F = sub(p2->F, F2);	

	assert(fabs(dot(a, F1) / length(a) / length(F1)) < 1e-5);
	assert(fabs(dot(b, F3) / length(b) / length(F3)) < 1e-5);
}
static void FangleP5SB(Particle *p, Particle *s, Particle *b)
{
	AngleBaseInfo info = getAngleBaseInfo(b->type);
	Fangle(p, s, b, info.P5SB);
}
static void FangleP3SB(Particle *p, Particle *s, Particle *b)
{
	AngleBaseInfo info = getAngleBaseInfo(b->type);
	Fangle(p, s, b, info.P3SB);
}


/* DIHEDRAL */
typedef struct {
	double sinDihedral;
	double cosDihedral;
} DihedralCache;

typedef struct {
	DihedralCache BS3P5S;
	DihedralCache S3P5SB;
} DihedralBaseInfo;

static DihedralBaseInfo dihedralsBases[NUM_BASE_TYPES];
static DihedralCache dihedralP5S3P5S;
static DihedralCache dihedralS3P5S3P;

static DihedralCache makeDihedralCache(double dihedralAngle)
{
	DihedralCache ret;
	ret.sinDihedral = sin(dihedralAngle);
	ret.cosDihedral = cos(dihedralAngle);
	return ret;
}
static void initDihedralCache(void)
{
	dihedralsBases[BASE_A].BS3P5S = makeDihedralCache(DIHEDRAL_A_S3_P_5S);
	dihedralsBases[BASE_A].S3P5SB = makeDihedralCache(DIHEDRAL_S3_P_5S_A);

	dihedralsBases[BASE_T].BS3P5S = makeDihedralCache(DIHEDRAL_T_S3_P_5S);
	dihedralsBases[BASE_T].S3P5SB = makeDihedralCache(DIHEDRAL_S3_P_5S_T);

	dihedralsBases[BASE_C].BS3P5S = makeDihedralCache(DIHEDRAL_C_S3_P_5S);
	dihedralsBases[BASE_C].S3P5SB = makeDihedralCache(DIHEDRAL_S3_P_5S_C);

	dihedralsBases[BASE_G].BS3P5S = makeDihedralCache(DIHEDRAL_G_S3_P_5S);
	dihedralsBases[BASE_G].S3P5SB = makeDihedralCache(DIHEDRAL_S3_P_5S_G);

	dihedralP5S3P5S = makeDihedralCache(DIHEDRAL_P_5S3_P_5S);
	dihedralS3P5S3P = makeDihedralCache(DIHEDRAL_S3_P_5S3_P);
}

/* V = k * (1 - cos(phi - phi0)) */
static double Vdihedral(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
							DihedralCache phi0)
{
	if (!interactions.enableDihedral)
		return 0;

	Vec3 r1 = nearestImageVector(p1->pos, p2->pos);
	Vec3 r2 = nearestImageVector(p2->pos, p3->pos);
	Vec3 r3 = nearestImageVector(p3->pos, p4->pos);

	/* Trigonometric magic:
	 * We need to compute cos(phi - phi0).
	 * We can easily compute cos(phi) and sin(phi) and we have:
	 *    cos(phi - phi0) = cos(phi)cos(phi0) + sin(phi)sin(phi0)
	 */
	double sinPhi, cosPhi;
	sinCosDihedral(r1, r2, r3, &sinPhi, &cosPhi);

	double sinPhi0 = phi0.sinDihedral;
	double cosPhi0 = phi0.cosDihedral;

	double cosPhiPhi0 = cosPhi*cosPhi0 + sinPhi*sinPhi0;
	return DIHEDRAL_COUPLING * (1 - cosPhiPhi0);
}
static double VdihedralBS3P5S(Particle *b, Particle *s1,
				Particle *p, Particle *s2)
{
	assert(0 <= b->type && b->type < 4);
	return Vdihedral(b, s1, p, s2, dihedralsBases[b->type].BS3P5S);
}
static double VdihedralS3P5SB(Particle *s1, Particle *p,
				Particle *s2, Particle *b)
{
	assert(0 <= b->type && b->type < 4);
	return Vdihedral(s1, p, s2, b, dihedralsBases[b->type].S3P5SB);
}
/* Return the dihedral force to the target particle. */
static Vec3 FdihedralNumDiffParticle(Particle *target, 
		Particle *p1, Particle *p2, Particle *p3, Particle *p4, 
		double Vorig, DihedralCache phi0)
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

	debugVectorSanity(F, "FdihedralNumDiffParticle");

	return F;
}
static void FdihedralNumericalDiff(Particle *p1, Particle *p2,
			Particle *p3, Particle *p4, DihedralCache phi0)
{
	if (!interactions.enableDihedral)
		return;

	/* This is a *mess* to do analytically, so we do a numerical 
	 * differentiation instead. */
	double Vorig = Vdihedral(p1, p2, p3, p4, phi0);
	Vec3 F1 = FdihedralNumDiffParticle(p1, p1, p2, p3, p4, Vorig, phi0);
	Vec3 F2 = FdihedralNumDiffParticle(p2, p1, p2, p3, p4, Vorig, phi0);
	Vec3 F3 = FdihedralNumDiffParticle(p3, p1, p2, p3, p4, Vorig, phi0);
	Vec3 negF4 = add(F1, add(F2, F3)); /* -F4 = F1 + F2 + F3 */
	p1->F = add(p1->F, F1);
	p2->F = add(p2->F, F2);
	p3->F = add(p3->F, F3);
	p4->F = sub(p4->F, negF4);
}

/*
 * Returns the *transposed* Jacobian matrix of the matrix cross product
 *    r12 x r23  =  (r2 - r1) x (r3 - r2)
 * where r1, r2 and r3 are vectors and differentiation is done with regards 
 * to ri (with i the given parameter) with all other rj (j != i) assumed to 
 * be fixed and independent of ri.
 *
 * The math is worked out in:
 * http://www-personal.umich.edu/~riboch/pubfiles/riboch-JacobianofCrossProduct.pdf
 * The result:
 *   J_ri(r12 x r23) = r12^x * J_ri(r23) - r23^x * J_ri(r12)
 * where r12^x and r23^x are the matrix form of the cross product:
 *         ( 0   -Az   Ay)
 *   A^x = ( Az   0   -Ax)
 *         (-Ay   Ax   0 )
 * and J_ri(r12) and J_ri(r23) are the jacobi matrices for r12 and r23, 
 * which are either 0, or +/- the identity matrix.
 * We have:
 *   J_r1(r12) = -1,    J_r2(r12) = +1,    J_r3(r12) =  0,
 *   J_r1(r23) =  0,    J_r2(r23) = -1,    J_r3(r23) = +1.
 */
static __inline__ Mat3 crossProdTransposedJacobian(Vec3 r12, Vec3 r23, int i)
{
	/* *Transposed* matrix form of the cross product. */
	Mat3 r12matrCrossProd = mat3(   0,    r12.z, -r12.y,
	                             -r12.z,    0,    r12.x,
	                              r12.y, -r12.x,    0   );
	/* *Transposed* matrix form of the cross product. */
	Mat3 r23matrCrossProd = mat3(   0,    r23.z, -r23.y,
	                             -r23.z,    0,    r23.x,
	                              r23.y, -r23.x,    0   );

	switch (i) {
	case 1:
		return r23matrCrossProd;
	case 2:
		return matScale(matAdd(r12matrCrossProd, r23matrCrossProd), -1);
	case 3:
		return r12matrCrossProd;
	default:
		assert(false); die("Internal error!\n");
		return mat3(0,0,0,0,0,0,0,0,0); /* To make compiler happy. */
	}
}

/* Derivative of the dihedral potential. This shit gets quite involved.
 *
 * We want to compute
 *   Fi = -grad_ri V(r1, r2, r3, r4)
 * Fi the force on the i'th particle and rj the position of the j'th 
 * particle.
 * 
 * We can rather easily write V as a function of
 *   A = r12 x r23  =  (r2 - r1) x (r3 - r2)
 * and
 *   B = r23 x r34  =  (r3 - r2) x (r4 - r3)
 *
 * Indeed, we have:
 *   V = DIHEDRAL_COUPLING * (1 - cos(phi - phi0))
 * where
 *   phi = acos(A dot B / |A| * |B|)
 *
 * We can then write Fi as:
 *   Fi = -grad_ri Vi(A(r1, r2, r3), B(r2, r3, r4))
 * and use the appropriate chain rule:
 *   Fi = - [J_ri(A)]^t * (grad_A V)
 *        - [J_ri(B)]^t * (grad_B V)
 * where [J_ri(X)]^t is the transpose of the Jacobian matrix of the vector 
 * function X with regards to ri.
 *
 * After some calculus and algebra, we get
 *   grad_A V = DIHEDRAL_COUPLING
 *                * [ sin(phi0)/sin(phi) * (A dot B)/(|A|^2 |B|^2)
 *                                  - cos(phi0) / (|A||B|) ]
 *                * [B - A (A dot B) / (|B||A|)]
 * and due to the symmetry in V of A and B, grad_B V is exactly the same, 
 * but with A and B interchanged.
 *   grad_B V = DIHEDRAL_COUPLING
 *                * [ sin(phi0)/sin(phi) * (A dot B)/(|A|^2 |B|^2)
 *                                  - cos(phi0) / (|A||B|) ]
 *                * [A - B (A dot B) / (|B||A|)]
 *
 *
 * NOTE: when debugging the dihedral force for energy conservation, you 
 * also need to enable angle and bond interactions (which should be 
 * debugged first), otherwise you get a badly conditioned problem and you 
 * will experience blow up.
 */
static void Fdihedral(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
							DihedralCache phi0)
{
	if (!interactions.enableDihedral)
		return;

	Vec3 r12 = nearestImageVector(p1->pos, p2->pos);
	Vec3 r23 = nearestImageVector(p2->pos, p3->pos);
	Vec3 r34 = nearestImageVector(p3->pos, p4->pos);

	double sinPhi, cosPhi; //TODO find better algorithm to only compute sinPhi
	sinCosDihedral(r12, r23, r34, &sinPhi, &cosPhi);
	if (UNLIKELY(fabs(sinPhi) < 1e-5)) //TODO just check for ==0? -> saves an fabs!
		return; /* (Unstable) equilibrium. */

	double sinPhi0 = phi0.sinDihedral;
	double cosPhi0 = phi0.cosDihedral;

	Vec3 A = cross(r12, r23);
	Vec3 B = cross(r23, r34);

	double lAl2 = length2(A);
	double lBl2 = length2(B);
	double lAl2lBl2 = lAl2 * lBl2;
	double lAllBl = sqrt(lAl2lBl2);
	double AdB = dot(A, B);
	double negGradPrefactor = DIHEDRAL_COUPLING * (cosPhi0 / lAllBl
					- sinPhi0/sinPhi * AdB/lAl2lBl2);
	/* -grad_A(V) */
	Vec3 negGrad_AV = scale(
			sub(B, scale(A, AdB/lAl2)),
			negGradPrefactor);
	/* -grad_B(V) */
	Vec3 negGrad_BV = scale(
			sub(A, scale(B, AdB/lBl2)),
			negGradPrefactor);

	/* F1 */
	Mat3 J_r1A = crossProdTransposedJacobian(r12, r23, 1); /* [J_r1(A)]^t */
	/* [J_r1(B)]^t = 0*/
	Vec3 F1 = matApply(J_r1A, negGrad_AV);

	/* F2 */
	Mat3 J_r2A = crossProdTransposedJacobian(r12, r23, 2); /* [J_r2(A)]^t */
	Mat3 J_r2B = crossProdTransposedJacobian(r23, r34, 1); /* [J_r2(B)]^t */
	Vec3 F2 = add(matApply(J_r2A, negGrad_AV), matApply(J_r2B, negGrad_BV));

	/* F4 (not F3 because that has two non-zero jacobi matrices, 
	 * whereas F4 has only one non-zero jacobi matrix, so it's 
	 * faster to compute)*/
	/* [J_r4(A)]^t = 0*/
	Mat3 J_r4B = crossProdTransposedJacobian(r23, r34, 3); /* [J_r4(B)]^t */
	Vec3 F4 = matApply(J_r4B, negGrad_BV);

	/* -F3 */
	Vec3 negF3 = add(add(F1, F2), F4);

	p1->F = add(p1->F, F1);
	p2->F = add(p2->F, F2);
	p3->F = sub(p3->F, negF3);
	p4->F = add(p4->F, F4);


#if 0
#ifdef DEBUG
	/* NOTE: lots of false positives if one component of the force is 
	 * small! */

	double eps = 0.01;

	/* Explicitly calculate F3 as well. */
	/* F3 */
	Mat3 J_r3A = crossProdTransposedJacobian(r12, r23, 3); /* [J_r3(A)]^t */
	Mat3 J_r3B = crossProdTransposedJacobian(r23, r34, 2); /* [J_r3(B)]^t */
	Vec3 F3 = add(matApply(J_r3A, negGrad_AV), matApply(J_r3B, negGrad_BV));
	assertVecEqualsEpsilon(scale(negF3, -1), F3, eps);

	/* Check with numerical differentiation code */
	double Vorig = Vdihedral(p1, p2, p3, p4, phi0);
	Vec3 correctF1 = FdihedralNumDiffParticle(p1, p1, p2, p3, p4, Vorig, phi0);
	Vec3 correctF2 = FdihedralNumDiffParticle(p2, p1, p2, p3, p4, Vorig, phi0);
	Vec3 correctF3 = FdihedralNumDiffParticle(p3, p1, p2, p3, p4, Vorig, phi0);
	Vec3 correctF4 = FdihedralNumDiffParticle(p4, p1, p2, p3, p4, Vorig, phi0);
	assertVecEqualsEpsilon(F1, correctF1, eps);
	assertVecEqualsEpsilon(F2, correctF2, eps);
	assertVecEqualsEpsilon(F3, correctF3, eps);
	assertVecEqualsEpsilon(F4, correctF4, eps);
#endif
#endif
}




static void FdihedralBS3P5S(Particle *b, Particle *s1,
				Particle *p, Particle *s2)
{
	Fdihedral(b, s1, p, s2, dihedralsBases[b->type].BS3P5S);
}
static void FdihedralS3P5SB(Particle *s1, Particle *p,
				Particle *s2, Particle *b)
{
	Fdihedral(s1, p, s2, b, dihedralsBases[b->type].S3P5SB);
}


/* STACKING */
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
 * Particle 1 is the 'bottom' particle in the helix (lower index in the 
 * strand), particle 2 is the 'top' particle of the helix (higher index in 
 * the strand).
 *
 * TODO WARNING, complementary strand should be of the correct 3'-5' 
 * "direction" as the "normal" strand, or this breaks down.
 *
 * monomerDistance:
 *   1 for immediate neighbours i and i+1
 *   2 for   next    neighbours i and i+2
 *   n for           neighbours i and i+n */
static double neighbourStackDistance2(ParticleType t1, ParticleType t2, int monomerDistance)
{
	HelixInfo hi1 = getHelixInfo(t1);
	HelixInfo hi2 = getHelixInfo(t2);
	return helixDistance2(hi2, hi1, /* TODO order of hi1 and hi2 like 
					   this should be wrong, but 
					   works... */
			monomerDistance * HELIX_DELTA_Z,
			monomerDistance * HELIX_DELTA_PHI);
}
/* monomerDistance:
 *   1 for immediate neighbours (i and i+1)
 *   2 for   next    neighbours (i and i+2) */
static double Vstack(Particle *p1, Particle *p2, int monomerDistance)
{
	assert(isBase(p1->type) && isBase(p2->type));

	if (!interactions.enableStack)
		return 0;

	double rSq = nearestImageDistance2(p1->pos, p2->pos);
	if (rSq > truncationLenSq)
		return 0;

	double rEqSq = neighbourStackDistance2(p1->type, p2->type,
						monomerDistance);
	double V = calcVLJ(STACK_COUPLING, rEqSq, rSq)
			- calcVLJ(STACK_COUPLING, rEqSq, truncationLenSq); //TODO cache correction!
	return V;
}
static void Fstack(Particle *p1, Particle *p2, int monomerDistance)
{
	assert(isBase(p1->type) && isBase(p2->type));

	if (!interactions.enableStack)
		return;

	Vec3 r = nearestImageVector(p1->pos, p2->pos);
	double rSq = length2(r);
	if (rSq > truncationLenSq)
		return;

	double rEqSq = neighbourStackDistance2(p1->type, p2->type,
						monomerDistance);
	double FperDist = calcFLJperDistance(STACK_COUPLING, rEqSq, rSq);
	Vec3 F = scale(r, FperDist);
	debugVectorSanity(F, "Fstack");
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
static BasePairInfo getBasePairInfo(Particle *p1, Particle *p2)
{
	ParticleType t1 = p1->type;
	ParticleType t2 = p2->type;

	BasePairInfo bpi;
	bpi.coupling = -1;
	bpi.distance2 = -1;

	switch (interactions.basePairInteraction) {
	case BASE_PAIR_ALL:
		break;
	case BASE_PAIR_HAIRPIN:
		if (p1->strand != p2->strand)
			return bpi;

		int n  = p1->strand->numMonomers;
		int i1 = p1->strandIndex;
		int i2 = p2->strandIndex;

		if (i1 != n - 1 - i2)
			return bpi;
		break;
	case BASE_PAIR_DOUBLE_STRAND:
		if (p1->strand == p2->strand)
			return bpi;
		if (p1->strandIndex != p2->strandIndex)
			return bpi;
		break;
	default:
		assert(false);
		return bpi;
	}

	if ((t1 == BASE_A  &&  t2 == BASE_T) 
			||  (t1 == BASE_T  &&  t2 == BASE_A)) {
		bpi.coupling = BASE_PAIR_COUPLING_A_T;
		bpi.distance2 = SQUARE(BASE_PAIR_DISTANCE_A_T);
	} else if ((t1 == BASE_C  &&  t2 == BASE_G)
			||  (t1 == BASE_G  &&  t2 == BASE_C)) {
		bpi.coupling = BASE_PAIR_COUPLING_G_C;
		bpi.distance2 = SQUARE(BASE_PAIR_DISTANCE_G_C);
	}
	return bpi;
}
/* Returns true for pairs of particles that can form a base pair binding. */
static bool isBondedBasePair(Particle *p1, Particle *p2)
{
	BasePairInfo bpi = getBasePairInfo(p1, p2);
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
	if (!interactions.enableBasePair)
		return 0;

	BasePairInfo bpi = getBasePairInfo(p1, p2);
	if (bpi.coupling < 0)
		return 0; /* Wrong pair */
	
	double rsq = nearestImageDistance2(p1->pos, p2->pos);
	if (rsq > truncationLenSq)
		return 0; /* Too far away */

	return calcVbasePair(bpi, rsq) - calcVbasePair(bpi, truncationLenSq);
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
	if (!interactions.enableBasePair)
		return;

	BasePairInfo bpi = getBasePairInfo(p1, p2);
	if (bpi.coupling < 0)
		return; /* Wrong pair */
	
	Vec3 rVec = nearestImageVector(p1->pos, p2->pos);
	double r = length(rVec);
	if (r > interactions.truncationLen)
		return; /* Too far away */

	Vec3 direction = scale(rVec, 1/r);

	double force = calcFbasePair(bpi, r);
	Vec3 forceVec = scale(direction, force);
	debugVectorSanity(forceVec, "FbasePair");

	/* Add force to particle objects */
	p1->F = sub(p1->F, forceVec);
	p2->F = add(p2->F, forceVec);
}


/* EXCLUSION */
/* Two particles feel exclusion forces if they are not connected by a 
 * direct bond. */
static bool feelExclusion(Particle *p1, Particle *p2)
{
	if (p1->strand != p2->strand)
		return true;

	if (ABS(p1->strandIndex - p2->strandIndex) > 1)
		return true;

	ParticleType t1 = p1->type;
	ParticleType t2 = p2->type;
	int i1 = p1->strandIndex;
	int i2 = p2->strandIndex;

	/* Use the ordering of the type enum:
	 * PHOSPHATE, SUGAR, BASE_X */
	if (t1 > t2)
		return feelExclusion(p2, p1);

	switch (t1) {
	case PHOSPHATE: /* t2 is PHOSPHATE, SUGAR or BASE */
		return t2 == PHOSPHATE || t2 != SUGAR || (i1 == i2 + 1);
	case SUGAR: /* t2 is a SUGAR or BASE */
		return t2 == SUGAR || i1 != i2;
	default: /* t1 and t2 are BASEs */
		if (interactions.mutuallyExclusivePairForces)
			assert(i1 != i2);
		return i1 != i2;
	}
}

/* Return (the exclusion cut off distance)^2. This is the distance where 
 * the repulsive part of the Lennard Jones potential stops. */
static double getExclusionCutOff2(ParticleType t1, ParticleType t2){
	if (isBase(t1) && isBase(t2))
		return SQUARE(EXCLUSION_DISTANCE_BASE);
	else
		return SQUARE(EXCLUSION_DISTANCE);
}
/* Only the repulsive part of a Lennard-Jones potential. (The potential 
 * gets lifted so it's zero at infinity.) */
static double Vexclusion(Particle *p1, Particle *p2)
{
	if (!interactions.enableExclusion)
		return 0;

	if (!feelExclusion(p1, p2))
		return 0;

	double cutOffSq = getExclusionCutOff2(p1->type, p2->type);
	double rSq = nearestImageDistance2(p1->pos, p2->pos);
	if (rSq > cutOffSq)
		return 0;

	return calcVLJ(EXCLUSION_COUPLING, cutOffSq, rSq) + EXCLUSION_COUPLING;
}
static void Fexclusion(Particle *p1, Particle *p2)
{
	if (!interactions.enableExclusion)
		return;

	if (!feelExclusion(p1, p2))
		return;

	double cutOffSq = getExclusionCutOff2(p1->type, p2->type);
	Vec3 r = nearestImageVector(p1->pos, p2->pos);
	double rSq = length2(r);
	if (rSq > cutOffSq)
		return;

	double FperDist = calcFLJperDistance(EXCLUSION_COUPLING, cutOffSq, rSq);
	Vec3 F = scale(r, FperDist);
	debugVectorSanity(F, "Fexclusion");
	p1->F = add(p1->F, F);
	p2->F = sub(p2->F, F);
}


/* COULOMB */
static double calcInvDebyeLength(void)
{
	double T = getHeatBathTemperature();
	double saltCon = interactions.saltConcentration;
	double lambdaBDenom, lambdaB;

	if (T == 0)
		return 1e100;

	lambdaBDenom = 4 * M_PI * H2O_PERMETTIVITY
			* BOLTZMANN_CONSTANT * T;
	lambdaB = SQUARE(ELECTRON_CHARGE) / lambdaBDenom;
	double invDebLength = sqrt(8 * M_PI * lambdaB * AVOGADRO * saltCon);
	return invDebLength;
}
static double calcVCoulomb(double r)
{
	double couplingConstant = SQUARE(ELECTRON_CHARGE)
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
	if (!interactions.enableCoulomb)
		return 0;

	if (!isChargedPair(p1->type, p2->type))
		return 0;

	double truncLength = interactions.truncationLen;
	double r = nearestImageDistance(p1->pos, p2->pos);

	if (r > truncLength)
		return 0; /* Too far away */
	
	return calcVCoulomb(r) - calcVCoulomb(truncLength);
}

static double calcFCoulomb(double r)
{
	double couplingConstant = SQUARE(ELECTRON_CHARGE)
				/ (4 * M_PI * H2O_PERMETTIVITY);
	double k0 = calcInvDebyeLength(); //TODO cache per iteration (or whenever T changes)
	double exponentialPart = exp(-r * k0);
	
	return couplingConstant * exponentialPart * (k0 + 1/r) / r;
}
static void FCoulomb(Particle *p1, Particle *p2)
{	
	if (!interactions.enableCoulomb)
		return;

	if (!isChargedPair(p1->type, p2->type))
		return;

	double truncLen = interactions.truncationLen;
	double r = nearestImageDistance(p2->pos, p1->pos);
	if (r > truncLen)
		return; /* Too far away */
	
	Vec3 direction = nearestImageUnitVector(p1->pos, p2->pos);
	Vec3 forceVec;
	
	forceVec = scale(direction, calcFCoulomb(r));
	debugVectorSanity(forceVec, "FCoulomb");
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
	FbondSB(&s->Ss[0], &s->Bs[0]);
	Fbond(&s->Ss[0], &s->Ps[0], BOND_S5_P);
	FangleP5SB(&s->Ps[0], &s->Ss[0], &s->Bs[0]);

	/* Rest of the monomers */
	for (int i = 1; i < s->numMonomers; i++) {
		FbondSB(&s->Ss[i], &s->Bs[i]);
		Fbond(&s->Ss[i], &s->Ps[i],   BOND_S5_P);
		Fbond(&s->Ss[i], &s->Ps[i-1], BOND_S3_P);

		Fstack(&s->Bs[i], &s->Bs[i-1], 1);

		FangleP5SB(&s->Ps[ i ], &s->Ss[ i ], &s->Bs[ i ]);
		Fangle    (&s->Ps[ i ], &s->Ss[ i ], &s->Ps[i-1], ANGLE_P_5S3_P);
		FangleP3SB(&s->Ps[i-1], &s->Ss[ i ], &s->Bs[ i ]);
		Fangle    (&s->Ss[i-1], &s->Ps[i-1], &s->Ss[ i ], ANGLE_S5_P_3S);

		FdihedralBS3P5S(&s->Bs[i], &s->Ss[ i ], &s->Ps[i-1], &s->Ss[i-1]);
		FdihedralS3P5SB(&s->Ss[i], &s->Ps[i-1], &s->Ss[i-1], &s->Bs[i-1]);
		Fdihedral(&s->Ps[i], &s->Ss[ i ], &s->Ps[i-1], &s->Ss[i-1],
							dihedralP5S3P5S);

		if (i < 2) continue;
		Fdihedral(&s->Ss[i], &s->Ps[i-1], &s->Ss[i-1], &s->Ps[i-2],
							dihedralS3P5S3P);
		Fstack(&s->Bs[i], &s->Bs[i-2], 2);
	}
}

static void mutuallyExclusivePairForces(Particle *p1, Particle *p2)
{
	/* Nonbonded pair interactions are mutually exclusive. See Knotts.
	 * Note that this screws up energy conservation!! */
	if (interactions.enableBond && isBondedBasePair(p1, p2))
		FbasePair(p1, p2); /* TODO MUST BE MUTUALLY EXCLUSIVE WITH 
				      STACKING, that is: two base pairs 
				      cannot form a bp bond if they are 
				      adjacent on the same strand! */
	else if (interactions.enableCoulomb && isChargedPair(p1->type, p2->type))
		FCoulomb(p1, p2);
	else
		Fexclusion(p1, p2);
}

static void pairForces(Particle *p1, Particle *p2)
{
	/* Apply all forces at once, this should conserve energy (for the 
	 * verlet integrator without thermal bath coupling). */
	FbasePair(p1, p2);
	FCoulomb(p1, p2);
	Fexclusion(p1, p2);
}

void calculateForces(void)
{
	/* Reset forces */
	forEveryParticle(&resetForce);

	/* Strand-based forces */
	for (int s = 0; s < world.numStrands; s++)
		strandForces(&world.strands[s]);

	/* Particle-based forces */
	if (interactions.mutuallyExclusivePairForces)
		forEveryPair(&mutuallyExclusivePairForces);
	else
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
double getKineticTemperature(void)
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
	Vb += VbondSB(&s->Ss[0], &s->Bs[0]);
	Vb += Vbond(&s->Ss[0], &s->Ps[0], BOND_S5_P);

	Va += VangleP5SB(&s->Ps[0], &s->Ss[0], &s->Bs[0]);
	for (int i = 1; i < s->numMonomers; i++) {
		Vb += VbondSB(&s->Ss[i], &s->Bs[i]);
		Vb += Vbond(&s->Ss[i], &s->Ps[i],   BOND_S5_P);
		Vb += Vbond(&s->Ss[i], &s->Ps[i-1], BOND_S3_P);

		Vs += Vstack(&s->Bs[i], &s->Bs[i-1], 1);

		Va += VangleP5SB(&s->Ps[ i ], &s->Ss[ i ], &s->Bs[ i ]);
		Va += Vangle    (&s->Ps[ i ], &s->Ss[ i ], &s->Ps[i-1], ANGLE_P_5S3_P);
		Va += VangleP3SB(&s->Ps[i-1], &s->Ss[ i ], &s->Bs[ i ]);
		Va += Vangle    (&s->Ss[i-1], &s->Ps[i-1], &s->Ss[ i ], ANGLE_S5_P_3S);

		Vd += VdihedralBS3P5S(&s->Bs[i], &s->Ss[ i ], &s->Ps[i-1], &s->Ss[i-1]);
		Vd += VdihedralS3P5SB(&s->Ss[i], &s->Ps[i-1], &s->Ss[i-1], &s->Bs[i-1]);
		Vd += Vdihedral(&s->Ps[i], &s->Ss[ i ], &s->Ps[i-1], &s->Ss[i-1],
							dihedralP5S3P5S);

		if (i < 2) continue;
		Vd += Vdihedral(&s->Ss[i], &s->Ps[i-1], &s->Ss[i-1], &s->Ps[i-2],
							dihedralS3P5S3P);
		Vs += Vstack(&s->Bs[i], &s->Bs[i-2], 2);
	}

	pe->bond     += Vb;
	pe->angle    += Va;
	pe->dihedral += Vd;
	pe->stack    += Vs;
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

static void addMomentum(Particle *p, void *data)
{
	Vec3 *deltaP = (Vec3*) data;
	Vec3 P = scale(p->vel, p->m);
	p->vel = scale(add(P, *deltaP), 1/p->m);
}

void killMomentum(void)
{
	Vec3 P = momentum();
	int n = numParticles();
	Vec3 Pcorr = scale(P, -1/n);
	forEveryParticleD(&addMomentum, (void*)&Pcorr);
}


bool physicsCheck(void)
{
	Vec3 P = momentum();
	double PPP = length(P) / numParticles();
	//printf("%e\n",PPP);
	if (PPP > 1e-20) {
		fprintf(stderr, "\nMOMENTUM CONSERVATION VIOLATED! "
				"Momentum per particle: |P| = %e\n", PPP);
		return false;
	}
	return true;
}



/* ===== INFORMATION FUNCTIONS ===== */

void dumpStats()
{
	PotentialEnergies pe = calcPotentialEnergies();
	double K = kineticEnergy() * ENERGY_FACTOR;
	double T = getKineticTemperature();
	double E = K + pe.bond + pe.angle + pe.dihedral + pe.stack + pe.basePair + pe.Coulomb + pe.exclusion;

	printf("E = %e, K = %e, Vb = %e, Va = %e, Vd = %e, Vs = %e, Vbp = %e, Vpp = %e, Ve = %e, T = %f\n",
			E, K, pe.bond, pe.angle, pe.dihedral, pe.stack, pe.basePair, pe.Coulomb, pe.exclusion, T);
}



/* ===== MISC FUNCTIONS ===== */

void initPhysics(void)
{
	initDihedralCache();
}

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

Vec3 getMonomerCOM(Strand *s, int monomer)
{
	Particle *base      = &s->Bs[monomer];
	Particle *sugar     = &s->Ss[monomer];
	Particle *phosphate = &s->Ps[monomer];

	Vec3 COM = add(scale(base->pos,      base->m),
	           add(scale(sugar->pos,     sugar->m),
	               scale(phosphate->pos, phosphate->m)));

	return scale(COM, 1/(base->m + sugar->m + phosphate->m));
}

double parseTemperature(const char *string)
{
	if (string == NULL || string[0] == '\0')
		return -1;

	/* Copy so we can modify characters */
	char *str = calloc(strlen(string) + 1, sizeof(*str));
	if (str == NULL)
		return -1;
	strcpy(str, string);

	int i = 1;
	while (str[i] != '\0') i++;
	i--;
	char unit = str[i];
	str[i] = '\0';
	double value = atof(str);
	free(str);

	double temperature;
	switch (unit) {
		case 'K':
			temperature = value;
			break;
		case 'C':
			temperature = CELSIUS(value);
			break;
		default:
			return -1;
	}
	return temperature;
}
