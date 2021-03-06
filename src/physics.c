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


/* ===== EXTRA INTERACTIONS ===== */
typedef struct extraIntNode {
	ExtraInteraction i;
	struct extraIntNode *next;
} ExtraIntNode;

static ExtraIntNode *extraIntList = NULL;

void registerExtraInteraction(ExtraInteraction *interaction)
{
	ExtraIntNode *node = malloc(sizeof(*node));
	node->i = *interaction; /* struct copy */
	node->next = extraIntList;
	extraIntList = node;
}



/* ===== CACHED STUFF FOR PERFORMANCE ===== */
static double invDebLength; /* Inverse Debye length for Coulomb screening */
static double calcInvDebyeLength(void);
void syncPhysics(void)
{
	invDebLength = calcInvDebyeLength();
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
static double sugarBaseBondDistanceLUT[NUM_BASE_TYPES];
static void initSugarBaseBondDistanceLUT(void) {
	for (int b = 0; b < NUM_BASE_TYPES; b++) {
		assert(isBase(b));
		switch (b) {
		case BASE_A:
			sugarBaseBondDistanceLUT[b] = BOND_S_A;
			break;
		case BASE_T:
			sugarBaseBondDistanceLUT[b] = BOND_S_T;
			break;
		case BASE_C:
			sugarBaseBondDistanceLUT[b] = BOND_S_C;
			break;
		case BASE_G:
			sugarBaseBondDistanceLUT[b] = BOND_S_G;
			break;
		case BASE_X:
			sugarBaseBondDistanceLUT[b] = BOND_S_X;
			break;
		case BASE_Y:
			sugarBaseBondDistanceLUT[b] = BOND_S_Y;
			break;
		default:
			die("Unknown base type %d in "
					"initSugarBaseBondDistanceLUT!\n", b);
		}
	}
}
static double getSugarBaseBondDistance(ParticleType base)
{
	assert(isBase(base));
	return sugarBaseBondDistanceLUT[base];
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
static AngleBaseInfo angleBaseInfoLUT[NUM_BASE_TYPES];
static void initAngleBaseInfoLUT(void) {
	for (int b = 0; b < NUM_BASE_TYPES; b++) {
		assert(isBase(b));
		switch (b) {
		case BASE_A:
			angleBaseInfoLUT[b].P5SB = ANGLE_P_5S_A;
			angleBaseInfoLUT[b].P3SB = ANGLE_P_3S_A;
			break;
		case BASE_T:
			angleBaseInfoLUT[b].P5SB = ANGLE_P_5S_T;
			angleBaseInfoLUT[b].P3SB = ANGLE_P_3S_T;
			break;
		case BASE_C:
			angleBaseInfoLUT[b].P5SB = ANGLE_P_5S_C;
			angleBaseInfoLUT[b].P3SB = ANGLE_P_3S_C;
			break;
		case BASE_G:
			angleBaseInfoLUT[b].P5SB = ANGLE_P_5S_G;
			angleBaseInfoLUT[b].P3SB = ANGLE_P_3S_G;
			break;
		case BASE_X:
			angleBaseInfoLUT[b].P5SB = ANGLE_P_5S_X;
			angleBaseInfoLUT[b].P3SB = ANGLE_P_3S_X;
			break;
		case BASE_Y:
			angleBaseInfoLUT[b].P5SB = ANGLE_P_5S_Y;
			angleBaseInfoLUT[b].P3SB = ANGLE_P_3S_Y;
			break;
		default:
			die("Unknown base type %d in initAngleBaseInfoLUT!\n", b);
		}
	}
}
static AngleBaseInfo getAngleBaseInfo(ParticleType base)
{
	assert(isBase(base));
	return angleBaseInfoLUT[base];
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

static DihedralBaseInfo dihedralsBasesLUT[NUM_BASE_TYPES];
static DihedralCache dihedralP5S3P5S;
static DihedralCache dihedralS3P5S3P;

static DihedralCache makeDihedralCache(double dihedralAngle)
{
	DihedralCache ret;
	ret.sinDihedral = sin(dihedralAngle);
	ret.cosDihedral = cos(dihedralAngle);
	return ret;
}
static void initDihedralLUT(void)
{
	/* Backbone */
	dihedralP5S3P5S = makeDihedralCache(DIHEDRAL_P_5S3_P_5S);
	dihedralS3P5S3P = makeDihedralCache(DIHEDRAL_S3_P_5S3_P);

	/* Everything involving bases */
	for (int b = 0; b < NUM_BASE_TYPES; b++) {
		assert(isBase(b));
		switch(b) {
		case BASE_A:
			dihedralsBasesLUT[b].BS3P5S = makeDihedralCache(DIHEDRAL_A_S3_P_5S);
			dihedralsBasesLUT[b].S3P5SB = makeDihedralCache(DIHEDRAL_S3_P_5S_A);
			break;
		case BASE_T:
			dihedralsBasesLUT[b].BS3P5S = makeDihedralCache(DIHEDRAL_T_S3_P_5S);
			dihedralsBasesLUT[b].S3P5SB = makeDihedralCache(DIHEDRAL_S3_P_5S_T);
			break;
		case BASE_C:
			dihedralsBasesLUT[b].BS3P5S = makeDihedralCache(DIHEDRAL_C_S3_P_5S);
			dihedralsBasesLUT[b].S3P5SB = makeDihedralCache(DIHEDRAL_S3_P_5S_C);
			break;
		case BASE_G:
			dihedralsBasesLUT[b].BS3P5S = makeDihedralCache(DIHEDRAL_G_S3_P_5S);
			dihedralsBasesLUT[b].S3P5SB = makeDihedralCache(DIHEDRAL_S3_P_5S_G);
			break;
		case BASE_X:
			dihedralsBasesLUT[b].BS3P5S = makeDihedralCache(DIHEDRAL_X_S3_P_5S);
			dihedralsBasesLUT[b].S3P5SB = makeDihedralCache(DIHEDRAL_S3_P_5S_X);
			break;
		case BASE_Y:
			dihedralsBasesLUT[b].BS3P5S = makeDihedralCache(DIHEDRAL_Y_S3_P_5S);
			dihedralsBasesLUT[b].S3P5SB = makeDihedralCache(DIHEDRAL_S3_P_5S_Y);
			break;
		default:
			die("Unknown base type %d in "
					"initDihedralLUT!\n", b);
		}
	}
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
	assert(b->type < NUM_BASE_TYPES);
	return Vdihedral(b, s1, p, s2, dihedralsBasesLUT[b->type].BS3P5S);
}
static double VdihedralS3P5SB(Particle *s1, Particle *p,
				Particle *s2, Particle *b)
{
	assert(b->type < NUM_BASE_TYPES);
	return Vdihedral(s1, p, s2, b, dihedralsBasesLUT[b->type].S3P5SB);
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
		die("Internal error!\n");
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
}

static void FdihedralBS3P5S(Particle *b, Particle *s1,
				Particle *p, Particle *s2)
{
	Fdihedral(b, s1, p, s2, dihedralsBasesLUT[b->type].BS3P5S);
}
static void FdihedralS3P5SB(Particle *s1, Particle *p,
				Particle *s2, Particle *b)
{
	Fdihedral(s1, p, s2, b, dihedralsBasesLUT[b->type].S3P5SB);
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
	case BASE_X:	info.r = X_R;  info.phi = X_PHI;  info.z = X_Z;  break;
	case BASE_Y:	info.r = Y_R;  info.phi = Y_PHI;  info.z = Y_Z;  break;
	default: die("Error in getHelixInfo()!\n");
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
static double neighbourStackDistance2LUT[NUM_BASE_TYPES][NUM_BASE_TYPES][2];
static double neighbourStackCorrectionLUT[NUM_BASE_TYPES][NUM_BASE_TYPES][2];
static void initNeighbourStackLUTs(void)
{
	assert(truncationLenSq != 0); // This must be set before initializing

	for (int b1 = 0; b1 < NUM_BASE_TYPES; b1++)
	for (int b2 = 0; b2 < NUM_BASE_TYPES; b2++) {
		HelixInfo hi1 = getHelixInfo(b1);
		HelixInfo hi2 = getHelixInfo(b2);
		double rEqSq1 = helixDistance2(hi2, hi1,
					HELIX_DELTA_Z, HELIX_DELTA_PHI);
		double rEqSq2 = helixDistance2(hi2, hi1,
					2 * HELIX_DELTA_Z, 2 * HELIX_DELTA_PHI);
		double Vcorr1 = calcVLJ(STACK_COUPLING, rEqSq1, truncationLenSq);
		double Vcorr2 = calcVLJ(STACK_COUPLING, rEqSq2, truncationLenSq);
		neighbourStackDistance2LUT[b1][b2][0] = rEqSq1;
		neighbourStackDistance2LUT[b1][b2][1] = rEqSq2;
		neighbourStackCorrectionLUT[b1][b2][0] = Vcorr1;
		neighbourStackCorrectionLUT[b1][b2][1] = Vcorr2;
	}
}

/* Distance between (next) neighbouring bases in the beta helix, ie base i 
 * and i+1 or i and i+2.
 * Base 1 is the 'bottom' base in the helix (lower index in the strand), 
 * base 2 is the 'top' base of the helix (higher index in the strand).
 *
 * TODO WARNING, complementary strand should be of the correct 3'-5' 
 * "direction" as the "normal" strand, or this breaks down.
 *
 * monomerDistance:
 *   0 for immediate neighbours i and i+1
 *   1 for   next    neighbours i and i+2 */
static double neighbourStackDistance2(ParticleType b1, ParticleType b2,
							int monomerDistance)
{
	assert(isBase(b1) && isBase(b2));
	assert(0 <= monomerDistance && monomerDistance <= 1);
	return neighbourStackDistance2LUT[b1][b2][monomerDistance];
}
static double neighbourStackCorrection(ParticleType b1, ParticleType b2,
							int monomerDistance)
{
	assert(isBase(b1) && isBase(b2));
	assert(0 <= monomerDistance && monomerDistance <= 1);
	return neighbourStackCorrectionLUT[b1][b2][monomerDistance];
}
/* monomerDistance:
 *   0 for immediate neighbours (i and i+1)
 *   1 for   next    neighbours (i and i+2) */
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
	double Vcorr = neighbourStackCorrection(p1->type, p2->type,
						monomerDistance);
	double V = calcVLJ(STACK_COUPLING, rEqSq, rSq) - Vcorr;
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
	bool onlyRepulsive; /* Use only repulsive part of interaction? */
} BasePairInfo;
/* Returns the BasePairInfo for the given types. If the given types do not 
 * constitute a valid base pair, then coupling and distance are set to -1. */
static BasePairInfo basePairInfoLUT[NUM_BASE_TYPES][NUM_BASE_TYPES];
static void initBasePairInfoLUT(void) {
	double factor = interactions.basePairFactor;

	BasePairInfo bpiGC, bpiAT, bpiXY, bpiNULL;
	bpiNULL.coupling = -1;
	bpiNULL.distance2 = -1;

	bpiGC.coupling = BASE_PAIR_COUPLING_G_C * factor;
	bpiGC.distance2 = SQUARE(BASE_PAIR_DISTANCE_G_C);
	bpiGC.onlyRepulsive = (interactions.onlyXYbasePairing ||
	                       interactions.onlyRepulsiveBasePairing);

	bpiAT.coupling = BASE_PAIR_COUPLING_A_T * factor;
	bpiAT.distance2 = SQUARE(BASE_PAIR_DISTANCE_A_T);
	bpiAT.onlyRepulsive = (interactions.onlyXYbasePairing ||
	                       interactions.onlyRepulsiveBasePairing);


	bpiXY.coupling = BASE_PAIR_COUPLING_X_Y * factor;
	bpiXY.distance2 = SQUARE(BASE_PAIR_DISTANCE_X_Y);
	bpiXY.onlyRepulsive = interactions.onlyRepulsiveBasePairing;

	/* Initialize to null */
	for (int b1 = 0; b1 < NUM_BASE_TYPES; b1++)
		for (int b2 = 0; b2 < NUM_BASE_TYPES; b2++)
			basePairInfoLUT[b1][b2] = bpiNULL;

	basePairInfoLUT[BASE_A][BASE_T] = bpiAT;
	basePairInfoLUT[BASE_T][BASE_A] = bpiAT;

	basePairInfoLUT[BASE_G][BASE_C] = bpiGC;
	basePairInfoLUT[BASE_C][BASE_G] = bpiGC;

	basePairInfoLUT[BASE_X][BASE_Y] = bpiXY;
	basePairInfoLUT[BASE_Y][BASE_X] = bpiXY;
}
static BasePairInfo getBasePairInfo(Particle *p1, Particle *p2)
{
	ParticleType t1 = p1->type;
	ParticleType t2 = p2->type;

	BasePairInfo bpiNULL;
	bpiNULL.coupling = -1;
	bpiNULL.distance2 = -1;

	if (!isBase(t1) || !isBase(t2))
		return bpiNULL;

	switch (interactions.basePairInteraction) {
	case BASE_PAIR_ALL:
		break;
	case BASE_PAIR_XY_HAIRPIN_REST_ALL:
		/* If we don't have an XY pair... */
		if ( !(    (t1 == BASE_X && t2 == BASE_Y)
		        || (t1 == BASE_Y && t2 == BASE_X)))
			break; /* ...then behave like BASE_PAIR_ALL */
		/* We have an XY pair -> behave like BASE_PAIR_HAIRPIN */
		/* <Intentional case fall through!> */
	case BASE_PAIR_HAIRPIN:
		if (p1->strand != p2->strand)
			return bpiNULL;

		int n  = p1->strand->numMonomers;
		int i1 = p1->strandIndex;
		int i2 = p2->strandIndex;

		if (i1 != n - 1 - i2)
			return bpiNULL;
		break;
	case BASE_PAIR_DOUBLE_STRAND:
		if (p1->strand == p2->strand)
			return bpiNULL;
		if (p1->strandIndex != p2->strandIndex)
			return bpiNULL;
		break;
	default:
		die("getBasePairInfo: Unknown basePairInteraction!\n");
	}

	assert(isBase(t1) && isBase(t2));
	return basePairInfoLUT[t1][t2];
}
/* Returns true for pairs of particles that can form a base pair binding. */
static bool canFormBasePair(Particle *p1, Particle *p2)
{
	BasePairInfo bpi = getBasePairInfo(p1, p2);
	return bpi.coupling > 0;
}
static double calcVbasePair(BasePairInfo bpi, double rsquared)
{
	/* Note: in VbasePair, we assume that this potential is zero for 
	 * rsquared == bpi.distance2 ! */
	double rfrac2 = bpi.distance2 / rsquared;
	double rfrac4 = rfrac2 * rfrac2;
	double rfrac8 = rfrac4 * rfrac4;
	double rfrac10 = rfrac8 * rfrac2;
	double rfrac12 = rfrac8 * rfrac4;
	
	return bpi.coupling * (5*rfrac12 - 6*rfrac10);
}
double VbasePair(Particle *p1, Particle *p2)
{
	if (!interactions.enableBasePair)
		return 0;

	BasePairInfo bpi = getBasePairInfo(p1, p2);
	if (bpi.coupling < 0)
		return 0; /* Wrong pair */

	double truncSq = truncationLenSq;
	if (bpi.onlyRepulsive)
		truncSq = MIN(truncSq, bpi.distance2);
	
	double rsq = nearestImageDistance2(p1->pos, p2->pos);
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
	if (!interactions.enableBasePair)
		return;

	BasePairInfo bpi = getBasePairInfo(p1, p2);
	if (bpi.coupling < 0)
		return; /* Wrong pair */
	
	Vec3 rVec = nearestImageVector(p1->pos, p2->pos);
	double r = length(rVec);
	if (   (r >= interactions.truncationLen)
	    || (bpi.onlyRepulsive && SQUARE(r) >= bpi.distance2))
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
 * direct bond.
 * Lookup table: areConnectedLUT[indexDiff][type1][type2],
 * where * type1 and type2 are the types of particles on the same strand, 
 *         separated at most 1 monomer from each other,
 *       * indexDiff is the difference in monomer-index of those particles 
 *         in the strand (can only be 0 or 1: particles should be ordered so 
 *         that type1 has a lower index than type2) */
static bool areConnectedLUT[2][NUM_PARTICLE_TYPES][NUM_PARTICLE_TYPES];

static void initAreConnectedLUT(void) {
	/* Initialize everything to false */
	for (int i = 0; i < NUM_PARTICLE_TYPES; i++) {
		for (int j = 0; j < NUM_PARTICLE_TYPES; j++) {
			areConnectedLUT[0][i][j] = false;
			areConnectedLUT[1][i][j] = false;
		}
	}

	/* Everything connected that doesn't involve bases */
	areConnectedLUT[0][SUGAR][PHOSPHATE] = true;
	areConnectedLUT[0][PHOSPHATE][SUGAR] = true;
	areConnectedLUT[1][PHOSPHATE][SUGAR] = true;

	/* Everything involving bases */
	for (int b = 0; b < NUM_BASE_TYPES; b++) {
		assert(isBase(b));

		/* Base is only connected to sugar on same monomer */
		areConnectedLUT[0][b][SUGAR] = true;
		areConnectedLUT[0][SUGAR][b] = true;
	}
}

static bool feelExclusion(Particle *p1, Particle *p2)
{
	/* If different strands: certainly exclusion */
	if (p1->strand != p2->strand)
		return true;

	/* Order so that p1->strandIndex <= p2->strandIndex */
	if (p1->strandIndex > p2->strandIndex) {
		Particle *tmp = p1;
		p1 = p2;
		p2 = tmp;
	}

	int indexDiff = p2->strandIndex - p1->strandIndex;

	/* If separated by more than 1 monomer: certainly exclusion */
	if (indexDiff > 1)
		return true;

	return !areConnectedLUT[indexDiff][p1->type][p2->type];
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
	return sqrt(8 * M_PI * lambdaB * AVOGADRO * saltCon);
}
static double calcVCoulomb(double r)
{
	double couplingConstant = SQUARE(ELECTRON_CHARGE)
				/ (4 * M_PI * H2O_PERMETTIVITY);
	double exponentialPart = exp(-r * invDebLength);
	
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
	double exponentialPart = exp(-r * invDebLength);
	
	return couplingConstant * exponentialPart * (invDebLength + 1/r) / r;
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

		Fstack(&s->Bs[i], &s->Bs[i-1], 0);

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
		Fstack(&s->Bs[i], &s->Bs[i-2], 1);
	}
}

static void mutuallyExclusivePairForces(Particle *p1, Particle *p2)
{
	/* Nonbonded pair interactions are mutually exclusive. See Knotts.
	 * Note that this screws up energy conservation!! */
	if (interactions.enableBasePair && canFormBasePair(p1, p2))
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

	/* Extra interactions */
	ExtraIntNode *node = extraIntList;
	while (node != NULL) {
		node->i.addForces(node->i.data);
		node = node->next;
	}
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

		Vs += Vstack(&s->Bs[i], &s->Bs[i-1], 0);

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
		Vs += Vstack(&s->Bs[i], &s->Bs[i-2], 1);
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
	
	return pe;
}

double getPotentialEnergy(void) {
	PotentialEnergies pe = calcPotentialEnergies();
	double V = pe.bond + pe.angle + pe.dihedral + pe.stack
			+ pe.basePair + pe.Coulomb + pe.exclusion;


	/* Extra potentials */
	ExtraIntNode *node = extraIntList;
	while (node != NULL) {
		V += node->i.potential(node->i.data);
		node = node->next;
	}

	return V;
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
	p->vel = add(p->vel, scale(*deltaP, 1/p->m));
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
	//printf("%le\n",PPP);
	if (PPP > 1e-20) {
		fprintf(stderr, "\nMOMENTUM CONSERVATION VIOLATED! "
				"Momentum per particle: |P| = %le\n", PPP);
		return false;
	}
	return true;
}



/* ===== INFORMATION FUNCTIONS ===== */

void dumpStats()
{
	PotentialEnergies pe = calcPotentialEnergies();
	/* Output in electron volt */
	pe.bond      /= ELECTRON_VOLT;
	pe.angle     /= ELECTRON_VOLT;
	pe.dihedral  /= ELECTRON_VOLT;
	pe.stack     /= ELECTRON_VOLT;
	pe.basePair  /= ELECTRON_VOLT;
	pe.Coulomb   /= ELECTRON_VOLT;
	pe.exclusion /= ELECTRON_VOLT;
	double K = kineticEnergy() / ELECTRON_VOLT;
	double T = getKineticTemperature();
	double E = K + pe.bond + pe.angle + pe.dihedral + pe.stack + pe.basePair + pe.Coulomb + pe.exclusion;

	printf("Vb = %le, Va = %le, Vd = %le, Vs = %le, Vbp = %le, Vpp = %le, Ve = %le",
			pe.bond, pe.angle, pe.dihedral, pe.stack, pe.basePair, pe.Coulomb, pe.exclusion);

	/* Extra interactions */
	ExtraIntNode *node = extraIntList;
	while (node != NULL) {
		double V = node->i.potential(node->i.data) / ELECTRON_VOLT;
		printf(", V%s = %le\n", node->i.symbol, V);
		E += V;
		node = node->next;
	}

	printf(", E = %le, K = %le, T = %lf\n", E, K, T);
}



/* ===== MISC FUNCTIONS ===== */
void registerInteractions(InteractionSettings interactionSettings)
{
	interactions = interactionSettings;
	truncationLenSq = SQUARE(interactionSettings.truncationLen);

	initDihedralLUT();
	initAreConnectedLUT();
	initSugarBaseBondDistanceLUT();
	initAngleBaseInfoLUT();
	initBasePairInfoLUT();
	initNeighbourStackLUTs();

	syncPhysics();
}
InteractionSettings getInteractionSettings(void)
{
	return interactions;
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

Vec3 getStrandCOM(Strand *s)
{
	return getCOM(s->all, 3*s->numMonomers);
}

double getMonomerMass(Strand *s, int monomer)
{
	Particle *base      = &s->Bs[monomer];
	Particle *sugar     = &s->Ss[monomer];
	Particle *phosphate = &s->Ps[monomer];
	return base->m + sugar->m + phosphate->m;
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

Vec3 getMonomerMomentum(Strand *s, int monomer)
{
	Particle *base      = &s->Bs[monomer];
	Particle *sugar     = &s->Ss[monomer];
	Particle *phosphate = &s->Ps[monomer];

	Vec3 P = add(scale(base->vel,      base->m),
	         add(scale(sugar->vel,     sugar->m),
	             scale(phosphate->vel, phosphate->m)));

	return P;
}

Vec3 getMonomerVelocity(Strand *s, int monomer)
{
	return scale(getMonomerMomentum(s, monomer),
			1/getMonomerMass(s, monomer));
}

Vec3 getMonomerForce(Strand *s, int monomer)
{
	Particle *base      = &s->Bs[monomer];
	Particle *sugar     = &s->Ss[monomer];
	Particle *phosphate = &s->Ps[monomer];

	return add(add(base->F, sugar->F), phosphate->F);
}

void distributeForceOverMonomer(Vec3 F, Strand *strand, int monomer)
{
	Particle *b = &strand->Bs[monomer];
	Particle *s = &strand->Ss[monomer];
	Particle *p = &strand->Ps[monomer];

	double M = b->m + s->m + p->m;

	b->F = add(b->F, scale(F, b->m / M));
	s->F = add(s->F, scale(F, s->m / M));
	p->F = add(p->F, scale(F, p->m / M));
}

Vec3 endToEndVector(Strand *s)
{
	return nearestImageVector(getMonomerCOM(s, 0),
	                          getMonomerCOM(s, s->numMonomers - 1));
}
Vec3 endToEndDirection(Strand *s)
{
	return normalize(endToEndVector(s));
}
double endToEndDistance(Strand *s)
{
	return length(endToEndVector(s));
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
