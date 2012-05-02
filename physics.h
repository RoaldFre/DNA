#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <stdbool.h>
#include "system.h"
#include "vmath.h"
#include "task.h"
#include "world.h"


/* Units */
#define RAD	(M_PI / 180)	/* Radian */
#define A	(1e-10)		/* Angstrom */
#define AU	(1.660539e-27)	/* Atomic unit */

/* Energy unit */
#define EPSILON 1.81e-21 /* 0.26kcal/mol == 1.81 * 10^-21 J (per particle) */


/* Structure of the helix. See table I in Knotts */

#define HELIX_DELTA_Z   (3.38 * A) 	/* vertical spacing between layers */
#define HELIX_DELTA_PHI (36 * RAD)	/* twist at each consecutive layer */

/* Angles */
#define P_PHI (94.038 * RAD)
#define S_PHI (70.197 * RAD)
#define A_PHI (41.905 * RAD)
#define T_PHI (86.119 * RAD)
#define C_PHI (85.027 * RAD)
#define G_PHI (40.691 * RAD)
/* Radial distance */
#define P_R (8.916 * A)
#define S_R (6.981 * A)
#define A_R (0.773 * A)
#define T_R (2.349 * A)
#define C_R (2.296 * A)
#define G_R (0.828 * A)
/* Heights */
#define P_Z (2.186 * A)
#define S_Z (1.280 * A)
#define A_Z (0.051 * A)
#define T_Z (0.191 * A)
#define C_Z (0.187 * A)
#define G_Z (0.053 * A)
/* Masses (in kg) */
#define P_M (94.97 * AU)
#define S_M (83.11 * AU)
#define A_M (134.1 * AU)
#define T_M (125.1 * AU)
#define C_M (110.1 * AU)
#define G_M (150.1 * AU)

/* Equilibrium distance of bonds (in m) */
#define D_S5P   3.899*A
#define D_S3P   3.559*A
#define D_SA    6.430*A



/* TODO rename coupling constants for consistency! */

/* Bond stretch */
#define BOND_K1      (10e20 * EPSILON) /* in J*A^-2 */
#define BOND_K2      (100e20 * EPSILON) /* in J*A^-2 */
/* Bond bend */
#define BOND_Ktheta  (400 * EPSILON) /* per radian^2 */
/* Bond dihedral */
#define BOND_Kphi    (4 * EPSILON)
/* Bond stack */
#define BOND_STACK   EPSILON

/* Bond angle */
#define ANGLE_S5_P_3S	( 94.49 * RAD)
#define ANGLE_P_5S3_P	(120.15 * RAD)
#define ANGLE_P_5S_A	(113.13 * RAD)
#define ANGLE_P_3S_A	(108.38 * RAD)

/* Dihedral angle */
#define DIHEDRAL_P_5S3_P_5S	(-154.80 * RAD)
#define DIHEDRAL_S3_P_5S3_P	(-179.17 * RAD)
#define DIHEDRAL_A_S3_P_5S	( -22.60 * RAD)
#define DIHEDRAL_S3_P_5S_A	(  50.69 * RAD)

/* Base-Pair couplings */
#define COUPLING_BP_AT 	1.928e-20  /* 2.77 kcal/mol (per particle) */
#define COUPLING_BP_GC	2.896e-20  /* 4.16 kcal/mol (per particle) */
#define DISTANCE_r0_AT	(2.9002*A) /* Knotts et al 2007, table III, 2.9002 A */
#define DISTANCE_r0_GC	(2.8694*A) /* Knotts et al 2007, table III, 2.8694 A */

/* Exlusion Force */
#define EXCLUSION_COUPLING	(0.1 * EPSILON)
#define D_CUT			(6.86*A)
#define D_CUT_BASE		(1.00*A)

/* Coulomb interaction between phosphates */
#define CHARGE_ELECTRON		1.602e-19  /* Coulomb */
#define VACUUM_PERMITTIVITY	8.8541e-12 /* Farad/m */
#define H2O_PERMETTIVITY	(80 * VACUUM_PERMITTIVITY)
#define DEBYE_LENGTH 		13.603e-10 /* 13.603 Angstrom for 50mM = [Na+]*/

/* Thermodynamics */
#define ENERGY_FACTOR	      (1/1.602177e-19) /* Energy in electronvolt */
#define BOLTZMANN_CONSTANT    1.38065e-23

#define AVOGADRO              6.023e23 /* particles per mol */




typedef enum
{
	LANGEVIN,
	VERLET
} Integrator;

typedef struct
{
	Integrator integrator;
	int numBoxes; /* For space partition grid */
} IntegratorConf;

void dumpStats(void);
bool physicsCheck(void);

double temperature(void);
/* Returns the position vector of the Center Of Mass. */
Vec3 getCOM(Particle *ps, int num);

double VbasePair(Particle *p1, Particle *p2);

double nearestLineDistance(Vec3 *pos1, Vec3 *pos2, Vec3 *dist1, Vec3 *dist2);
double getExclusionCutOff(ParticleType t1, ParticleType t2);

Task makeIntegratorTask(IntegratorConf *conf);

#endif
