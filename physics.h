#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <stdbool.h>
#include "system.h"
#include "vmath.h"
#include "task.h"
#include "world.h"


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
#define COUPLING_BP_AT 	1.928e-20  /* 2.77 kcal/mol (per particle) */
#define COUPLING_BP_GC	2.896e-20  /* 4.16 kcal/mol (per particle) */
#define DISTANCE_r0_AT	2.9002e-10 /* Knotts et al 2007, table III, 2.9002 A */
#define DISTANCE_r0_GC	2.8694e-10 /* Knotts et al 2007, table III, 2.8694 A */

/* Coulomb interaction between phosphates */
#define CHARGE_ELECTRON		1.602e-19  /* 1.602 Coulomb */
#define VACUUM_PERMITTIVITY	8.8541e-12 /* 8.854e-12 Farads/m */
#define COUPLING_EPS_H2O	(78*VACUUM_PERMITTIVITY) /* 78 epsilon_0 */
#define DEBYE_LENGTH 		13.603e-10 /* 13.603 Angstrom for 50mM = [Na+]*/

/* Rope dynamics */
#define ROPE_COUPLING		(10 * EPSILON) /* TODO sane value? */
//#define ROPE_COUPLING		0
#define ROPE_DIST		(1.0e-10) /* characteristic length TODO sane? */
#define ROPE_TRUNCATION		3e-10 /* 3 Angstrom */

/* Thermodynamics */
#define ENERGY_FACTOR	      (1/1.602177e-19) /* Energy in electronvolt */
#define BOLTZMANN_CONSTANT    1.38065e-23
#define FROM_ANGSTROM_SQUARED 1e20 /* Bond constants are given for Angstrom */
#define TO_RADIANS	      (M_PI / 180)

#define DIELECTRIC_CST_H20    80
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

Task makeIntegratorTask(IntegratorConf *conf);

#endif
