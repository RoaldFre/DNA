#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <stdbool.h>
#include "system.h"
#include "math.h"
#include "task.h"
#include "world.h"


/* Units */
#define DEGREE		(M_PI / 180)	/* Convert to radians */
#define A		(1e-10)		/* Angstrom */
#define AU		(1.660539e-27)	/* Atomic unit */
#define MILLISECONDS	1e-3
#define MICROSECONDS	1e-6
#define NANOSECONDS	1e-9
#define PICOSECONDS	1e-12
#define FEMTOSECONDS	1e-15
#define EPSILON		1.81e-21 	/* 0.26 kcal/mol = 1.81e-21 J/particle */
#define ELECTRON_CHARGE	1.602177e-19
#define ENERGY_FACTOR	(1/ELECTRON_CHARGE) /* Energies in electron volt */
#define CELSIUS_TO_KELVIN 273.15
#define CELSIUS(T)	(T + CELSIUS_TO_KELVIN)


/* STRUCTURE OF THE HELIX. SEE TABLE I IN KNOTTS */
/* Vertical spacing between layers */
#define HELIX_DELTA_Z   (3.38 * A)
/* Twist at each consecutive layer */
#define HELIX_DELTA_PHI (36 * DEGREE)
/* Angles */
#define P_PHI (94.038 * DEGREE)
#define S_PHI (70.197 * DEGREE)
#define A_PHI (41.905 * DEGREE)
#define T_PHI (86.119 * DEGREE)
#define C_PHI (85.027 * DEGREE)
#define G_PHI (40.691 * DEGREE)
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

/* BOND INTERACTION */
#define BOND_K1   ( 10e20 * EPSILON) /* in J/A^2 */
#define BOND_K2   (100e20 * EPSILON) /* in J/A^4 */
#define BOND_S5_P (3.899*A)
#define BOND_S3_P (3.559*A)
#define BOND_S_A  (6.430*A)
#define BOND_S_T  (4.880*A)
#define BOND_S_C  (4.921*A)
#define BOND_S_G  (6.392*A)

/* ANGLE INTERACTION */
#define ANGLE_COUPLING	(400 * EPSILON) /* per radian^2 */
#define ANGLE_S5_P_3S	( 94.49 * DEGREE)
#define ANGLE_P_5S3_P	(120.15 * DEGREE)
#define ANGLE_P_5S_A	(113.13 * DEGREE)
#define ANGLE_P_3S_A	(108.38 * DEGREE)
#define ANGLE_P_5S_T	(102.79 * DEGREE)
#define ANGLE_P_3S_T	(112.72 * DEGREE)
#define ANGLE_P_5S_C	(103.49 * DEGREE)
#define ANGLE_P_3S_C	(112.39 * DEGREE)
#define ANGLE_P_5S_G	(113.52 * DEGREE)
#define ANGLE_P_3S_G	(108.12 * DEGREE)

/* DIHEDRAL INTERACTION */
#define DIHEDRAL_COUPLING	(4 * EPSILON)
#define DIHEDRAL_P_5S3_P_5S	(-154.80 * DEGREE)
#define DIHEDRAL_S3_P_5S3_P	(-179.17 * DEGREE)
#define DIHEDRAL_A_S3_P_5S	( -22.60 * DEGREE)
#define DIHEDRAL_S3_P_5S_A	(  50.69 * DEGREE)
#define DIHEDRAL_T_S3_P_5S	( -33.42 * DEGREE)
#define DIHEDRAL_S3_P_5S_T	(  54.69 * DEGREE)
#define DIHEDRAL_C_S3_P_5S	( -32.72 * DEGREE)
#define DIHEDRAL_S3_P_5S_C	(  54.50 * DEGREE)
#define DIHEDRAL_G_S3_P_5S	( -22.30 * DEGREE)
#define DIHEDRAL_S3_P_5S_G	(  50.66 * DEGREE)

/* BASE PAIR INTERACTION */
#define BASE_PAIR_COUPLING_A_T 	1.928e-20
#define BASE_PAIR_COUPLING_G_C	2.896e-20
#define BASE_PAIR_DISTANCE_A_T	(2.9002*A)
#define BASE_PAIR_DISTANCE_G_C	(2.8694*A)

/* STACKING INTERACTION */
#define STACK_COUPLING	EPSILON
/* Equilibrium distances get derived from the structure of the equilibrium 
 * helix above */

/* EXLUSION INTERACTION */
#define EXCLUSION_COUPLING	(0.1 * EPSILON)
#define EXCLUSION_DISTANCE	(5.50*A) /* 6.86A in Knotts */
#define EXCLUSION_DISTANCE_BASE	(1.00*A)

/* COULOMB INTERACTION (between phosphates) */
#define VACUUM_PERMITTIVITY	8.8541e-12 /* Farad/m */
#define H2O_PERMETTIVITY	(80 * VACUUM_PERMITTIVITY)

/* THERMODYNAMICS */
#define BOLTZMANN_CONSTANT	1.38065e-23
#define AVOGADRO		6.023e23


typedef enum {
	LANGEVIN,
	VERLET
} Integrator;

typedef struct {
	bool enableBond;
	bool enableAngle;
	bool enableDihedral;
	bool enableStack;
	bool enableExclusion;
	bool enableBasePair;
	bool enableCoulomb;

	/* True for Knotts' model */
	bool mutuallyExclusivePairForces;

	enum {
		/* Every pair of bases in the world participates in base 
		 * pairing. */
		BASE_PAIR_ALL,

		/* Bases at position i and n-i of the same strand 
		 * participate in base pairing, where n is the length of 
		 * the strand. */
		BASE_PAIR_HAIRPIN,

		/* Bases at position i of one strand pair with bases of 
		 * position i in another strand participate in base 
		 * pairing. */
		BASE_PAIR_DOUBLE_STRAND,
	} basePairInteraction;
} InteractionSettings;

typedef struct {
	Integrator integrator;
	InteractionSettings interactionSettings;
	int numBoxes; /* For space partition grid */
	double reboxInterval; /* time interval after which to rebox the particles */
} IntegratorConf;

void dumpStats(void);
bool physicsCheck(void);

double temperature(void);
/* Returns the position vector of the Center Of Mass. */
Vec3 getCOM(Particle *ps, int num);
/* Sets the total momentum of the world to zero by shifting the velocities. */
void killMomentum(void);

double VbasePair(Particle *p1, Particle *p2);

double nearestLineDistance(Vec3 *pos1, Vec3 *pos2, Vec3 *dist1, Vec3 *dist2);
double getExclusionCutOff(ParticleType t1, ParticleType t2);

Task makeIntegratorTask(IntegratorConf *conf);

#endif
