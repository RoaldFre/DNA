#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <stdbool.h>
#include "system.h"
#include "math.h"
#include "task.h"
#include "world.h"


/* Units */
#define DEGREE		(M_PI / 180)	/* Convert to radians */
#define ANGSTROM	(1e-10)
#define AU		(1.660539e-27)	/* Atomic unit */
#define MILLISECONDS	1e-3
#define MICROSECONDS	1e-6
#define NANOSECONDS	1e-9
#define PICOSECONDS	1e-12
#define FEMTOSECONDS	1e-15
#define KCAL_PER_MOL	6.947695e-21	/* kcal/mol in Joule(/particle) */
#define EPSILON		(0.26 * KCAL_PER_MOL)
#define ELECTRON_CHARGE	1.602177e-19
#define ELECTRON_VOLT	(ELECTRON_CHARGE)
#define CELSIUS_TO_KELVIN 273.15
#define CELSIUS(T)	(T + CELSIUS_TO_KELVIN)
#define TO_CELSIUS(T)	(T - CELSIUS_TO_KELVIN)


/* STRUCTURE OF THE HELIX. SEE TABLE I IN KNOTTS */
/* Vertical spacing between layers */
#define HELIX_DELTA_Z   (3.38 * ANGSTROM)
/* Twist at each consecutive layer */
#define HELIX_DELTA_PHI (36 * DEGREE)
/* Angles */
#define P_PHI (94.038 * DEGREE)
#define S_PHI (70.197 * DEGREE)
#define A_PHI (41.905 * DEGREE)
#define T_PHI (86.119 * DEGREE)
#define C_PHI (85.027 * DEGREE)
#define G_PHI (40.691 * DEGREE)
#define X_PHI A_PHI
#define Y_PHI T_PHI
/* Radial distance */
#define P_R (8.916 * ANGSTROM)
#define S_R (6.981 * ANGSTROM)
#define A_R (0.773 * ANGSTROM)
#define T_R (2.349 * ANGSTROM)
#define C_R (2.296 * ANGSTROM)
#define G_R (0.828 * ANGSTROM)
#define X_R A_R
#define Y_R T_R
/* Heights */
#define P_Z (2.186 * ANGSTROM)
#define S_Z (1.280 * ANGSTROM)
#define A_Z (0.051 * ANGSTROM)
#define T_Z (0.191 * ANGSTROM)
#define C_Z (0.187 * ANGSTROM)
#define G_Z (0.053 * ANGSTROM)
#define X_Z A_Z
#define Y_Z T_Z
/* Masses (in kg) */
#define P_M (94.97 * AU)
#define S_M (83.11 * AU)
#define A_M (134.1 * AU)
#define T_M (125.1 * AU)
#define C_M (110.1 * AU)
#define G_M (150.1 * AU)
#define X_M A_M
#define Y_M T_M

/* BOND INTERACTION */
#define BOND_K1   (  1e20 * EPSILON) /* in J/ANGSTROM^2 */
#define BOND_K2   (100e20 * EPSILON) /* in J/ANGSTROM^4 */
#define BOND_S5_P (3.899*ANGSTROM)
#define BOND_S3_P (3.559*ANGSTROM)
#define BOND_S_A  (6.430*ANGSTROM)
#define BOND_S_T  (4.880*ANGSTROM)
#define BOND_S_C  (4.921*ANGSTROM)
#define BOND_S_G  (6.392*ANGSTROM)
#define BOND_S_X BOND_S_A
#define BOND_S_Y BOND_S_T

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
#define ANGLE_P_5S_X ANGLE_P_5S_A
#define ANGLE_P_3S_X ANGLE_P_3S_A
#define ANGLE_P_5S_Y ANGLE_P_5S_T
#define ANGLE_P_3S_Y ANGLE_P_3S_T

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
#define DIHEDRAL_X_S3_P_5S DIHEDRAL_A_S3_P_5S
#define DIHEDRAL_S3_P_5S_X DIHEDRAL_S3_P_5S_A
#define DIHEDRAL_Y_S3_P_5S DIHEDRAL_T_S3_P_5S
#define DIHEDRAL_S3_P_5S_Y DIHEDRAL_S3_P_5S_T

/* BASE PAIR INTERACTION */
#define BASE_PAIR_COUPLING_A_T (2.77 * KCAL_PER_MOL) /* Knotts */
#define BASE_PAIR_COUPLING_G_C (4.16 * KCAL_PER_MOL) /* Knotts */
#define BASE_PAIR_COUPLING_X_Y (10 * BASE_PAIR_COUPLING_A_T)
//#define BASE_PAIR_COUPLING_A_T (3.90 * KCAL_PER_MOL) /* Florescu & Joyeux */
//#define BASE_PAIR_COUPLING_G_C (4.37 * KCAL_PER_MOL) /* Florescu & Joyeux */
#define BASE_PAIR_DISTANCE_A_T (2.9002*ANGSTROM)
#define BASE_PAIR_DISTANCE_G_C (2.8694*ANGSTROM)
#define BASE_PAIR_DISTANCE_X_Y BASE_PAIR_DISTANCE_A_T

/* STACKING INTERACTION */
#define STACK_COUPLING	EPSILON
/* Equilibrium distances get derived from the structure of the equilibrium 
 * helix above */

/* EXLUSION INTERACTION */
#define EXCLUSION_COUPLING	(0.1 * EPSILON)
#define EXCLUSION_DISTANCE	(5.50*ANGSTROM) /* 6.86A in Knotts */
#define EXCLUSION_DISTANCE_BASE	(1.00*ANGSTROM)

/* COULOMB INTERACTION (between phosphates) */
#define VACUUM_PERMITTIVITY	8.8541e-12 /* Farad/m */
#define H2O_PERMETTIVITY	(80 * VACUUM_PERMITTIVITY)

/* THERMODYNAMICS */
#define BOLTZMANN_CONSTANT	1.38065e-23
#define AVOGADRO		6.023e23

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

	enum BasePairInteraction {
		/* Every pair of bases in the world participates in base 
		 * pairing. */
		BASE_PAIR_ALL,

		/* Bases at position i and n-i of the same strand 
		 * participate in base pairing, where n is the length of 
		 * the strand. */
		BASE_PAIR_HAIRPIN,

		/* Base pairing for XY pairs is as BASE_PAIR_HAIRPIN, base 
		 * pairing for the other bases is as BASE_PAIR_ALL. */
		BASE_PAIR_XY_HAIRPIN_REST_ALL,

		/* Bases at position i of one strand pair with bases of 
		 * position i in another strand participate in base 
		 * pairing. */
		BASE_PAIR_DOUBLE_STRAND,
	} basePairInteraction;

	/* Scale de default base pairing strengths with this factor. This 
	 * is useful when disabling, for instance, the dihedral 
	 * interaction, as that also reduces the base pairing stability and 
	 * hence modifies the melting temperatures. */
	double basePairFactor;

	/* Only enable base pairing for XY pairs. The interaction between 
	 * other base pairs that would normally form bonds is now reduced 
	 * to only the repulsive part. Effectively, this behaves as if 
	 * truncationLen would be the distance of the potential minimum (or 
	 * the original truncationLen, if this is smaller than the 
	 * potential minimum distance -- but you probably don't want that 
	 * to happen). */
	bool onlyXYbasePairing;

	double saltConcentration; /* Na+ concentration, in mol/m^3 */
	double truncationLen; /* Length at which potentials are truncated */
} InteractionSettings;

/* This needs to be called before doing any physics calculations! */
void registerInteractions(InteractionSettings interactionSettings);

typedef struct {
	/* Short (one word) description of the interaction. */
	const char *name;

	/* Symbol/suffix for this interaction. This will be appended after 
	 * V to designate the potential energy. For example, for a dihedral 
	 * interaction, this could be 'd', to give a potential Vd. */
	const char *symbol;

	/* Data pointer that holds configuration and/or state that gets 
	 * passed along to the functions below. */
	void *data;

	/* Return the potential energy associated with this interaction. */
	double (*potential)(void *data);

	/* Add the forces of this interaction to the relevant particles. */
	void (*addForces)(void *data);
} ExtraInteraction;

void registerExtraInteraction(ExtraInteraction *interaction);


/* This needs to be called whenever some relevant external setting (like 
 * the heat bath temperature) is changed. */
void syncPhysics(void);

void calculateForces(void);

double getKineticTemperature(void);
double getPotentialEnergy(void);

/* Parse a string of the form "<number><C|K>" and return the corresponding 
 * temperature (in Kelvin). Returns -1 in case of error. */
double parseTemperature(const char *str);

Vec3 endToEndVector(Strand *s);
Vec3 endToEndDirection(Strand *s);
double endToEndDistance(Strand *s);
/* Returns the position vector of the Center Of Mass. */
Vec3 getCOM(Particle *ps, int num);
/* Returns the total mass of the given monomer in the given strand. */
double getMonomerMass(Strand *s, int monomer);
/* Returns the position vector of the Center Of Mass of the given monomer 
 * in the given strand. */
Vec3 getMonomerCOM(Strand *s, int monomer);
/* Returns the momentum vector of the Center Of Mass of the given monomer 
 * in the given strand. */
Vec3 getMonomerMomentum(Strand *s, int monomer);
/* Returns the velocity vector of the Center Of Mass of the given monomer 
 * in the given strand. */
Vec3 getMonomerVelocity(Strand *s, int monomer);
/* Returns the force vector acting on the Center Of Mass of the given 
 * monomer in the given strand. */
Vec3 getMonomerForce(Strand *s, int monomer);
/* Returns the position vector of the Center Of Mass of the given strand */ 
Vec3 getStrandCOM(Strand *s);
/* Distributes the force over the monomer. The force is assumed to be 
 * applied to the COM of the monomer. */
void distributeForceOverMonomer(Vec3 F, Strand *s, int monomer);

/* Sets the total momentum of the world to zero by shifting the velocities. */
void killMomentum(void);

double VbasePair(Particle *p1, Particle *p2);

void dumpStats(void);
bool physicsCheck(void);

#endif
