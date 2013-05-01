#ifndef _SPGRID_H_
#define _SPGRID_H_

/* SPgrid: Space Partition grid */

#include "world.h"
#include "physics.h" //TODO separate in particle.h (= all we need)?


/* Allocates a grid with 'numBoxes' boxes in each dimension.
 * Precondition: grid can't already be allocated (unless it was freed 
 * afterwards).
 * Returns true on succes, false on failure. */
bool allocGrid(int numBoxes);

/* Allocates a grid (cfr allocGrid()) and adds every particle of the world. */
bool initGrid(int numBoxes);

/* All particles are removed from the grid and the memory gets freed. */
void freeGrid(void);

/* Adds the given particle to the grid. In the case that the particle is 
 * outside of the grid, periodic boundary conditions are used to force its 
 * position to be within the grid.
 * Precondition: The particle can't already be added to the grid (ie, 
 * p->prev == p->next == NULL).
 */
void addToGrid(Particle *p);

/* Put particles back in their correct boxes in case they escaped. This 
 * also forces periodic boundary conditions on the particle positions in 
 * case the particles escaped from the grid.
 * NOTE: This uses the 'forEveryParticle()' function of world.c. Hence, 
 * this only is useful if ALL particles of the world are already added to 
 * the grid!
 *
 * TODO: Possible safe alternative: Only loop over occupied boxes (i.e. 
 * only particles that are truly added to the grid) -- but that's probably 
 * a tad slower. (There is no problem with 'occupiedboxes' getting 
 * reshuffled while we are reboxing, because particles can only *leave* 
 * boxes that we haven't checked yet, and new boxes are pushed in front of 
 * the linked list, so we won't visit them again when walking the 
 * occupiedboxes list.)*/
void reboxParticles(void);
/* Same as above but for a single particle. This particle must already have 
 * been added to the grid! */
void reboxParticle(Particle *p);

/* Execute a given function for every discinct pair of particles that are 
 * within the same box, or in adjacent boxes (taking into account periodic 
 * boundary conditions).
 * Arguments:
 *  - Function pointer to function that will be fed all the particle pairs.
 *  - Pointer to data that will be supplied to said function.
 */
void forEveryPairD(void (*f)(Particle *p1, Particle *p2, void *data), void *data);
void forEveryPair(void (*f)(Particle *p1, Particle *p2));

/* Runs over all connections of particles that are within the same box, or 
 * in adjacent boxes. The first connection is from particle a1 to a2, the 
 * second is from b1 to b2. It is guaranteed that
 *    a2 == getConnectedParticle(a1)
 *    b2 == getConnectedParticle(b1)
 */
void forEveryConnectionPairD(void (*f)(Particle *a1, Particle *a2,
		Particle *b1, Particle *b2, void *data), void *data);
void forEveryConnectionPair(void (*f)(Particle *a1, Particle *a2, 
		Particle *b1, Particle *b2));

/* Check whether internal structure is still consistent. If checkCorrectBox 
 * is true, then also check if all particles are in their correct boxes.
 * This check also does a forEveryPairCheck.
 * If checkConnections is true, it also performs a forEveryConnectionCheck. */
bool spgridSanityCheck(bool checkCorrectBox, bool checkConnections);

/* Test to see if we iterate over the correct number of pairs, and see if 
 * we don't give a single same particle as two constituents of a pair. Does 
 * not explicitly test if we do *all* pairs, nor that we don't do doubles.
 * Returns true if everything is OK, false otherwise. */
bool forEveryPairCheck(void);
/* Analogous to forEveryPairCheck() */
bool forEveryConnectionPairCheck(void);

/* Returns the shortest vector that points from v1 to v2, taking into 
 * account the periodic boundary conditions. 
 * Precondition: The given vectors are allowed to break out of the grid, 
 * but they must be within one 'world-size' of the grid (ie in [-L, 2L] if 
 * the grid is [0, L] in each dimension.)
 *
 * NOTE: We can't fully define these here with an __inline__ attribute in 
 * the header because we need the dimensions of the world, which we hide in 
 * the .c file.  This shouldn't be a performance hit if your compiler can 
 * do Link Time Optimization (LTO), though.
 * If not: you'd probably want to bring those variables into this header 
 * and move the implementation from the .c file to here. */
Vec3 nearestImageVector(Vec3 v1, Vec3 v2);
Vec3 nearestImageUnitVector(Vec3 v1, Vec3 v2);
double nearestImageDistance(Vec3 v1, Vec3 v2);
double nearestImageDistance2(Vec3 v1, Vec3 v2);
/* Like nearestImageVector, but also works if the vectors are off by more 
 * than a 'grid length'. */
Vec3 nearestImageVectorSafe(Vec3 v1, Vec3 v2);

#endif
