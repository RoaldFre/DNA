#ifndef _SPGRID_H_
#define _SPGRID_H_

/* SPgrid: Space Partition grid */

#include "system.h" //TODO separate in particle.h (= all we need)?

/* Allocates a grid [0, size] x [0, size] x [0, size] with 'numBoxes' boxes 
 * in each dimension.
 * Precondition: grid can't already be allocated (unless it was freed 
 * afterwards).
 * Returns true on succes, false on failure. */
bool allocGrid(int numBoxes, double size);

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
 * case the particles escaped from the grid. */
void reboxParticles(void);

/* Execute a given function for every discinct pair of particles that are 
 * within the same box, or in adjacent boxes (taking into account periodic 
 * boundary conditions).
 * Arguments:
 *  - Function pointer to function that will be fed all the particle pairs.
 *  - Pointer to data that will be supplied to said function.
 */
void forEveryPair(void (*f)(Particle *p1, Particle *p2, void *data), void *data);

/* Check whether internal structure is still consistent. */
bool sanityCheck(void);

#endif
