#include "spgrid.h"
#include "vmath.h"
#include <stdio.h>
#include <stdlib.h>

typedef struct box
{
	Particle *p;
	int n;
} Box;


static void addToBox(Particle *p, Box *b);
static void removeFromBox(Particle *p, Box *b);
static Box *boxFromIndex(int ix, int iy, int iz);
static Box *boxFromParticle(const Particle *p);
static Box *boxFromNonPeriodicIndex(int ix, int iy, int iz);
static void forEVERYpair(void (*f)(Particle *p1, Particle *p2, void *data), 
			 void *data);
static void pairWrapper(Particle *p1, Particle *p2, void *data);


/* Globals */
static Box *grid;
static double boxSize = 0; /* Linear length of one box. */
static int nb = 0; /* Number of Boxes in one dimension. Total number of 
		      boxes is nb^3. Total volume of the grid is 
		      (nb*boxSize)^3. */
static int numParticles = 0; /* Total number of particles in the grid. For 
				consistency checking only! */

bool allocGrid(int numBoxes, double size)
{
	assert(grid == NULL && nb == 0);
	grid = calloc(numBoxes * numBoxes * numBoxes, sizeof(*grid));
	if (grid == NULL)
		return false;
	nb = numBoxes;
	boxSize = size / numBoxes;
	assert(sanityCheck());
	return true;
}

void freeGrid()
{
	if (grid == NULL) {
		assert(nb == 0);
		return;
	}

	for (int i = 0; i < nb*nb*nb; i++) {
		Box *box = &grid[i];
		Particle *p = box->p;

		if (p == NULL)
			continue; /* Empty box */

		int n = box->n; /* Save original number, it will decrease 
				   during the loop as we remove particles! */
		for (int j = 0; j < n; j++) {
			Particle *next = p->next;
			removeFromBox(p, box);
			numParticles--;
			p = next;
		}
		assert(box->n == 0  &&  box->p == NULL);
	}
	assert(sanityCheck());
	assert(numParticles == 0);

	nb = 0;
	free(grid);
}

void addToGrid(Particle *p) {
	periodic(nb * boxSize, &p->pos, &p->pos);
	Box *box = boxFromParticle(p);
	addToBox(p, box);
	numParticles++;
	assert(sanityCheck());
}

void reboxParticles(void)
{
	double gridSize = nb * boxSize;

	for (int i = 0; i < nb*nb*nb; i++) {
		Box *currentBox = &grid[i];
		Particle *p = currentBox->p;

		if (p == NULL)
			continue; /* Empty box */

		int n = currentBox->n; /* Save original number, it can 
					  decrease during the loop if we 
					  swap out particles! */
		for (int j = 0; j < n; j++) {
			/* Force periodic boundary condition. */
			periodic(gridSize, &p->pos, &p->pos);

			/* Since p might be removed from the current box, 
			 * we keep a pointer to its successor. */
			Particle *next = p->next;

			Box *correctBox = boxFromParticle(p);
			if (currentBox != correctBox) {
				removeFromBox(p, currentBox);
				addToBox(p, correctBox);
			}

			p = next;
		}
	}
	assert(sanityCheck());
}

/* Precondition: particle must be within the grid. */
static Box *boxFromParticle(const Particle *p)
{
	int nx, ny, nz;

	assert(p != NULL);
	assert(!isnan(p->pos.x) && !isnan(p->pos.y) && !isnan(p->pos.z));
	assert(0 <= p->pos.x  &&  p->pos.x < boxSize * nb);
	assert(0 <= p->pos.y  &&  p->pos.y < boxSize * nb);
	assert(0 <= p->pos.z  &&  p->pos.z < boxSize * nb);

	nx = p->pos.x / boxSize;
	ny = p->pos.y / boxSize;
	nz = p->pos.z / boxSize;

	return boxFromIndex(nx, ny, nz);
}

static Box *boxFromNonPeriodicIndex(int ix, int iy, int iz)
{
	ix = ix % nb;
	if (ix < 0) ix += nb;

	iy = iy % nb;
	if (iy < 0) iy += nb;

	iz = iz % nb;
	if (iz < 0) iz += nb;

	return boxFromIndex (ix, iy, iz);
}

static Box *boxFromIndex(int ix, int iy, int iz)
{
	assert(0 <= ix && ix < nb);
	assert(0 <= iy && iy < nb);
	assert(0 <= iz && iz < nb);

	return grid + ix*nb*nb + iy*nb + iz;
}

static void removeFromBox(Particle *p, Box *b)
{
	assert(p != NULL);
	assert(b->n != 0);

	if (b->n == 1) {
		assert(p->prev == p);
		assert(p->next == p);
		b->p = NULL;
	} else {
		assert(p->prev->next == p);
		assert(p->next->prev == p);
		p->prev->next = p->next;
		p->next->prev = p->prev;

		if (b->p == p)
			b->p = p->next;
	}

	p->prev = NULL;
	p->next = NULL;

	b->n--;
}

static void addToBox(Particle *p, Box *b)
{
	assert(p->prev == NULL);
	assert(p->next == NULL);

	if (b->p == NULL) {
		assert(b->n == 0);
		b->p = p;
		p->prev = p;
		p->next = p;
	} else {
		assert(b->n > 0);
		p->next = b->p;
		p->prev = b->p->prev;
		p->prev->next = p;
		p->next->prev = p;
	}

	b->n++;
}

void forEveryPairD(void (*f)(Particle *p1, Particle *p2, void *data), void *data)
{
	int n1, n2;

	assert(sanityCheck());

	if (nb < 3) {
		/* Brute force. Reason: see comment below */
		forEVERYpair(f, data);
		return;
	}

	/* Loop over all boxes */
	for (int ix = 0; ix < nb; ix++)
	for (int iy = 0; iy < nb; iy++)
	for (int iz = 0; iz < nb; iz++) {
		/* Loop over all particles in this box */
		Box *box = boxFromIndex(ix, iy, iz);
		Particle *p = box->p;

		/* Loop over every partner of the i'th particle 'p' from the 
		 * box 'box' */
		n1 = box->n;
		for (int i = 0; i < n1; i++) { /* i'th particle in box */
			Particle *p2 = p->next;
			for (int j = i + 1; j < n1; j++) {
				(*f)(p, p2, data);
				p2 = p2->next;
			}
			/* Loop over particles in adjacent boxes to the box 
			 * of p. We need a total ordering on the boxes so 
			 * we don't check the same box twice. We use the 
			 * pointer value for this.
			 * However, due to periodic boundary conditions, 
			 * this ONLY works when there are AT LEAST 3 boxes 
			 * in each dimension! */
			for (int dix = -1; dix <= 1; dix++)
			for (int diy = -1; diy <= 1; diy++)
			for (int diz = -1; diz <= 1; diz++) {
				Box *b = boxFromNonPeriodicIndex(
						ix+dix, iy+diy, iz+diz);
				if (b <= box)
					continue;
					/* if b == box: it's our own box!
					 * else: only check boxes that have 
					 * a strictly larger pointer value 
					 * to avoid double work. */
				p2 = b->p;
				n2 = b->n;
				for (int j = 0; j < n2; j++) {
					(*f)(p, p2, data);
					p2 = p2->next;
				}
				assert(p2 == b->p);
			}

			p = p->next;
		}

		assert(p == box->p);
	}
}

void forEveryPair(void (*f)(Particle *p1, Particle *p2))
{
	/* I *hope* the compiler can optimize this deep chain of function 
	 * pointer magic. TODO: Check this and deal with the ISO C warnings 
	 * somehow. */
	forEveryPairD(&pairWrapper, (void*) f);
}

/* This is a bit of a hack, but it works and avoids code duplication. Ask 
 * the non-data function poiner as the data argument. */
static void pairWrapper(Particle *p1, Particle *p2, void *data)
{
	/* Black function pointer casting magic */
	void (*f)(Particle *p1, Particle *p2) = (void (*)(Particle *p1, Particle *p2)) data;
	(*f)(p1, p2);
}
	



/* Brute force over *every single* pair, including those that are more than 
 * a boxlength apart. */
static void forEVERYpair(void (*f)(Particle *p1, Particle *p2, void *data), 
			 void *data)
{

	for (int b1 = 0; b1 < nb*nb*nb; b1++) {
		Particle *p1 = grid[b1].p;
		for (int i1 = 0; i1 < grid[b1].n; i1++) {
			for (int b2 = 0; b2 < nb*nb*nb; b2++) {
				Particle *p2 = grid[b1].p;
				for (int i2 = 0; i2 < grid[b1].n; i2++) {
					(*f)(p1, p2, data);
					p2 = p2->next;
				}
			}
			p1 = p1->next;
		}
	}
}

bool sanityCheck(void)
{
	int nParts1 = 0;
	int nParts2 = 0;
	bool OK = true;

	/* Check linked list consistency */
	for (int i = 0; i < nb*nb*nb; i++) {
		Box *box = &grid[i];
		Particle *p = box->p;

		if (p == NULL)
			continue; /* Empty box */

		int n = box->n; 
		for (int j = 0; j < n; j++) {
			if (p->next->prev != p || p->prev->next != p) {
				fprintf(stderr, "%p is in a broken list\n",
						(const void *) p);
				OK = false;
			}
			p = p->next;
		}
	}

	/* 1) Check if each particle is in the box it should be in, given 
	 * its coordinates.
	 * 2) Count the number of particles and check them with
	 * the total we got when initially adding the particles. */
	for (int i = 0; i < nb * nb * nb; i++)
	{
		Box *b = &grid[i];
		const Particle *p, *first;

		if (b->p == NULL) {
			if (b->n != 0) {
				fprintf(stderr, "Box %d: found zero particles, "
						"expected %d\n", i, b->n);
				OK = false;
			}
			continue;
		}

		first = b->p;
		p = first;
		int j = 0;
		do {
			Box *correctBox = boxFromParticle(p);
			if (correctBox != b) {
				int c = (correctBox - grid)/sizeof(*correctBox);
				fprintf(stderr, "Particle is in box %d, "
					"should be in %ld\n", i, (correctBox -
					grid)/sizeof(*correctBox));
				fprintf(stderr, "numBox per dim: %d\n", nb);
				fprintf(stderr, "Pos:\t");
				fprintVector(stderr, &p->pos);
				fprintf(stderr, "\n");
				fprintf(stderr, "Actual box coords:  %d %d %d\n",
						i/nb/nb, (i/nb)%nb, i%nb);
				fprintf(stderr, "Correct box coords: %d %d %d\n",
						c/nb/nb, (c/nb)%nb, c%nb);
				OK = false;
			}
			j++;
			nParts1++;
			p = p->next;
		} while (p != first);

		if (j != b->n) {
			fprintf(stderr, "Box %d: found %d particles, "
					"expected %d\n", i, j, b->n);
			OK = false;
		}
		nParts2 += b->n;
	}

	if (nParts1 != numParticles)
	{
		fprintf(stderr, "1: Found a total of %d particles, "
			"should be %d\n", nParts1, numParticles);
		OK = false;
	}

	if (nParts2 != numParticles)
	{
		fprintf(stderr, "2: Found a total of %d particles, "
			"should be %d\n", nParts2, numParticles);
		OK = false;
	}

	return OK;
}

Vec3 nearestImageVector(Vec3 *v1, Vec3 *v2)
{
	Vec3 diff;
	double L = nb * boxSize;
	/* +5*L/2 instead of simply +L/2 to make sure that the first 
	 * argument of fmod is positive. Otherwise, this won't work! */
	diff.x = fmod(v2->x - v1->x + 5*L/2, L)  -  L/2;
	diff.y = fmod(v2->y - v1->y + 5*L/2, L)  -  L/2;
	diff.z = fmod(v2->z - v1->z + 5*L/2, L)  -  L/2;
	return diff;
}

double nearestImageDistance(Vec3 *v1, Vec3 *v2)
{
	Vec3 rijVec = nearestImageVector(v1, v2);
	double rijLength = length(&rijVec);
	
	return rijLength;
}


Vec3 nearestImageUnitVector(Vec3 *v1, Vec3 *v2)
{
	Vec3 rijUnit;
	
	Vec3 rijVec = nearestImageVector(v1, v2);
	double rijLength = length(&rijVec);
	
	scale(&rijVec, rijLength, &rijUnit);
	return rijUnit;
}


