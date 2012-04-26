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
static void forEveryPairBruteForce(void (*f)(Particle *p1, Particle *p2, void *data), 
			 void *data);


/* Globals */
static Box *grid;
static double boxSize = 0; /* Linear length of one box. */
static int nb = 0; /* Number of Boxes in one dimension. Total number of 
		      boxes is nb^3. Total volume of the grid is 
		      (nb*boxSize)^3. */
static int gridNumParticles = 0; /* Total number of particles in the grid. For 
				consistency checking only! */


bool allocGrid(int numBoxes, double size)
{
	assert(grid == NULL && nb == 0);
	grid = calloc(numBoxes * numBoxes * numBoxes, sizeof(*grid));
	if (grid == NULL)
		return false;
	nb = numBoxes;
	boxSize = size / numBoxes;
	assert(spgridSanityCheck(true, false));
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
			gridNumParticles--;
			p = next;
		}
		assert(box->n == 0  &&  box->p == NULL);
	}
	assert(spgridSanityCheck(true, true));
	assert(gridNumParticles == 0);

	nb = 0;
	free(grid);
}

void addToGrid(Particle *p) {
	p->pos = periodic(nb * boxSize, p->pos);
	Box *box = boxFromParticle(p);
	addToBox(p, box);
	gridNumParticles++;

	assert(spgridSanityCheck(false, false));
}

void reboxParticles(void)
{
	double gridSize = nb * boxSize;

	assert(spgridSanityCheck(false, true));

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
			p->pos = periodic(gridSize, p->pos);

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
	assert(spgridSanityCheck(true, true));
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
/* Precondition: particle must be within one worldlength distance from 
 * grid. */
static Box *boxFromNonPeriodicParticle(const Particle *p)
{
	int nx, ny, nz;

	assert(p != NULL);
	assert(!isnan(p->pos.x) && !isnan(p->pos.y) && !isnan(p->pos.z));

	nx = p->pos.x / boxSize;
	ny = p->pos.y / boxSize;
	nz = p->pos.z / boxSize;

	return boxFromNonPeriodicIndex(nx, ny, nz);
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



/* ITERATION OVER PAIRS */
void forEveryPairD(void (*f)(Particle *p1, Particle *p2, void *data), void *data)
{
	int n1, n2;

	if (nb < 3) {
		/* Brute force. Reason: see comment below */
		forEveryPairBruteForce(f, data);
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
				f(p, p2, data);
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
					f(p, p2, data);
					p2 = p2->next;
				}
				assert(p2 == b->p);
			}

			p = p->next;
		}

		assert(p == box->p);
	}
}

static void pairWrapper(Particle *p1, Particle *p2, void *data)
{
	void (**f)(Particle *p1, Particle *p2) =
			(void (**)(Particle *p1, Particle *p2)) data;
	(*f)(p1, p2);
}
void forEveryPair(void (*f)(Particle *p1, Particle *p2))
{
	/* I *hope* the compiler can optimize this deep chain of (function) 
	 * pointer magic. TODO: Check this! */
	forEveryPairD(&pairWrapper, (void*) &f);
}

/* Brute force over *every single* pair, including those that are more than 
 * a boxlength apart. */
static void forEveryPairBruteForce(void (*f)(Particle*, Particle*, void*), 
			 void *data)
{
	/* Pairs within the same box */
	for (int b = 0; b < nb*nb*nb; b++) {
		Particle *p1 = grid[b].p;
		for (int i = 0; i < grid[b].n; i++) { /* i'th particle in box */
			Particle *p2 = p1->next;
			for (int j = i + 1; j < grid[b].n; j++) {
				f(p1, p2, data);
				p2 = p2->next;
			}
			p1 = p1->next;
		}
	}

	/* Pairs in different boxes */
	for (int b1 = 0; b1 < nb*nb*nb; b1++) {
		Particle *p1 = grid[b1].p;
		for (int i1 = 0; i1 < grid[b1].n; i1++) {
			for (int b2 = b1 + 1; b2 < nb*nb*nb; b2++) {
				Particle *p2 = grid[b2].p;
				for (int i2 = 0; i2 < grid[b2].n; i2++) {
					f(p1, p2, data);
					p2 = p2->next;
				}
			}
			p1 = p1->next;
		}
	}
}




/* ITERATION OVER CONNECTION PAIRS */
typedef struct {
	void (*f)(Particle*, Particle*, Particle*, Particle*, void*);
	void *data;
} ForEveryConnectionPairData;
static void forEveryConnectionPairHelper(Particle *p1, Particle *p2, void *data)
{
	ForEveryConnectionPairData *fecpd = (ForEveryConnectionPairData*) data;

	Particle *q1 = getConnectedParticle(p1);
	Particle *q2 = getConnectedParticle(p2);
	if (q1 == NULL  ||  q2 == NULL)
		return;

	fecpd->f(p1, q1, p2, q2, fecpd->data);
}

void forEveryConnectionPairD(void (*f)(Particle *p1, Particle *p2,
		Particle *p3, Particle *p4, void *data), void *data)
{
	ForEveryConnectionPairData fecpd;
	fecpd.f = f;
	fecpd.data = data;
	forEveryPairD(&forEveryConnectionPairHelper, &fecpd);
}

static void connectionPairWrapper(Particle *p1, Particle *p2,
		Particle *p3, Particle *p4, void *data)
{
	void (**f)(Particle*, Particle*, Particle*, Particle*) =
			(void (**)(Particle*, Particle*, Particle*, Particle*)) data;
	(*f)(p1, p2, p3, p4);
}
void forEveryConnectionPair(void (*f)(Particle*, Particle*, Particle*, Particle*))
{
	/* I *hope* the compiler can optimize this deep chain of (function) 
	 * pointer magic. TODO: Check this! */
	forEveryConnectionPairD(&connectionPairWrapper, (void*) &f);
}




/* PERIODIC VECTOR FUNCTIONS */
Vec3 nearestImageVector(Vec3 v1, Vec3 v2)
{
	Vec3 res;
	double L = nb * boxSize;
	/* +5*L/2 instead of simply +L/2 to make sure that the first 
	 * argument of fmod is positive. Otherwise, this won't work! */
	res.x = fmod(v2.x - v1.x + 5*L/2, L)  -  L/2;
	res.y = fmod(v2.y - v1.y + 5*L/2, L)  -  L/2;
	res.z = fmod(v2.z - v1.z + 5*L/2, L)  -  L/2;
	return res;
}

double nearestImageDistance(Vec3 v1, Vec3 v2)
{
	return length(nearestImageVector(v1, v2));
}
double nearestImageDistance2(Vec3 v1, Vec3 v2)
{
	return length2(nearestImageVector(v1, v2));
}
Vec3 nearestImageUnitVector(Vec3 v1, Vec3 v2)
{
	return normalize(nearestImageVector(v1, v2));
}










/* TEST ROUTINES */

typedef struct
{
	int count;
	bool error;
} ForEveryCheckData;
static void forEveryPairCheckHelper(Particle *p1, Particle *p2, void *data)
{
	ForEveryCheckData *fecd = (ForEveryCheckData*) data;
	fecd->count++;
	if (p1 == p2) {
		fecd->error = true;
		fprintf(stderr, "forEveryPair gave illegal pair with "
				"particle %p\n", (void*)p1);
	}
}
bool forEveryPairCheck(void)
{
	ForEveryCheckData data;
	data.count = 0;
	data.error = false;

	forEveryPairD(&forEveryPairCheckHelper, &data);

	int correctCount = 0;
	if (nb < 3) {
		int n = gridNumParticles;
		correctCount = n * (n - 1) / 2;
	} else {
		for (int ix = 0; ix < nb; ix++)
		for (int iy = 0; iy < nb; iy++)
		for (int iz = 0; iz < nb; iz++) {
			Box *box = boxFromIndex(ix, iy, iz);
			/* Pairs in this box */
			int n1 = box->n;
			correctCount += n1 * (n1 - 1) / 2;

			/* Loop over pair with adjacent boxes to the box 
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
					 * to avoid double counting. */
				int n2 = b->n;
				correctCount += n1 * n2;
			}
		}
	}

	if (data.count != correctCount) {
		fprintf(stderr, "forEveryPair ran over %d pairs, but should "
				"be %d\n", data.count, correctCount);
		return false;
	}

	return !data.error;
}

static void forEveryConnectionPairCheckHelper(Particle *a1, Particle *a2, 
		Particle *b1, Particle *b2, void *data)
{
	ForEveryCheckData *fecd = (ForEveryCheckData*) data;
	fecd->count++;
	if (a1 == NULL || a2 == NULL || b1 == NULL || b2 == NULL
			|| a1 == a2 || b1 == b2
			|| (a1 == b1 && a2 == b2)
			|| (a1 == b2 && a2 == b1)
			|| (a2 != getConnectedParticle(a1))
			|| (b2 != getConnectedParticle(b1))) {
		fecd->error = false;
		fprintf(stderr, "forEveryConnectionPair gave illegal "
				"connection pair! %p %p %p %p\n",
				(void*)a1, (void*)a2, (void*)b1, (void*)b2);
	}
}
bool forEveryConnectionPairCheck(void)
{
	ForEveryCheckData data;
	data.count = 0;
	data.error = false;

	forEveryConnectionPairD(&forEveryConnectionPairCheckHelper, &data);

	//TODO determine the correct count to check.
	//printf("Connection pairs: %d\n", data.count);

	return !data.error;
}




bool spgridSanityCheck(bool checkCorrectBox, bool checkConnections)
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
			Box *correctBox = boxFromNonPeriodicParticle(p);
			if (checkCorrectBox && correctBox != b) {
				int c = (correctBox - grid)/sizeof(*correctBox);
				fprintf(stderr, "Particle is in box %d, "
					"should be in %ld\n", i, (correctBox -
					grid)/sizeof(*correctBox));
				fprintf(stderr, "numBox per dim: %d\n", nb);
				fprintf(stderr, "Pos:\t");
				fprintVector(stderr, p->pos);
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

	if (nParts1 != gridNumParticles)
	{
		fprintf(stderr, "1: Found a total of %d particles, "
			"should be %d\n", nParts1, gridNumParticles);
		OK = false;
	}

	if (nParts2 != gridNumParticles)
	{
		fprintf(stderr, "2: Found a total of %d particles, "
			"should be %d\n", nParts2, gridNumParticles);
		OK = false;
	}

	OK = forEveryPairCheck() && OK;

	if (checkConnections)
		OK = forEveryConnectionPairCheck() && OK;

	return OK;
}
