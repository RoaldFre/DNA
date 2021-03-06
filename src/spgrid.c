#include "spgrid.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>

/* !! NOTE !!
 * Currently, particles can have a position at (or beyond!) the 'last' box, 
 * yet be assigned to the 'first' box in a dimension. This is to assure 
 * that we don't access the wrong memory, in case of numerical rounding 
 * problems!
 * The result: we can't count on the box indices to know where a particle 
 * is. Hence, always use nearestImageVector() and friends.
 * Also, all asserts regarding numerical accuracy have been commented.
 * TODO: if particles are guaranteed to be within the boxes (with all the 
 * numerical problems that this brings in), then we can quickly do a 
 * nearestImageVector() based on the box indices! Maybe that's faster than 
 * the current way. (we can pass the difference vector to functions in a 
 * forEveryPair() loop...) */

struct box
{
	Particle *p;	/* (Linked list of) particles in this box */
	int n;		/* number of particles in this box */

	/* Circular linked list of boxes that contain particles. Both are 
	 * NULL if this box has no particles, and hence is not part of the 
	 * list. */
	struct box *prevOccupied;
	struct box *nextOccupied;

	/* TODO be more smart so we don't need this */
	struct box *prevX;
	struct box *nextX;
	struct box *prevY;
	struct box *nextY;
	struct box *prevZ;
	struct box *nextZ;
};
typedef struct box Box;


static void addToBox(Particle *p, Box *b);
static void removeFromBox(Particle *p, Box *b);
static Box *boxFromIndex(int ix, int iy, int iz);
static Box *boxFromNonPeriodicParticle(const Particle *p);
static Box *boxFromNonPeriodicIndex(int ix, int iy, int iz);
static void forEveryPairBruteForce(void (*f)(Particle *p1, Particle *p2, void *data), 
			 void *data);


/* Globals */
static Box *grid; /* All boxes in the spgrid. */
static Box *occupiedBoxes; /* First element in linked list of boxes that 
			      contain particles, or NULL if all boxes are 
			      empty. */
static int nb = 0; /* Number of Boxes in one dimension. Total number of 
		      boxes is nb^3. Total volume of the grid is 
		      (nb*boxSize)^3. */
static int gridNumParticles = 0; /* Total number of particles in the grid. For 
				consistency checking only! */

/* Add the (newly) occupied box to the list. It cannot already be part of 
 * the list. */
static void addOccupiedBox(Box *box)
{
	assert(box != NULL);
	assert(box->n > 0);
	assert(box->p != NULL);
	assert(box->nextOccupied == NULL);
	assert(box->prevOccupied == NULL);

	if (occupiedBoxes == NULL) {
		/* Start a new list */
		box->nextOccupied = box;
		box->prevOccupied = box;
		occupiedBoxes = box;
	} else {
		box->nextOccupied = occupiedBoxes;
		box->prevOccupied = occupiedBoxes->prevOccupied;
		box->prevOccupied->nextOccupied = box;
		box->nextOccupied->prevOccupied = box;
	}
}
/* Remove the empty box from the list */
static void removeNonOccupiedBox(Box *emptyBox)
{
	assert(emptyBox != NULL);
	assert(occupiedBoxes != NULL);
	assert(emptyBox->n == 0);
	assert(emptyBox->p == NULL);

	if (emptyBox->nextOccupied == emptyBox) {
		/* List contains only emptyBox itself */
		assert(emptyBox->prevOccupied == emptyBox);
		assert(emptyBox == occupiedBoxes);

		occupiedBoxes = NULL;
	} else {
		/* List contains at least one other box */
		assert(emptyBox->prevOccupied->nextOccupied == emptyBox);
		assert(emptyBox->nextOccupied->prevOccupied == emptyBox);

		emptyBox->prevOccupied->nextOccupied = emptyBox->nextOccupied;
		emptyBox->nextOccupied->prevOccupied = emptyBox->prevOccupied;

		/* If the list uses emptyBox as head, pick another head */
		if (occupiedBoxes == emptyBox)
			occupiedBoxes = emptyBox->nextOccupied;
	}

	emptyBox->prevOccupied = NULL;
	emptyBox->nextOccupied = NULL;
}


void allocGrid(int numBoxes)
{
	assert(grid == NULL && nb == 0);
	grid = calloc(numBoxes * numBoxes * numBoxes, sizeof(*grid));
	if (grid == NULL)
		die("Couldn't allocate memory for space partitioning grid!\n");
	nb = numBoxes;

	/* set the prev/nextXYZ pointers */
	for (int ix = 0; ix < nb; ix++)
	for (int iy = 0; iy < nb; iy++)
	for (int iz = 0; iz < nb; iz++) {
		Box *box = boxFromIndex(ix, iy, iz);

		box->nextX = boxFromNonPeriodicIndex(ix+1, iy,   iz  );
		box->prevX = boxFromNonPeriodicIndex(ix-1, iy,   iz  );
		box->nextY = boxFromNonPeriodicIndex(ix,   iy+1, iz  );
		box->prevY = boxFromNonPeriodicIndex(ix,   iy-1, iz  );
		box->nextZ = boxFromNonPeriodicIndex(ix,   iy,   iz+1);
		box->prevZ = boxFromNonPeriodicIndex(ix,   iy,   iz-1);
	}

	assert(spgridSanityCheck(true, false));
}

void initGrid(int numBoxes)
{
	allocGrid(numBoxes);
	forEveryParticle(&addToGrid);
}

void freeGrid()
{
	if (grid == NULL) {
		assert(nb == 0);
		assert(occupiedBoxes == NULL);
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

	occupiedBoxes = NULL;

	nb = 0;
	free(grid);
}

void addToGrid(Particle *p) {
	p->pos = periodic(world.worldSize, p->pos);
	Box *box = boxFromNonPeriodicParticle(p);
	addToBox(p, box);
	gridNumParticles++;

	assert(spgridSanityCheck(false, false));
}

void reboxParticle(Particle *p)
{
	assert(p->myBox != NULL);

	periodicPosition(p);

	Box *correctBox = boxFromNonPeriodicParticle(p);
	if (correctBox == p->myBox)
		return;

	removeFromBox(p, p->myBox);
	addToBox(p, correctBox);
}
void reboxParticles(void)
{
	if (nb == 1) {
		forEveryParticle(&periodicPosition);
		return;
	}

	assert(spgridSanityCheck(false, true));

	forEveryParticle(&reboxParticle);

	assert(spgridSanityCheck(true, true));
}

/* Particle may be outside the grid */
static Box *boxFromNonPeriodicParticle(const Particle *p)
{
	assert(p != NULL);
	assert(!isnan(p->pos.x) && !isnan(p->pos.y) && !isnan(p->pos.z));

	double ws = world.worldSize;
	double boxSize = ws / nb;
	Vec3 shifted = add(p->pos, (Vec3) {ws/2.0, ws/2.0, ws/2.0});

	int ix = shifted.x / boxSize;
	int iy = shifted.y / boxSize;
	int iz = shifted.z / boxSize;

	return boxFromNonPeriodicIndex(ix, iy, iz);
}

static Box *boxFromNonPeriodicIndex(int ix, int iy, int iz)
{
	ix = ix % nb;
	if (UNLIKELY(ix < 0)) ix += nb;

	iy = iy % nb;
	if (UNLIKELY(iy < 0)) iy += nb;

	iz = iz % nb;
	if (UNLIKELY(iz < 0)) iz += nb;

	return boxFromIndex(ix, iy, iz);
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
	assert(p != NULL && b != NULL);
	assert(p->myBox == b);
	assert(b->n >= 0);

	b->n--;

	if (b->n == 0) {
		assert(p->prev == p);
		assert(p->next == p);
		b->p = NULL;
		removeNonOccupiedBox(b);
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
	p->myBox = NULL;
}

static void addToBox(Particle *p, Box *b)
{
	assert(p->prev == NULL);
	assert(p->next == NULL);
	assert(p->myBox == NULL);

	b->n++;

	if (b->p == NULL) {
		assert(b->n == 1);
		b->p = p;
		p->prev = p;
		p->next = p;
		addOccupiedBox(b);
	} else {
		assert(b->n > 1);
		p->next = b->p;
		p->prev = b->p->prev;
		p->prev->next = p;
		p->next->prev = p;
	}
	p->myBox = b;
}






/* ITERATIONS OVER NEIGHBOURS OF A GIVEN PARTICLE */
/* Loop over all pairs (p,q), where q is a particle in the given box, and p 
 * is the given particle. The case p==q is ignored. */
static void matchAgainstBox(Particle *p, Box *box,
		void (*f)(Particle *p1, Particle *p2, void *data), void *data)
{
	if (box->n == 0)
		return;

	int n = box->n;
	Particle *q = box->p;
	for (int j = 0; j < n; j++) {
		if (p != q)
			f(p, q, data);
		q = q->next;
	}
	assert(q == box->p);
}

typedef struct {
	Particle *p;
	void (*f)(Particle *p1, Particle *p2, void *data);
	void *data;
} FenbfData;
static void forEveryNeighBruteForceWrapper(Particle *p, void *data)
{
	FenbfData *bfData = (FenbfData*) data;
	bfData->f(p, bfData->p, bfData->data);
}

/* Visit all neighbouring particles of the given particle. */
void forEveryNeighbourOf(Particle *p,
		void (*f)(Particle *p1, Particle *p2, void *data), void *data)
{
	if (nb < 3) {
		/* Reason: see forEveryPairD() */
		/* Anonimous lambda functions would be really nice here */
		FenbfData bfData;
		bfData.p = p;
		bfData.f = f;
		bfData.data = data;
		forEveryParticleD(&forEveryNeighBruteForceWrapper, &bfData);
	}

	Box *box = p->myBox;

	/* x-1 */
	matchAgainstBox(p, box->prevX->prevY->nextZ, f, data);
	matchAgainstBox(p, box->prevX->prevY,        f, data);
	matchAgainstBox(p, box->prevX->prevY->prevZ, f, data);

	matchAgainstBox(p, box->prevX->nextZ,        f, data);
	matchAgainstBox(p, box->prevX,               f, data);
	matchAgainstBox(p, box->prevX->prevZ,        f, data);

	matchAgainstBox(p, box->prevX->nextY->nextZ, f, data);
	matchAgainstBox(p, box->prevX->nextY,        f, data);
	matchAgainstBox(p, box->prevX->nextY->prevZ, f, data);


	/* x */
	matchAgainstBox(p, box->prevY->nextZ, f, data);
	matchAgainstBox(p, box->prevY,        f, data);
	matchAgainstBox(p, box->prevY->prevZ, f, data);

	matchAgainstBox(p, box->nextZ,        f, data);
	matchAgainstBox(p, box,               f, data);
	matchAgainstBox(p, box->prevZ,        f, data);

	matchAgainstBox(p, box->nextY->nextZ, f, data);
	matchAgainstBox(p, box->nextY,        f, data);
	matchAgainstBox(p, box->nextY->prevZ, f, data);


	/* x+1 */
	matchAgainstBox(p, box->nextX->prevY->nextZ, f, data);
	matchAgainstBox(p, box->nextX->prevY,        f, data);
	matchAgainstBox(p, box->nextX->prevY->prevZ, f, data);

	matchAgainstBox(p, box->nextX->nextZ,        f, data);
	matchAgainstBox(p, box->nextX,               f, data);
	matchAgainstBox(p, box->nextX->prevZ,        f, data);

	matchAgainstBox(p, box->nextX->nextY->nextZ, f, data);
	matchAgainstBox(p, box->nextX->nextY,        f, data);
	matchAgainstBox(p, box->nextX->nextY->prevZ, f, data);
}





/* Loop between particles of adjacent boxes of box.  We need a total 
 * ordering on the boxes so we don't check the same box twice. We use the 
 * pointer value for this.
 * However, due to periodic boundary conditions, this ONLY works when there 
 * are AT LEAST 3 boxes in each dimension!  */
static void visitNeighbourPairs(Box *box, Box *neighbour,
		void (*f)(Particle *p1, Particle *p2, void *data), void *data)
{
	if (neighbour <= box)
		return;
		/* if neighbour == box: it's our own box!
		 * else: only check boxes that have a 
		 * strictly larger pointer value to avoid 
		 * double work. */

	if (neighbour->n == 0)
		return;

	/* Match up all particles from box and neighbour. */
	int n1 = box->n;
	int n2 = neighbour->n;
	Particle *p1 = box->p;
	for (int i = 0; i < n1; i++) {
		Particle *p2 = neighbour->p;
		for (int j = 0; j < n2; j++) {
			f(p1, p2, data);
			p2 = p2->next;
		}
		assert(p2 == neighbour->p);
		p1 = p1->next;
	}
	assert(p1 == box->p);
}

static void visitNeigbourPairsOf(Box *box,
		void (*f)(Particle *p1, Particle *p2, void *data), void *data)
{
	//TODO: be more smart/elegant

	/* x-1 */
	visitNeighbourPairs(box, box->prevX->prevY->nextZ, f, data);
	visitNeighbourPairs(box, box->prevX->prevY,        f, data);
	visitNeighbourPairs(box, box->prevX->prevY->prevZ, f, data);

	visitNeighbourPairs(box, box->prevX->nextZ,        f, data);
	visitNeighbourPairs(box, box->prevX,               f, data);
	visitNeighbourPairs(box, box->prevX->prevZ,        f, data);

	visitNeighbourPairs(box, box->prevX->nextY->nextZ, f, data);
	visitNeighbourPairs(box, box->prevX->nextY,        f, data);
	visitNeighbourPairs(box, box->prevX->nextY->prevZ, f, data);


	/* x */
	visitNeighbourPairs(box, box->prevY->nextZ, f, data);
	visitNeighbourPairs(box, box->prevY,        f, data);
	visitNeighbourPairs(box, box->prevY->prevZ, f, data);

	visitNeighbourPairs(box, box->nextZ,        f, data);
	visitNeighbourPairs(box, box->prevZ,        f, data);

	visitNeighbourPairs(box, box->nextY->nextZ, f, data);
	visitNeighbourPairs(box, box->nextY,        f, data);
	visitNeighbourPairs(box, box->nextY->prevZ, f, data);


	/* x+1 */
	visitNeighbourPairs(box, box->nextX->prevY->nextZ, f, data);
	visitNeighbourPairs(box, box->nextX->prevY,        f, data);
	visitNeighbourPairs(box, box->nextX->prevY->prevZ, f, data);

	visitNeighbourPairs(box, box->nextX->nextZ,        f, data);
	visitNeighbourPairs(box, box->nextX,               f, data);
	visitNeighbourPairs(box, box->nextX->prevZ,        f, data);

	visitNeighbourPairs(box, box->nextX->nextY->nextZ, f, data);
	visitNeighbourPairs(box, box->nextX->nextY,        f, data);
	visitNeighbourPairs(box, box->nextX->nextY->prevZ, f, data);
}

/* ITERATION OVER PAIRS */
void forEveryPairD(void (*f)(Particle *p1, Particle *p2, void *data), void *data)
{
	if (nb < 3) {
		/* Brute force. We can't do anything better anyway with 2 
		 * boxes or less per dimension, as all particles will have 
		 * to be checked against each other.
		 * If we use the method below (with visitNeighbourPairs()), 
		 * we get into trouble when nb==2. In that case, the 'prev' 
		 * and the 'next' per dimension will map to the same box, 
		 * which would lead to double visits per pair! */
		forEveryPairBruteForce(f, data);
		return;
	}

	/* Loop over all occupied boxes */
	if (occupiedBoxes == NULL)
		return;

	Box *box = occupiedBoxes;
	do {
		/* Loop over all i'th particles 'p' from the box 'box' and 
		 * match them with the j'th praticle in the same box */
		Particle *p = box->p;
		int n = box->n;
		for (int i = 0; i < n; i++) {
			Particle *p2 = p->next;
			for (int j = i + 1; j < n; j++) {
				f(p, p2, data);
				p2 = p2->next;
			}

			p = p->next;
		}
		assert(p == box->p); /* We went 'full circle' */

		visitNeigbourPairsOf(box, f, data);

		box = box->nextOccupied;
	} while (box != occupiedBoxes);
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
 * a boxlength apart.
 * TODO: Possible improvement: Only loop over occupiedBoxes. But we really 
 * only call this in case of 1 or 2 boxes per dimension atm, so it's not 
 * really much of a speed-up then. */
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
	return fastPeriodic(world.worldSize, sub(v2, v1));
}
Vec3 nearestImageVectorSafe(Vec3 v1, Vec3 v2)
{
	//return periodic(world.worldSize, sub(v2, v1));
	return closePeriodic(world.worldSize, sub(v2, v1));
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
		fprintf(stderr, "number of particles in grid: %d\n", 
				gridNumParticles);
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

	/* BOXES */

	/* Check linked list consistency of boxes */
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



	/* OCCUPIED BOXES */

	/* Check linked list consistency of occupiedboxes and count them */
	int numOccupiedBoxes = 0;
	if (occupiedBoxes != NULL) {
		Box *box = occupiedBoxes;
		do {
			if (box->nextOccupied->prevOccupied != box
					|| box->prevOccupied->nextOccupied != box) {
				fprintf(stderr, "numOccupiedBoxes is broken at %p\n",
						(const void *) box);
				OK = false;
			}
			numOccupiedBoxes++;
			if (numOccupiedBoxes > nb*nb*nb) {
				fprintf(stderr, "numOccupiedBoxes contains more boxes "
						"than there are boxes in the "
						"grid!\n");
				OK = false;
				break; /* Avoid possible infinite loops */
			}
			box = box->nextOccupied;
		} while (box != occupiedBoxes);
	}

	/* Check whether all occupied boxes (n > 0) are in the list, and 
	 * that we don't have any extra in the list (count them) */
	int correctNumOccupiedBoxes = 0;
	for (int ix = 0; ix < nb; ix++)
	for (int iy = 0; iy < nb; iy++)
	for (int iz = 0; iz < nb; iz++) {
		Box *box = boxFromIndex(ix, iy, iz);

		if (box->n == 0) {
			/* empty box */
			if (box->nextOccupied != NULL
					|| box->prevOccupied != NULL) {
				fprintf(stderr, "box %p has no particles "
						"but appears to be in the "
						"occupiedBoxes list!\n",
						(void *) box);
				OK = false;
			}
		} else {
			/* occupied box */
			correctNumOccupiedBoxes++;

			if (box->nextOccupied == NULL
					|| box->prevOccupied == NULL) {
				fprintf(stderr, "box %p has particles but "
						"doesn't appear to be in the "
						"occupiedBoxes list!\n",
						(void *) box);
				OK = false;
			}
		}
	}

	if (numOccupiedBoxes != correctNumOccupiedBoxes) {
		fprintf(stderr, "Found %d boxes in the numOccupiedBoxes "
				"list, but counted %d occupied boxes!\n", 
				numOccupiedBoxes, correctNumOccupiedBoxes);
		OK = false;
	}



	/* PAIRS AND CONNECTIONS */

	OK = forEveryPairCheck() && OK;

	if (checkConnections)
		OK = forEveryConnectionPairCheck() && OK;

	return OK;
}
