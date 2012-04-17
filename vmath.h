#ifndef _VMATH_H_
#define _VMATH_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct Vec3
{
	double x, y, z;
} Vec3;

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

/* In header for inlining */

static __inline__ void printVector(const Vec3 *v);
static __inline__ void fprintVector(FILE *stream, const Vec3 *v);
static __inline__ void add(const Vec3 *a, const Vec3 *b, Vec3 *dest);
static __inline__ void sub(const Vec3 *a, const Vec3 *b, Vec3 *dest);
static __inline__ void scale(const Vec3 *v, double lambda, Vec3 *dest);
static __inline__ void normalize(const Vec3 *v, Vec3 *w);
static __inline__ double dot(const Vec3 *v, const Vec3 *w);
static __inline__ double length2(const Vec3 *v);
static __inline__ double length(const Vec3 *v);
static __inline__ double distance2(const Vec3 *a, const Vec3 *b);
static __inline__ double distance(const Vec3 *a, const Vec3 *b);

static __inline__ void fprintVector(FILE *stream, const Vec3 *v)
{
	fprintf(stream, "%10f\t%10f\t%10f\t", v->x, v->y, v->z);
}

static __inline__ void printVector(const Vec3 *v)
{
	printf("%10f\t%10f\t%10f\t", v->x, v->y, v->z);
}

static __inline__ void printVectorExp(const Vec3 *v)
{
	printf("%15e %15e %15e", v->x, v->y, v->z);
}

static __inline__ void add(const Vec3 *a, const Vec3 *b, Vec3 *dest)
{
	dest->x = a->x + b->x;
	dest->y = a->y + b->y;
	dest->z = a->z + b->z;
}

static __inline__ void sub(const Vec3 *a, const Vec3 *b, Vec3 *dest)
{
	dest->x = a->x - b->x;
	dest->y = a->y - b->y;
	dest->z = a->z - b->z;
}

static __inline__ void scale(const Vec3 *v, double lambda, Vec3 *dest)
{
	dest->x = lambda * v->x;
	dest->y = lambda * v->y;
	dest->z = lambda * v->z;
}

static __inline__ void normalize(const Vec3 *v, Vec3 *w)
{
	double l = length(v);
	scale(v, 1/l, w);
}

static __inline__ Vec3 cross(const Vec3 *v, const Vec3 *w)
{
	Vec3 x;
	x.x = v->y * w->z  -  v->z * w->y;
	x.y = v->z * w->x  -  v->x * w->z;
	x.z = v->x * w->y  -  v->y * w->x;

	assert(fabs(dot(v, &x) / length(v) / length(&x)) < 1e-10);
	assert(fabs(dot(w, &x) / length(w) / length(&x)) < 1e-10);

	return x;
}

static __inline__ double dot(const Vec3 *v, const Vec3 *w)
{
	return v->x * w->x + v->y * w->y + v->z * w->z;
}

static __inline__ double length2(const Vec3 *v)
{
	return dot(v, v);
}

static __inline__ double length(const Vec3 *v)
{
	return sqrt(length2(v));
}

static __inline__ double distance2(const Vec3 *a, const Vec3 *b)
{
	Vec3 c;

	sub(a, b, &c);
	return length2(&c);
}

static __inline__ double distance(const Vec3 *a, const Vec3 *b)
{
	return sqrt(distance2(a, b));
}

static __inline__ double cosAngle(const Vec3 *v, const Vec3 *w)
{
	return dot(v, w) / (length(v) * length (w));
}

static __inline__ double angle(const Vec3 *v, const Vec3 *w)
{
	return acos(cosAngle(v, w));
}

static __inline__ double dihedral(const Vec3 *v1, const Vec3 *v2, const Vec3 *v3)
{
	Vec3 v1xv2 = cross(v1, v2);
	Vec3 v2xv3 = cross(v2, v3);
	return atan2(length(v2) * dot(v1, &v2xv3), dot(&v1xv2, &v2xv3));
}

static __inline__ void periodic(double period, const Vec3 *v, Vec3 *dest)
{
	/* Fmod doesn't handle negative values the way we want it to, so we 
	 * need an extra check for the sign.
	 * TODO: see if we can do this without a branch! */
	dest->x = fmod(v->x, period);
	if (dest->x < 0) dest->x += period;
	dest->y = fmod(v->y, period);
	if (dest->y < 0) dest->y += period;
	dest->z = fmod(v->z, period);
	if (dest->z < 0) dest->z += period;
	
	assert(0 <= v->x  &&  v->x < period);
	assert(0 <= v->y  &&  v->y < period);
	assert(0 <= v->z  &&  v->z < period);
}

/* Returns a number sampled from a standard normal distribution. */
static __inline__ double randNorm(void)
{
	/* Box-Muller transform */
	double u1 = ((double) rand()) / RAND_MAX;
	double u2 = ((double) rand()) / RAND_MAX;

	return sqrt(-2*log(u1)) * cos(2*M_PI*u2);
}

/* Returns a vector with components sampled from a standard normal 
 * distribution. */
static __inline__ Vec3 randNormVec(void)
{
	Vec3 res;
	res.x = randNorm();
	res.y = randNorm();
	res.z = randNorm();
	return res;
}

#endif
