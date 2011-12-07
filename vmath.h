#ifndef _VMATH_H_
#define _VMATH_H_

#include <stdio.h>
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

#endif
