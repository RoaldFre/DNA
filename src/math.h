#ifndef _VMATH_H_
#define _VMATH_H_

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "tinymt/tinymt64.h"

/* tiny MT state */
extern tinymt64_t tinymt;
void seedRandomWith(uint64_t seed);
/* Automatically seeds random number generator based on current time and 
 * PID of process */
void seedRandom(void);

typedef struct Vec3
{
	double x, y, z;
} Vec3;

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

/* Warning, you can't 'call' the macros below with functions that have side 
 * effects as arguments!
 * Ie, don't do "ABS(functionWithSideEffectsThatReturnsADouble())"
 * In fact, just never use these with any function calls! */
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABS(x) ((x) > 0 ? (x) : -(x))
#define SQUARE(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))

/* For micro optimizations. Remove these if you aren't compiling with gcc 
 * and your compiler doesn't support them. */
#define LIKELY(x)       __builtin_expect((x),1)
#define UNLIKELY(x)     __builtin_expect((x),0)

/* In header for inlining */

static __inline__ bool isSaneNumber(double x)
{
	return !isnan(x) && !isinf(x);
}

static __inline__ bool isSaneVector(Vec3 v)
{
	return isSaneNumber(v.x) && isSaneNumber(v.y) && isSaneNumber(v.z);
}

/* Enable for debugging purposes. If an interaction generates an invalid 
 * vector, we will trigger a segfault. Only usefull if you run the code 
 * from a debugger or enable core dumps.
 *
 * This is useful for rare bugs because it only checks for vector sanity, 
 * whereas compiling with assertions checks all assertions and is therefore 
 * slower. */
#define DEBUG_VECTOR_SANITY true

static __inline__ void debugVectorSanity(Vec3 v, const char *location)
{
	if (!DEBUG_VECTOR_SANITY)
		return;

	if (isSaneVector(v))
		return;

	fprintf(stderr, "Found invalid vector at '%s'!\n"
			"Triggering segfault!\n", location);
	int *nil = (int*)NULL;
	*nil = 1; /* segfaults */
}

/* Check for equality of doubles (up to some small error). */
static __inline__ bool equalsEpsilon(double a, double b)
{
	if (a + b == 0)
		return a == 0;

	return fabs((a-b) / (a+b)) < 1e-5;
}

static __inline__ void fprintVector(FILE *stream, Vec3 v)
{
	fprintf(stream, "%10f\t%10f\t%10f\t", v.x, v.y, v.z);
}

static __inline__ void printVector(Vec3 v)
{
	fprintVector(stdout, v);
}

static __inline__ void printVectorExp(Vec3 v)
{
	printf("%15e %15e %15e", v.x, v.y, v.z);
}

static __inline__ Vec3 add(Vec3 a, Vec3 b)
{
	Vec3 res;
	res.x = a.x + b.x;
	res.y = a.y + b.y;
	res.z = a.z + b.z;
	return res;
}

static __inline__ Vec3 sub(Vec3 a, Vec3 b)
{
	Vec3 res;
	res.x = a.x - b.x;
	res.y = a.y - b.y;
	res.z = a.z - b.z;
	return res;
}

static __inline__ Vec3 scale(Vec3 v, double lambda)
{
	Vec3 res;
	res.x = lambda * v.x;
	res.y = lambda * v.y;
	res.z = lambda * v.z;
	return res;
}

static __inline__ double dot(Vec3 v, Vec3 w)
{
	return v.x * w.x + v.y * w.y + v.z * w.z;
}

static __inline__ double length2(Vec3 v)
{
	return dot(v, v);
}

static __inline__ double length(Vec3 v)
{
	return sqrt(length2(v));
}

static __inline__ Vec3 cross(Vec3 v, Vec3 w)
{
	Vec3 x;
	x.x = v.y * w.z  -  v.z * w.y;
	x.y = v.z * w.x  -  v.x * w.z;
	x.z = v.x * w.y  -  v.y * w.x;

	//assert(fabs(dot(v, x) / length(v) / length(x)) < 1e-5);
	//assert(fabs(dot(w, x) / length(w) / length(x)) < 1e-5);

	return x;
}

static __inline__ Vec3 normalize(Vec3 v)
{
	double l = length(v);
	assert(l != 0);
	return scale(v, 1/l);
}

static __inline__ double distance2(Vec3 a, Vec3 b)
{
	return(length2(sub(a, b)));
}

static __inline__ double distance(Vec3 a, Vec3 b)
{
	return sqrt(distance2(a, b));
}

static __inline__ double cosAngle(Vec3 v, Vec3 w)
{
	double c = dot(v, w) / (length(v) * length (w));
	/* Known to be strictly larger than 1 due to numerical errors! */
	return MIN(1, c);
}

static __inline__ double angle(Vec3 v, Vec3 w)
{
	return acos(cosAngle(v, w));
}

static __inline__ double dihedral(Vec3 v1, Vec3 v2, Vec3 v3)
{
	Vec3 v1xv2 = cross(v1, v2);
	Vec3 v2xv3 = cross(v2, v3);
	return atan2(length(v2) * dot(v1, v2xv3), dot(v1xv2, v2xv3));
}
/* Returns fills the argument pointers with sin(phi) and cos(phi) where phi 
 * is the dihedral angle. */
static __inline__ void sinCosDihedral(Vec3 v1, Vec3 v2, Vec3 v3, 
					double *sinPhi, double *cosPhi)
{
	/*
	 *    sin(phi) = sin(atan(tan(phi))
	 *             = tan(phi) / sqrt(tan(phi)^2 + 1)
	 * and
	 *    cos(phi) = cos(atan(tan(phi))
	 *             =     1    / sqrt(tan(phi)^2 + 1)
	 */

	Vec3 v1xv2 = cross(v1, v2);
	Vec3 v2xv3 = cross(v2, v3);
	double y = length(v2) * dot(v1, v2xv3);
	double x = dot(v1xv2, v2xv3);
	double tanPhi = y/x;
	double denom = sqrt(SQUARE(tanPhi) + 1);

	/* Because tanPhi is degenerate over additions of pi, we need to 
	 * make this distinction to get back to the full 2pi range of phi 
	 * and the correct sin or cos values change by a sign. */
	if (x > 0) {
		*sinPhi = tanPhi / denom;
		*cosPhi =    1   / denom;
	} else {
		*sinPhi = -tanPhi / denom;
		*cosPhi =    -1   / denom;
	}
}

/* Returns the vector clamped to periodic boundary conditions.
 * This is slower, but works every time. Use the faster functions below for 
 * specific cases.
 * PostConditions:
 *     -period/2.0 <= res.x   &&   res.x < period/2.0
 *     -period/2.0 <= res.y   &&   res.y < period/2.0
 *     -period/2.0 <= res.z   &&   res.z < period/2.0
 */
static __inline__ Vec3 periodic(double period, Vec3 v)
{
	if (LIKELY(length2(v) < SQUARE(period)/4.0))
		return v;

	Vec3 res;
	double hp = period / 2.0; /* Half Period */
	res.x = hp + floor((v.x - hp) / period) * period;
	res.y = hp + floor((v.y - hp) / period) * period;
	res.z = hp + floor((v.z - hp) / period) * period;

	return res;
}

/* Helper for function below */
static __inline__ double _closePeriodic(double period, double val)
{
	if (UNLIKELY(2*val < -period)) {
		do val += period; while (UNLIKELY(2*val < -period));
		return val;
	}
	if (UNLIKELY(2*val >= period)) {
		do val -= period; while (UNLIKELY(2*val >= period));
		return val;
	}
	return val;
}
/* Returns the vector clamped to periodic boundary conditions.
 * Only use this if the vector is "only a couple of times" outside of the 
 * range. If it is far out, use periodic(). If it is at most 1.5 periods 
 * out, use fastPeriodic().
 * PostConditions:
 *     -period/2.0 <= res.x   &&   res.x < period/2.0
 *     -period/2.0 <= res.y   &&   res.y < period/2.0
 *     -period/2.0 <= res.z   &&   res.z < period/2.0
 */
static __inline__ Vec3 closePeriodic(double period, Vec3 v)
{
	// This saves a couple of percents of time...
	if (LIKELY(length2(v) < SQUARE(period)/4.0))
		return v;

	Vec3 res;
	res.x = _closePeriodic(period, v.x);
	res.y = _closePeriodic(period, v.y);
	res.z = _closePeriodic(period, v.z);
	
	assert(-period/2.0 <= res.x  &&  res.x < period/2.0);
	assert(-period/2.0 <= res.y  &&  res.y < period/2.0);
	assert(-period/2.0 <= res.z  &&  res.z < period/2.0);

	return res;
}

/* Helper for function below */
static __inline__ double _fastPeriodic(double period, double val)
{
	//attempt to minimise branches, but still slower:
	//return val + period * ((2*val >= period) - (2*val < -period));

	if (UNLIKELY(2*val < -period))
		return val + period;

	if (UNLIKELY(2*val >= period))
		return val - period;

	return val;
}
/* Returns the vector clamped to periodic boundary conditions and is faster 
 * than the general period() function above.
 * However, we require that (for each dimension) the length of the vector 
 * is at most 1.5*'period', i.e.:
 * PreConditions:
 *     -1.5 * period <= v.x   &&   v.x < 1.5 * period
 *     -1.5 * period <= v.y   &&   v.y < 1.5 * period
 *     -1.5 * period <= v.z   &&   v.z < 1.5 * period
 * PostConditions:
 *     -period/2.0 <= res.x   &&   res.x < period/2.0
 *     -period/2.0 <= res.y   &&   res.y < period/2.0
 *     -period/2.0 <= res.z   &&   res.z < period/2.0
 */
static __inline__ Vec3 fastPeriodic(double period, Vec3 v)
{
	assert(-1.5 * period <= v.x  &&  v.x < 1.5 * period);
	assert(-1.5 * period <= v.y  &&  v.y < 1.5 * period);
	assert(-1.5 * period <= v.z  &&  v.z < 1.5 * period);

	// This saves a couple of percents of time...
	if (LIKELY(length2(v) < SQUARE(period)/4.0))
		return v;

	Vec3 res;
	res.x = _fastPeriodic(period, v.x);
	res.y = _fastPeriodic(period, v.y);
	res.z = _fastPeriodic(period, v.z);

	assert(-period/2.0 <= res.x  &&  res.x < period/2.0);
	assert(-period/2.0 <= res.y  &&  res.y < period/2.0);
	assert(-period/2.0 <= res.z  &&  res.z < period/2.0);

	return res;
}

/* y axis is the vertical axis */
static __inline__ Vec3 fromCilindrical(double r, double phi, double height)
{
	Vec3 res;
	res.x = r * cos(phi);
	res.z = r * sin(phi);
	res.y = height;
	return res;
}

/* These are static *globals* so that inlining randNorm multiple times in a 
 * single function can optimize out the caching of the results! 
 * TODO: verify this */
static double _randNormCache;
static bool _randNormCached = false;
/* Returns a number sampled from a standard normal distribution. */
static __inline__ double randNorm(void)
{
	if (_randNormCached) {
		_randNormCached = false;
		return _randNormCache;
	}

	/* Box-Muller transform */
	/* 0 < u1,u2 <= 1 */
	double u1 = tinymt64_generate_doubleOC(&tinymt);
	double u2 = tinymt64_generate_doubleOC(&tinymt);

	double sqrtLog = sqrt(-2 * log(u1));
	double c = cos(2*M_PI * u2);
	/* TODO: This is faster than a sin()? */
	double s = (u2 < 0.5 ? 1 : -1) * sqrt(1 - c*c);

	_randNormCache = sqrtLog * s;
	_randNormCached = true;
	return sqrtLog * c;
}

/* Returns a vector with components sampled from a standard normal 
 * distribution. */
static __inline__ Vec3 randNormVec(double stdDev)
{
	Vec3 res;
	res.x = randNorm() * stdDev;
	res.y = randNorm() * stdDev;
	res.z = randNorm() * stdDev;
	return res;
}

static __inline__ Vec3 rotate(Vec3 v, Vec3 axis, double theta)
{
	Vec3 u = normalize(axis);
	double ux = u.x;
	double uy = u.y;
	double uz = u.z;
	double c = cos(theta);
	double s = sin(theta);
	Vec3 res;

	/* Rotation matrix fully written out */
	res.x = (c   +    ux*ux*(1 - c)) * v.x
	      + (ux*uy*(1 - c)  -  uz*c) * v.y
	      + (ux*uz*(1 - c)  +  uy*c) * v.z;

	res.y = (ux*uy*(1 - c)  +  uz*s) * v.x
	      + (c   +    uy*uy*(1 - c)) * v.y
	      + (uy*uz*(1 - c)  -  ux*s) * v.z;

	res.z = (ux*uz*(1 - c)  -  uy*s) * v.x
	      + (uy*uz*(1 - c)  +  ux*s) * v.y
	      + (c   +    uz*uz*(1 - c)) * v.z;

	return res;
}

#endif
