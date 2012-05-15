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


/* In header for inlining */

static __inline__ bool isSaneNumber(double x)
{
	return !isnan(x) && !isinf(x);
}

static __inline__ bool isSaneVector(Vec3 v)
{
	return isSaneNumber(v.x) && isSaneNumber(v.y) && isSaneNumber(v.z);
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

/* Helper for function below */
static double _periodic(double period, double val)
{
	if (-period/2.0 > val) {
		do val += period; while (-period/2.0 > val);
		return val;
	}
	if (val >= period/2.0) {
		do val -= period; while (val >= period/2.0);
		return val;
	}
	return val;
}
/* Returns the vector clamped to periodic boundary conditions.
 * PostConditions:
 *     -period/2.0 <= res.x   &&   res.x < period/2.0
 *     -period/2.0 <= res.y   &&   res.y < period/2.0
 *     -period/2.0 <= res.z   &&   res.z < period/2.0
 */
static __inline__ Vec3 periodic(double period, Vec3 v)
{
	Vec3 res;
	res.x = _periodic(period, v.x);
	res.y = _periodic(period, v.y);
	res.z = _periodic(period, v.z);
	
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