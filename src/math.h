#ifndef _VMATH_H_
#define _VMATH_H_

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "tinymt/tinymt64.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

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
static __inline__ Vec3 vec3(double x, double y, double z)
{
	return (Vec3) {x, y, z};
}

static __inline__ void fprintVector(FILE *stream, Vec3 v)
{
	fprintf(stream, "%10lf\t%10lf\t%10lf", v.x, v.y, v.z);
}

static __inline__ void printVector(Vec3 v)
{
	fprintVector(stdout, v);
}

static __inline__ void fprintVectorExp(FILE *stream, Vec3 v)
{
	fprintf(stream, "%15le %15le %15le", v.x, v.y, v.z);
}
static __inline__ void printVectorExp(Vec3 v)
{
	fprintVectorExp(stdout, v);
}

static __inline__ bool fscanVector(FILE *stream, Vec3 *v)
{
	int n = fscanf(stream, "%lf %lf %lf", &v->x, &v->y, &v->z);
	return (n == 3);
}
static __inline__ bool fscanVectorExp(FILE *stream, Vec3 *v)
{
	int n = fscanf(stream, "%le %le %le", &v->x, &v->y, &v->z);
	return (n == 3);
}

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

/* To stop compiler from issuing unused warnings when it is intentional. */
#define UNUSED(x) ((void) x)

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
#define DEBUG_VECTOR_SANITY false

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
static __inline__ bool equalsEpsilon(double a, double b, double eps)
{
	if (a + b == 0)
		return a == 0;

	return fabs((a-b) / (a+b)) < eps;
}

/* Check for equality of vectors (up to some small error). */
static __inline__ bool vecEqualsEpsilon(Vec3 a, Vec3 b, double eps)
{
	return equalsEpsilon(a.x, b.x, eps)
	    && equalsEpsilon(a.y, b.y, eps)
	    && equalsEpsilon(a.z, b.z, eps);
}
static __inline__ void assertVecEqualsEpsilon(Vec3 a, Vec3 b, double eps)
{
#ifndef DEBUG
	UNUSED(a);
	UNUSED(b);
	UNUSED(eps);
	return;
#else
	if (vecEqualsEpsilon(a, b, eps))
		return;

	fprintf(stderr, "Different vectors!\n");
	fprintVectorExp(stderr, a);
	fprintf(stderr, "\n");
	fprintVectorExp(stderr, b);
	fprintf(stderr, "\n");
	assert(false);
#endif
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
/* Fills the argument pointers with sin(phi) and cos(phi) where phi 
 * is the dihedral angle. */
static __inline__ void sinCosDihedral(Vec3 v1, Vec3 v2, Vec3 v3, 
					double *sinPhi, double *cosPhi)
{
	double theSin, theCos;
	Vec3 v1xv2 = cross(v1, v2);
	Vec3 v2xv3 = cross(v2, v3);
	double lv1xv2llv2xv3l = sqrt(length2(v1xv2) * length2(v2xv3));
	theCos = dot(v1xv2, v2xv3) / lv1xv2llv2xv3l;
	if (UNLIKELY(SQUARE(theCos) >= 1)) {
		theCos = 1;
		theSin = 0;
	} else {
		if (dot(v1xv2, v3) > 0)
			theSin = -sqrt(1 - SQUARE(theCos));
		else
			theSin = sqrt(1 - SQUARE(theCos));
	}
	*sinPhi = theSin;
	*cosPhi = theCos;
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

/* Returns a uniform random number x, where 0 <= x < 1. */
static __inline__ double rand01(void)
{
	return tinymt64_generate_double01(&tinymt);
}

/* Returns a uniform random index x, where 0 <= x < numElements. */
static __inline__ int randIndex(int numElements)
{
	return (int) (numElements * rand01());
}

static __inline__ Vec3 randUniformVec(double low, double high)
{
	Vec3 res;
	res.x = low + (high - low) * rand01();
	res.y = low + (high - low) * rand01();
	res.z = low + (high - low) * rand01();
	return res;
}

/* Return a unit vector uniformly sampled over the surface of a sphere */
static __inline__ Vec3 randomDirection(void)
{
	double theta = 2*M_PI * rand01();
	double s,c;
	sincos(theta, &s, &c);
	double z = 2*rand01() - 1;
	double r = sqrt(1 - SQUARE(z));
	return vec3(r*c, r*s, z);
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
	      + (ux*uy*(1 - c)  -  uz*s) * v.y
	      + (ux*uz*(1 - c)  +  uy*s) * v.z;

	res.y = (ux*uy*(1 - c)  +  uz*s) * v.x
	      + (c   +    uy*uy*(1 - c)) * v.y
	      + (uy*uz*(1 - c)  -  ux*s) * v.z;

	res.z = (ux*uz*(1 - c)  -  uy*s) * v.x
	      + (uy*uz*(1 - c)  +  ux*s) * v.y
	      + (c   +    uz*uz*(1 - c)) * v.z;

	return res;
}
static __inline__ Vec3 rotateAround(Vec3 v, Vec3 origin, Vec3 axis, double theta)
{
	return add(origin, rotate(sub(v, origin), axis, theta));
}


typedef struct Mat3
{
	Vec3 r1, r2, r3; /* Rows of the matrix */
} Mat3;

static __inline__ Mat3 mat3(double m11, double m12, double m13,
                            double m21, double m22, double m23,
                            double m31, double m32, double m33)
{
	Mat3 m;
	m.r1 = (Vec3) {m11, m12, m13};
	m.r2 = (Vec3) {m21, m22, m23};
	m.r3 = (Vec3) {m31, m32, m33};
	return m;
}

static __inline__ Mat3 matScale(Mat3 m, double lambda)
{
	Mat3 s;
	s.r1 = scale(m.r1, lambda);
	s.r2 = scale(m.r2, lambda);
	s.r3 = scale(m.r3, lambda);
	return s;
}

static __inline__ Mat3 matAdd(Mat3 a, Mat3 b)
{
	Mat3 c;
	c.r1 = add(a.r1, b.r1);
	c.r2 = add(a.r2, b.r2);
	c.r3 = add(a.r3, b.r3);
	return c;
}

static __inline__ Mat3 matSub(Mat3 a, Mat3 b)
{
	Mat3 c;
	c.r1 = sub(a.r1, b.r1);
	c.r2 = sub(a.r2, b.r2);
	c.r3 = sub(a.r3, b.r3);
	return c;
}

/* Multiply the matrix m with the column vector v. */
static __inline__ Vec3 matApply(Mat3 m, Vec3 v)
{
	Vec3 w;
	w.x = dot(m.r1, v);
	w.y = dot(m.r2, v);
	w.z = dot(m.r3, v);
	return w;
}
#endif
