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

/* We default to using reals throughout, but you can choose for single 
 * precision floats by defining USE_SINGLE_PRECISION. */
#ifdef USE_SINGLE_PRECISION
	typedef float real;
#       define REAL_EPSILON (1e-3)
#       define SQRT(x)   sqrtf(x)
#       define SIN(x)    sinf(x)
#       define COS(x)    cosf(x)
#       define ASIN(x)   asinf(x)
#       define ACOS(x)   acosf(x)
#       define LOG(x)    logf(x)
#       define EXP(x)    expf(x)
#       define SINCOS(x,s,c) sincosf(x,s,c)
#else
	typedef double real;
#       define REAL_EPSILON (1e-5)
#       define SQRT(x)   sqrt(x)
#       define SIN(x)    sin(x)
#       define COS(x)    cos(x)
#       define ASIN(x)   asin(x)
#       define ACOS(x)   acos(x)
#       define LOG(x)    log(x)
#       define EXP(x)    exp(x)
#       define SINCOS(x,s,c) sincos(x,s,c)
#endif


typedef struct Vec3
{
#ifdef USE_SINGLE_PRECISION
	float xyz[4] __attribute__((aligned(16)));
	/* Loop over 4 elements, in addition/multiplication/..., so the 
	 * compiler (gcc) can automatically vectorize it into a _single_ 
	 * addps/mulps/... instruction. */
#       define NUM_TO_LOOP_OVER 4
#else
	real xyz[3] __attribute__((aligned(16)));
	/* We can only pack two reals in an SSE register, so we need two 
	 * instructions anyway. When looping over 3, the compiler can issue 
	 * one addpd and one additional addsd. */
#       define NUM_TO_LOOP_OVER 3
#endif
} Vec3;

/* Useful for indexing elements in a Vec3 */
enum {
	X = 0,
	Y = 1,
	Z = 2,
};

static __inline__ Vec3 vec3(real x, real y, real z)
{
	Vec3 res;
	res.xyz[X] = x;
	res.xyz[Y] = y;
	res.xyz[Z] = z;
	return res;
}

static __inline__ void fprintVector(FILE *stream, Vec3 v)
{
	fprintf(stream, "%10f\t%10f\t%10f\t", v.xyz[X], v.xyz[Y], v.xyz[Z]);
}

static __inline__ void printVector(Vec3 v)
{
	fprintVector(stdout, v);
}

static __inline__ void fprintVectorExp(FILE *stream, Vec3 v)
{
	fprintf(stream, "%15e %15e %15e", v.xyz[X], v.xyz[Y], v.xyz[Z]);
}
static __inline__ void printVectorExp(Vec3 v)
{
	fprintVectorExp(stdout, v);
}

/* Warning, you can't 'call' the macros below with functions that have side 
 * effects as arguments!
 * Ie, don't do "ABS(functionWithSideEffectsThatReturnsAReal())"
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

static __inline__ bool isSaneNumber(real x)
{
	return !isnan(x) && !isinf(x);
}

static __inline__ bool isSaneVector(Vec3 v)
{
	return isSaneNumber(v.xyz[X]) && isSaneNumber(v.xyz[Y])
			&& isSaneNumber(v.xyz[Z]);
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

static __inline__ Vec3 add(Vec3 a, Vec3 b)
{
	Vec3 res;

	for (int i = 0; i < NUM_TO_LOOP_OVER; i++)
		res.xyz[i] = a.xyz[i] + b.xyz[i];

	return res;
}

static __inline__ Vec3 sub(Vec3 a, Vec3 b)
{
	Vec3 res;

	for (int i = 0; i < NUM_TO_LOOP_OVER; i++)
		res.xyz[i] = a.xyz[i] - b.xyz[i];

	return res;
}

static __inline__ Vec3 scale(Vec3 v, real lambda)
{
	Vec3 res;

	for (int i = 0; i < NUM_TO_LOOP_OVER; i++)
		res.xyz[i] = v.xyz[i] * lambda;

	return res;
}

static __inline__ real dot(Vec3 v, Vec3 w)
{
	float res = 0;

	/* TODO: is this faster (given that we have to set something to 0) 
	 * than just looping over only three floats without packed 
	 * instructions?
	 * NOTE: does setting the "fourth" component to zero make the 
	 * compiler copy the entire original vector somewhere for safety 
	 * when inlining, or can he do that completely in registers??? */
#ifdef USE_SINGLE_PRECISION
	v.xyz[3] = 0; /* Otherwise we get a wrong result! */
#endif
	for (int i = 0; i < NUM_TO_LOOP_OVER; i++)
		res += v.xyz[i] * w.xyz[i];

	return res;
}

static __inline__ real length2(Vec3 v)
{
	return dot(v, v);
}

static __inline__ real length(Vec3 v)
{
	return sqrt(length2(v));
}

static __inline__ Vec3 cross(Vec3 v, Vec3 w)
{
	Vec3 a,b,x;
	/*
	x.x = v.y * w.z  -  v.z * w.y;
	x.y = v.z * w.x  -  v.x * w.z;
	x.z = v.x * w.y  -  v.y * w.x;
	*/
	a.xyz[X] = v.xyz[Y] * w.xyz[Z];
	a.xyz[Y] = v.xyz[Z] * w.xyz[X];
	a.xyz[Z] = v.xyz[X] * w.xyz[Y];

	b.xyz[X] = v.xyz[Z] * w.xyz[Y];
	b.xyz[Y] = v.xyz[X] * w.xyz[Z];
	b.xyz[Z] = v.xyz[Y] * w.xyz[X];

	x = sub(a, b);

	assert(fabs(dot(v, x) / length(v) / length(x)) < 1e-5);
	assert(fabs(dot(w, x) / length(w) / length(x)) < 1e-5);

	return x;
}

static __inline__ Vec3 normalize(Vec3 v)
{
	real l = length(v);
	assert(l != 0);
	return scale(v, 1/l);
}

static __inline__ real distance2(Vec3 a, Vec3 b)
{
	return(length2(sub(a, b)));
}

static __inline__ real distance(Vec3 a, Vec3 b)
{
	return sqrt(distance2(a, b));
}

static __inline__ real cosAngle(Vec3 v, Vec3 w)
{
	real c = dot(v, w) / (length(v) * length (w));
	/* Known to be strictly larger than 1 due to numerical errors! */
	return MIN(MAX(-1, c), 1);
}

static __inline__ real angle(Vec3 v, Vec3 w)
{
	return ACOS(cosAngle(v, w));
}

static __inline__ real dihedral(Vec3 v1, Vec3 v2, Vec3 v3)
{
	Vec3 v1xv2 = cross(v1, v2);
	Vec3 v2xv3 = cross(v2, v3);
	return atan2(length(v2) * dot(v1, v2xv3), dot(v1xv2, v2xv3));
}
/* Returns fills the argument pointers with sin(phi) and cos(phi) where phi 
 * is the dihedral angle. */
static __inline__ void sinCosDihedral(Vec3 v1, Vec3 v2, Vec3 v3, 
					real *sinPhi, real *cosPhi)
{
	real theSin, theCos;
	Vec3 v1xv2 = cross(v1, v2);
	Vec3 v2xv3 = cross(v2, v3);
	real lv1xv2llv2xv3l = sqrt(length2(v1xv2) * length2(v2xv3));
	theCos = dot(v1xv2, v2xv3) / lv1xv2llv2xv3l;
	if (UNLIKELY(SQUARE(theCos) >= 1)) {
		theCos = 1;
		theSin = 0;
	} else {
		if (dot(v1xv2, v3) > 0)
			theSin = sqrt(1 - SQUARE(theCos));
		else
			theSin = -sqrt(1 - SQUARE(theCos));
	}
	*sinPhi = theSin;
	*cosPhi = theCos;
}

/* Returns the vector clamped to periodic boundary conditions.
 * This is slower, but works every time. Use the faster functions below for 
 * specific cases.
 * PostConditions:
 *     -period/2.0 <= res.xyz[X]   &&   res.xyz[X] < period/2.0
 *     -period/2.0 <= res.xyz[Y]   &&   res.xyz[Y] < period/2.0
 *     -period/2.0 <= res.xyz[Z]   &&   res.xyz[Z] < period/2.0
 */
static __inline__ Vec3 periodic(real period, Vec3 v)
{
	if (LIKELY(length2(v) < SQUARE(period)/4.0))
		return v;

	Vec3 res;
	real hp = period / 2.0; /* Half Period */
	res.xyz[X] = hp + floor((v.xyz[X] - hp) / period) * period;
	res.xyz[Y] = hp + floor((v.xyz[Y] - hp) / period) * period;
	res.xyz[Z] = hp + floor((v.xyz[Z] - hp) / period) * period;

	return res;
}

/* Helper for function below */
static __inline__ real _closePeriodic(real period, real val)
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
 *     -period/2.0 <= res.xyz[X]   &&   res.xyz[X] < period/2.0
 *     -period/2.0 <= res.xyz[Y]   &&   res.xyz[Y] < period/2.0
 *     -period/2.0 <= res.xyz[Z]   &&   res.xyz[Z] < period/2.0
 */
static __inline__ Vec3 closePeriodic(real period, Vec3 v)
{
	// This saves a couple of percents of time...
	if (LIKELY(length2(v) < SQUARE(period)/4.0))
		return v;

	Vec3 res;
	res.xyz[X] = _closePeriodic(period, v.xyz[X]);
	res.xyz[Y] = _closePeriodic(period, v.xyz[Y]);
	res.xyz[Z] = _closePeriodic(period, v.xyz[Z]);
	
	assert(-period/2.0 <= res.xyz[X]  &&  res.xyz[X] < period/2.0);
	assert(-period/2.0 <= res.xyz[Y]  &&  res.xyz[Y] < period/2.0);
	assert(-period/2.0 <= res.xyz[Z]  &&  res.xyz[Z] < period/2.0);

	return res;
}

/* Helper for function below */
static __inline__ real _fastPeriodic(real period, real val)
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
 *     -1.5 * period <= v.xyz[X]   &&   v.xyz[X] < 1.5 * period
 *     -1.5 * period <= v.xyz[Y]   &&   v.xyz[Y] < 1.5 * period
 *     -1.5 * period <= v.xyz[Z]   &&   v.xyz[Z] < 1.5 * period
 * PostConditions:
 *     -period/2.0 <= res.xyz[X]   &&   res.xyz[X] < period/2.0
 *     -period/2.0 <= res.xyz[Y]   &&   res.xyz[Y] < period/2.0
 *     -period/2.0 <= res.xyz[Z]   &&   res.xyz[Z] < period/2.0
 */
static __inline__ Vec3 fastPeriodic(real period, Vec3 v)
{
	assert(-1.5 * period <= v.xyz[X]  &&  v.xyz[X] < 1.5 * period);
	assert(-1.5 * period <= v.xyz[Y]  &&  v.xyz[Y] < 1.5 * period);
	assert(-1.5 * period <= v.xyz[Z]  &&  v.xyz[Z] < 1.5 * period);

	// This saves a couple of percents of time...
	if (LIKELY(length2(v) < SQUARE(period)/4.0))
		return v;

	Vec3 res;
	res.xyz[X] = _fastPeriodic(period, v.xyz[X]);
	res.xyz[Y] = _fastPeriodic(period, v.xyz[Y]);
	res.xyz[Z] = _fastPeriodic(period, v.xyz[Z]);

	assert(-period/2.0 <= res.xyz[X]  &&  res.xyz[X] < period/2.0);
	assert(-period/2.0 <= res.xyz[Y]  &&  res.xyz[Y] < period/2.0);
	assert(-period/2.0 <= res.xyz[Z]  &&  res.xyz[Z] < period/2.0);

	return res;
}

/* y axis is the vertical axis */
static __inline__ Vec3 fromCilindrical(real r, real phi, real height)
{
	Vec3 res;
	real s, c;
	SINCOS(phi, &s, &c);
	res.xyz[X] = r * c;
	res.xyz[Z] = r * s;
	res.xyz[Y] = height;
	return res;
}

/* These are static *globals* so that inlining randNorm multiple times in a 
 * single function can optimize out the caching of the results! 
 * TODO: verify this */
static real _randNormCache;
static bool _randNormCached = false;


/* Returns a number sampled from a standard normal distribution. */
static __inline__ real randNorm(void)
{
	if (_randNormCached) {
		_randNormCached = false;
		return _randNormCache;
	}

	/* Box-Muller transform */
	/* 0 < u1,u2 <= 1 */
	real u1 = tinymt64_generate_doubleOC(&tinymt);
	real u2 = tinymt64_generate_doubleOC(&tinymt);

	assert(u1 != 0);
	assert(u2 != 0);

	real sqrtLog = sqrt(-2 * log(u1));
	real c = cos(2*M_PI * u2);
	/* TODO: This is faster than a sin()? */
	real s = (u2 < 0.5 ? 1 : -1) * sqrt(1 - c*c);

	_randNormCache = sqrtLog * s;
	_randNormCached = true;
	return sqrtLog * c;
}

/* Returns a vector with components sampled from a standard normal 
 * distribution. */
static __inline__ Vec3 randNormVec(real stdDev)
{
	Vec3 rnd, res;

	for (int i = 0; i < 3; i++)
		rnd.xyz[i] = randNorm();

#ifdef USE_SINGLE_PRECISION
	rnd.xyz[3] = 0;
#endif

	for (int i = 0; i < NUM_TO_LOOP_OVER; i++)
		res.xyz[i] = rnd.xyz[i] * stdDev;

	return res;
}

static __inline__ Vec3 rotate(Vec3 v, Vec3 axis, real theta)
{
	Vec3 u = normalize(axis);
	real ux = u.xyz[X];
	real uy = u.xyz[Y];
	real uz = u.xyz[Z];
	real s, c;
	SINCOS(theta, &s, &c);
	Vec3 res;

	/* Rotation matrix fully written out */
	res.xyz[X] = (c   +    ux*ux*(1 - c)) * v.xyz[X]
	      + (ux*uy*(1 - c)  -  uz*c) * v.xyz[Y]
	      + (ux*uz*(1 - c)  +  uy*c) * v.xyz[Z];

	res.xyz[Y] = (ux*uy*(1 - c)  +  uz*s) * v.xyz[X]
	      + (c   +    uy*uy*(1 - c)) * v.xyz[Y]
	      + (uy*uz*(1 - c)  -  ux*s) * v.xyz[Z];

	res.xyz[Z] = (ux*uz*(1 - c)  -  uy*s) * v.xyz[X]
	      + (uy*uz*(1 - c)  +  ux*s) * v.xyz[Y]
	      + (c   +    uz*uz*(1 - c)) * v.xyz[Z];

	return res;
}


typedef struct Mat3
{
	Vec3 r1, r2, r3; /* Rows of the matrix */
} Mat3;

static __inline__ Mat3 mat3(real m11, real m12, real m13,
                            real m21, real m22, real m23,
                            real m31, real m32, real m33)
{
	Mat3 m;
	m.r1 = vec3(m11, m12, m13);
	m.r2 = vec3(m21, m22, m23);
	m.r3 = vec3(m31, m32, m33);
	return m;
}

static __inline__ Mat3 matScale(Mat3 m, real lambda)
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
	w.xyz[X] = dot(m.r1, v);
	w.xyz[Y] = dot(m.r2, v);
	w.xyz[Z] = dot(m.r3, v);
	return w;
}




/* Check for equality of reals (up to some small error). */
static __inline__ bool equalsEpsilon(real a, real b)
{
	if (a + b == 0)
		return a == 0;

	return fabs((a-b) / (a+b)) < REAL_EPSILON;
}

/* Check for equality of vectors (up to some small error). */
static __inline__ bool vecEqualsEpsilon(Vec3 a, Vec3 b)
{
	return length(sub(a,b)) / length(add(a,b)) < REAL_EPSILON;
}
static __inline__ void assertVecEqualsEpsilon(Vec3 a, Vec3 b)
{
#ifndef DEBUG
	UNUSED(a);
	UNUSED(b);
	return;
#else
	if (vecEqualsEpsilon(a, b))
		return;

	fprintf(stderr, "Different vectors!\n");
	fprintVectorExp(stderr, a);
	fprintf(stderr, "\n");
	fprintVectorExp(stderr, b);
	fprintf(stderr, "\n");
	assert(false);
#endif
}

#endif
