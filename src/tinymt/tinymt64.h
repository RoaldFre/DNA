#ifndef TINYMT64_H
#define TINYMT64_H
/**
 * @file tinymt64.h
 *
 * @brief Tiny Mersenne Twister only 127 bit internal state
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2011 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <stdint.h>
#include <inttypes.h>

#define TINYMT64_MEXP 127
#define TINYMT64_SH0 12
#define TINYMT64_SH1 11
#define TINYMT64_SH8 8
#define TINYMT64_MASK UINT64_C(0x7fffffffffffffff)
#define TINYMT64_MUL (1.0 / 18446744073709551616.0)

/*
 * tinymt64 internal state vector and parameters
 */
struct TINYMT64_T {
    uint64_t status[2];
    uint32_t mat1;
    uint32_t mat2;
    uint64_t tmat;
};

typedef struct TINYMT64_T tinymt64_t;

void tinymt64_init(tinymt64_t * rnd, uint64_t seed);
void tinymt64_init_by_array(tinymt64_t * rnd, const uint64_t init_key[],
			    int key_length);

#if defined(__GNUC__)
/**
 * This function always returns 127
 * @param rnd not used
 * @return always 127
 */
inline static int tinymt64_get_mexp(
    tinymt64_t * rnd  __attribute__((unused))) {
    return TINYMT64_MEXP;
}
#else
inline static int tinymt64_get_mexp(tinymt64_t * rnd) {
    return TINYMT64_MEXP;
}
#endif

/**
 * This function changes internal state of tinymt64.
 * Users should not call this function directly.
 * @param rnd tinymt internal status
 */
inline static void tinymt64_next_state(tinymt64_t * rnd) {
    uint64_t x;

    rnd->status[0] &= TINYMT64_MASK;
    x = rnd->status[0] ^ rnd->status[1];
    x ^= x << TINYMT64_SH0;
    x ^= x >> 32;
    x ^= x << 32;
    x ^= x << TINYMT64_SH1;
    rnd->status[0] = rnd->status[1];
    rnd->status[1] = x;
    rnd->status[0] ^= -((int64_t)(x & 1)) & rnd->mat1;
    rnd->status[1] ^= -((int64_t)(x & 1)) & (((uint64_t)rnd->mat2) << 32);
}

/**
 * This function outputs 64-bit unsigned integer from internal state.
 * Users should not call this function directly.
 * @param rnd tinymt internal status
 * @return 64-bit unsigned pseudorandom number
 */
inline static uint64_t tinymt64_temper(tinymt64_t * rnd) {
    uint64_t x;
#if defined(LINEARITY_CHECK)
    x = rnd->status[0] ^ rnd->status[1];
#else
    x = rnd->status[0] + rnd->status[1];
#endif
    x ^= rnd->status[0] >> TINYMT64_SH8;
    x ^= -((int64_t)(x & 1)) & rnd->tmat;
    return x;
}

/**
 * This function outputs floating point number from internal state.
 * Users should not call this function directly.
 * @param rnd tinymt internal status
 * @return floating point number r (1.0 <= r < 2.0)
 */
inline static double tinymt64_temper_conv(tinymt64_t * rnd) {
    uint64_t x;
    union {
	uint64_t u;
	double d;
    } conv;
#if defined(LINEARITY_CHECK)
    x = rnd->status[0] ^ rnd->status[1];
#else
    x = rnd->status[0] + rnd->status[1];
#endif
    x ^= rnd->status[0] >> TINYMT64_SH8;
    conv.u = ((x ^ (-((int64_t)(x & 1)) & rnd->tmat)) >> 12)
	| UINT64_C(0x3ff0000000000000);
    return conv.d;
}

/**
 * This function outputs floating point number from internal state.
 * Users should not call this function directly.
 * @param rnd tinymt internal status
 * @return floating point number r (1.0 < r < 2.0)
 */
inline static double tinymt64_temper_conv_open(tinymt64_t * rnd) {
    uint64_t x;
    union {
	uint64_t u;
	double d;
    } conv;
#if defined(LINEARITY_CHECK)
    x = rnd->status[0] ^ rnd->status[1];
#else
    x = rnd->status[0] + rnd->status[1];
#endif
    x ^= rnd->status[0] >> TINYMT64_SH8;
    conv.u = ((x ^ (-((int64_t)(x & 1)) & rnd->tmat)) >> 12)
	| UINT64_C(0x3ff0000000000001);
    return conv.d;
}

/**
 * This function outputs 64-bit unsigned integer from internal state.
 * @param rnd tinymt internal status
 * @return 64-bit unsigned integer r (0 <= r < 2^64)
 */
inline static uint64_t tinymt64_generate_uint64(tinymt64_t * rnd) {
    tinymt64_next_state(rnd);
    return tinymt64_temper(rnd);
}

/**
 * This function outputs floating point number from internal state.
 * This function is implemented using multiplying by 1 / 2^64.
 * @param rnd tinymt internal status
 * @return floating point number r (0.0 <= r < 1.0)
 */
inline static double tinymt64_generate_double(tinymt64_t * rnd) {
    tinymt64_next_state(rnd);
    return tinymt64_temper(rnd) * TINYMT64_MUL;
}

/**
 * This function outputs floating point number from internal state.
 * This function is implemented using union trick.
 * @param rnd tinymt internal status
 * @return floating point number r (0.0 <= r < 1.0)
 */
inline static double tinymt64_generate_double01(tinymt64_t * rnd) {
    tinymt64_next_state(rnd);
    return tinymt64_temper_conv(rnd) - 1.0;
}

/**
 * This function outputs floating point number from internal state.
 * This function is implemented using union trick.
 * @param rnd tinymt internal status
 * @return floating point number r (1.0 <= r < 2.0)
 */
inline static double tinymt64_generate_double12(tinymt64_t * rnd) {
    tinymt64_next_state(rnd);
    return tinymt64_temper_conv(rnd);
}

/**
 * This function outputs floating point number from internal state.
 * This function is implemented using union trick.
 * @param rnd tinymt internal status
 * @return floating point number r (0.0 < r <= 1.0)
 */
inline static double tinymt64_generate_doubleOC(tinymt64_t * rnd) {
    tinymt64_next_state(rnd);
    return 2.0 - tinymt64_temper_conv(rnd);
}

/**
 * This function outputs floating point number from internal state.
 * This function is implemented using union trick.
 * @param rnd tinymt internal status
 * @return floating point number r (0.0 < r < 1.0)
 */
inline static double tinymt64_generate_doubleOO(tinymt64_t * rnd) {
    tinymt64_next_state(rnd);
    return tinymt64_temper_conv_open(rnd) - 1.0;
}

#endif
