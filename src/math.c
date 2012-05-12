#include "math.h"
#include <unistd.h>
#include <sys/time.h>

/* This is the parameter set from the check example of tinymt64. */
tinymt64_t tinymt = {
	.mat1 = 0xfa051f40,
	.mat2 = 0xffd0fff4,
	.tmat = 0x58d02ffeffbfffbc,
};

void seedRandomWith(uint64_t seed)
{
	tinymt64_init(&tinymt, seed);
}
void seedRandom(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	seedRandomWith(tv.tv_usec ^ tv.tv_sec 
			^ ((uint64_t)getpid() * 123456789));
}
