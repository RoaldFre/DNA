#include "math.h"
#include "system.h"
#include <unistd.h>

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
uint64_t seedRandom(void)
{
	FILE *stream = fopen("/dev/random", "r");
	uint64_t seed;
	int numRead = fread(&seed, sizeof(seed), 1, stream);
	if (numRead != 1)
		die("Couldn't seed random number generator\n");
	fclose(stream);

	seedRandomWith(seed);
	return seed;
}

bool readSeed(const char *file, uint64_t *seed)
{
	FILE *in = fopen(file, "r");
	if (in == NULL)
		return false;

	long long unsigned theSeed;
	int n = fscanf(in, "%llx", &theSeed);
	*seed = theSeed;
	fclose(in);

	return n == 1;
}

bool writeSeed(const char *file, uint64_t seed)
{
	FILE *out = fopen(file, "w");
	if (out == NULL) {
		return false;
	}

	fprintf(out, "%llx\n%s\n", (long long unsigned) seed, VERSION);
	fclose(out);

	return true;
}

