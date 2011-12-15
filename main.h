#ifndef _MAIN_H_
#define _MAIN_H_

/* Use these when printing data, as that's how the user inputs them */
#define LENGTH_FACTOR	 		1e-10 /* lengths are in angstrom */
#define TIME_FACTOR	 		1e-15  /* time is in femtoseconds */

void die(const char *fmt, ...);

#endif
