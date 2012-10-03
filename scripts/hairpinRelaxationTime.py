#!/usr/bin/env python

import sys

try:
    S = int(sys.argv[1])
except:
    print("Returns the relaxation time (in nanoseconds) for a hairpin with "
            "<S> monomers in the Stem, and 4 monomers in the loop at a "
            "temperature of 90C.")
    print "Usage: "+sys.argv[0]+" <S>"
    exit(1)

N = 2*S + 4 # total number of monomers in the strand
tau1 = 1.07e-7 * (N / 84.0) ** 2.088
tau2 = 3.7e-12 * (N ** 2.29)

print (5*tau2) * 1e9
