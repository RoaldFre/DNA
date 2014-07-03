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

#tau1C =
#
#   2.7908e-12   4.9614e-12
#
#tau1E =
#
#   2.37458   0.34041
tau = 2.79e-12 * (N ** 2.37)

relaxtime = (1*tau) * 1e9
if relaxtime < 10:
    relaxtime = 10

print relaxtime
