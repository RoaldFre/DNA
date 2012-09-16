#!/usr/bin/env python

import random, sys

try:
    N = int(sys.argv[1])
except:
    print("Generate a random DNA base sequence of length <N> where every "
            "base in the sequence has a uniform probability of being A, C, T "
            "or G.")
    print "Usage: "+sys.argv[0]+" <N>"
    exit(1)

alphabet = 'ACTG'
string = ''

for x in range(N):
    string += random.choice(alphabet)
    
print string
