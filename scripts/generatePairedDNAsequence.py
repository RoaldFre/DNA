#!/usr/bin/env python

import random, sys

try:
    N = int(sys.argv[1])
except:
    print("Generate a random DNA base sequence of length 2<N> of the form "
            "XYXYXY... where every X has an equal probability of being C or "
            "G, and every Y has an equal probability of being A or T. Hence " 
            "there are exactly <N> bases of the C/G type, and <N> bases of "
            "the A/T type.")
    print "Usage: "+sys.argv[0]+" <N>"
    exit(1)

string = ''

for x in range(N):
    string += random.choice("CG")
    string += random.choice("AT")
    
print string

