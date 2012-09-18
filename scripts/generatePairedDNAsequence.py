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

seq = ''

for x in range(N):
    seq += random.choice("CG")
    seq += random.choice("AT")
    
print seq
print


# Check length of longest possible mismatch in a hairpin structure

longestMismatch = 0;
mismatchStart = -1;
for i in range(1, len(seq)):
    n = 0;
    while i+n < len(seq)  and  seq[i+n] == seq[n]:
        n += 1
    if (n > longestMismatch):
        longestMismatch = n
        mismatchStart = i

if longestMismatch == 0:
    print "No mismatches possible in a hairpin!"
else:
    print "Longest mismatch in a hairpin: ",longestMismatch," bases, starting at base ",mismatchStart,"."
    print seq
    print " " * (mismatchStart - 1), seq
    print " " * (mismatchStart - 1), "~" * longestMismatch

