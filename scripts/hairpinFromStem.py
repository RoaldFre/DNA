#!/usr/bin/python

import sys

if len(sys.argv) < 2:
    print "Usage: "+sys.argv[0]+" <full stem sequence> [num in stem]"
    sys.exit(1)

seq = sys.argv[1]
if len(sys.argv) > 2:
    seq = seq[0 : int(sys.argv[2])]

middle = 'AAAA'

complementaryBase = {
        'A' : 'T',
        'T' : 'A',
        'C' : 'G',
        'G' : 'C'}

complSeq = ''
for i in range(len(seq) - 1, -1, -1):
    complSeq += complementaryBase.get(seq[i])
    
print(seq + middle + complSeq)

