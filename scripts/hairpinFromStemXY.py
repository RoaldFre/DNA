#!/usr/bin/python

import sys

if len(sys.argv) < 3:
    print "Usage: "+sys.argv[0]+" <full stem sequence> <num in stem> <num X-Y pairs>"
    sys.exit(1)

seq = sys.argv[1][0 : int(sys.argv[2])]
NXY = int(sys.argv[3])

middle = 'X'*NXY + 'AAAA' + 'Y'*NXY

complementaryBase = {
        'A' : 'T',
        'T' : 'A',
        'C' : 'G',
        'G' : 'C'}

complSeq = ''
for i in range(len(seq) - 1, -1, -1):
    complSeq += complementaryBase.get(seq[i])
    
print(seq + middle + complSeq)

