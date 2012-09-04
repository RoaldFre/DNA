#!/usr/bin/python

import random, sys

try: N = int(sys.argv[1])
except IndexError: N = 100

alphabet = 'ACTG'
string = ''

for x in range(N):
    string += random.choice(alphabet)
    
print(string)
