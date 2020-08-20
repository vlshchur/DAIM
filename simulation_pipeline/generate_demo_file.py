#!/usr/bin/env python

import sys

if len(sys.argv) != 3:
	    exit("Usage: generate_demo_file.py admix_fraction population_size")

impulse = float(sys.argv[1])
pop = sys.argv[2]

print("pop1	pop2	sex	0	1	2\n0	0	F	" + pop + "\t" + pop + "	1")
print("0	a0	F	" + str(1-impulse) + "	0	0")
print("0	a1	F	" + str(impulse) + "	0	0")
print("1	0	F	0	1	0\n1	1	F	" + pop + "	" + pop + "	" + pop + "\n1	a0	F	1	0	0")
