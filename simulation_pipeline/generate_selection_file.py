#!/usr/bin/env python

import sys

output = "P	F	0	0.5	1	1	1"

if len(sys.argv) != 3:
	exit("Usage: generate_selection_file.py ploidy selective_coefficient \n Ploidy: d for diploid, h for haploid.")

ploidy = sys.argv[1]
select = float(sys.argv[2])

if ploidy == "d":
	sel = select + 1
	sel2 = select/2.0 + 1
elif ploidy == "h":
	sel = 2.0*select + 1
	sel2 = select + 1
else:
	exit("Wrong ploidy. Use d for diploid or h for haploid.")

output += "\t" + str(1.0/sel)
output += "\t" + str(sel2/sel)
output += "\t" + str(1)

print(output)
