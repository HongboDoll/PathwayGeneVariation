#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # sweep_region.v6.txt
i2 = open(sys.argv[2])  # DM_v6.1_starch_pathway_gene.xls

sweep = {}
for line in i1:
	line = line.strip().split()
	sweep[line[1]+"~"+line[2]+'~'+line[3]] = ''

for line in i2:
	if '#Chr' not in line:
		line = line.strip().split()
		n = 0
		for k in sweep:
			Chr = k.split('~')[0]
			start = int(k.split('~')[1])
			end = int(k.split('~')[2])
			if line[0] == Chr:
				if int(line[1]) <= start and int(line[2]) >= start and int(line[2]) <= end:
					n += 1
				elif int(line[1]) >= start and int(line[1]) <= end and int(line[2]) >= end:
					n += 1
				elif int(line[1]) >= start and int(line[2]) <= end:
					n += 1
		if n:
			print(line[3], 'domestication sweeps', sep='\t')
		else:
			print(line[3], 'normal', sep='\t')
			
			
			

