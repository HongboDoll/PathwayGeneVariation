#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # vcf


for line in i1:
	line = line.strip().split()
	n_count1 = line[3].count('N')+line[3].count('n')
	n_count2 = line[4].count('N')+line[4].count('n')
	if n_count1/len(line[3]) < 0.2 and n_count2/len(line[4]) < 0.2:
		print('\t'.join(line))

