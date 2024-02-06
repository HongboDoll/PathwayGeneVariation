#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # SGA_gene_6_genomes_number_variant_each_haplotype_POR06V12_t.xls
i2 = open(sys.argv[2])  # SGA_gene_ID.xls

coores = {}
for line in i2:
	line = line.strip().split()
	coores[line[0]] = line[1]

n = 1
for line in i1:
	if n == 1:
		n += 1
		for k in line.strip().split():
			print(coores[k]+"_"+k, end='\t')
		print()
	else:
		print(line.strip())

