#!/usr/bin/env python3

import sys

i1 = sys.stdin  # Soltu.DM.12G028820.2_merge4HomoChr_6_genomes.vcf
spe = str(sys.argv[1]) # Altus

snp = 0
indel = 0
sv = 0
for line in i1:
	if '#' not in line:
		line = line.strip().split()
		if ',' not in line[4]:
			if len(line[3]) == len(line[4]) == 1:
				snp += 1
			elif abs(len(line[3]) - len(line[4])) <= 50:
				indel += 1
			else:
				sv += 1
		else:
			for k in line[4].split(','):
				if len(line[3]) == len(k) == 1:
					snp += 1
				elif abs(len(line[3]) - len(k)) <= 50:
					indel += 1
				else:
					sv += 1
print(spe, snp, indel, sv, sep='\t')
			
			
			
			
			
