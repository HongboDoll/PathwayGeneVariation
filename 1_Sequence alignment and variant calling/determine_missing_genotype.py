#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # Soltu.DM.01G049590.2_merge4HomoChr_Spunta.vcf
i2 = sys.stdin  # cat Soltu.DM.01G049590.2*coords

pos = {}
for line in i2:
	line = line.strip().split()
	chr_name = line[-1].split(':')[0]
	if chr_name not in pos:
		pos[chr_name] = []
		pos[chr_name].append([int(line[0]), int(line[1])])
	else:
		pos[chr_name].append([int(line[0]), int(line[1])])
	
hchr = {}
for line in i1:
	if '##' == line[0:2]:
		print(line.strip())
	elif '#CHROM' in line:
		print(line.strip())
		line = line.strip().split()
		for n in range(9, len(line)):
			hchr[n] = line[n]
	else:
		line = line.strip().split()
		p = int(line[1])
		print('\t'.join(line[0:9]), end='\t')
		for n in range(9, len(line)):
			if n != len(line) - 1:
				if line[n] == "0/0":
					f = 0
					for k in pos[hchr[n]]:
						if not f and p >= k[0] and p <= k[1]:
							f += 1
							print("0/0", end='\t')
					else:
						if not f:
							print("./.", end='\t')
				else:
					print(line[n], end='\t')
			else:
				if line[n] == "0/0":
					f = 0
					for k in pos[hchr[n]]:
						if not f and p >= k[0] and p <= k[1]:
							f += 1
							print("0/0")
					else:
						if not f:
							print("./.")
				else:
					print(line[n])

