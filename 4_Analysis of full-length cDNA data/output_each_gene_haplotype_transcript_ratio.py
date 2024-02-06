#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # starch_avenger_gene_merged.fa_info.xls

gener = []
for line in i1:
	line = line.strip().split('_')
	gener.append(line[0])
i1.seek(0)

gene = sorted(list(set(gener)))

gene_ratio = {}
for line in i1:
	line = line.strip().split()
	g = line[0].split('_')[0]
	if g not in gene_ratio:
		gene_ratio[g] = [0, 0]
		if int(line[1]) == 0:
			gene_ratio[g][0] += 1
		elif int(line[1]) == 2:
			gene_ratio[g][1] += 1
	else:
		if int(line[1]) == 0:
			gene_ratio[g][0] += 1
		elif int(line[1]) == 2:
			gene_ratio[g][1] += 1

print('gene\tnon-transcrib transcrib')	
for k in gene:
	if k in gene_ratio:
		print(k, '\t', gene_ratio[k][0], ' ', gene_ratio[k][1], sep='')
	
	
