#!/usr/bin/env python3

import sys

i1 = sys.stdin

gt = {}
corres = {}
for line in i1:
	if '#CHROM' in line:
		line = line.strip().split()
		for k in range(9, len(line)):
			corres[k] = line[k]
	else:
		line = line.strip().split()
		for j in range(9, len(line)):
			k = corres[j]
			if k not in gt:
				gt[k] = ''
				if line[j] == '0/0':
					gt[k] += '0'
				elif line[j] == './.':
					gt[k] += '.'
				else:
					gt[k] += '1'
			else:
				if line[j] == '0/0':
					gt[k] += '0'
				elif line[j] == './.':
					gt[k] += '.'
				else:
					gt[k] += '1'

for k in gt:
	print(k, gt[k])

