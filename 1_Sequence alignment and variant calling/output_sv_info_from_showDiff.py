#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # Soltu.DM.01G051390.1_301_Colomba.sv
i2 = open(sys.argv[2])  # Soltu.DM.01G051390.1_301_Colomba.filter.1coords

brk = {}
gap = {}
gap_cor = {}
for line in i1:
	line =  line.strip().split()
	if line[1] == 'BRK' and line[2] != line[3]:
		brk[line[2]+'~'+line[3]] = ''
	elif line[1] == 'GAP':
		if int(line[2]) <= int(line[3]):
			gap[line[2]+'~'+line[3]] = ''
		else:
			gap[line[2]+'~'+line[3]] = ''
		gap_cor[line[2]+'~'+line[3]] = []
if brk:
	for k in brk:
		for line in i2:
			line =  line.strip().split()
			if int(k.split('~')[0]) == 1 and int(k.split('~')[1])+1 == int(line[0]):
				print(line[-2], k.split('~')[0], k.split('~')[1], int(k.split('~')[1])-int(k.split('~')[0])+1, line[-1], line[2], int(line[2]), sep='\t')
			elif int(k.split('~')[0]) != 1 and int(k.split('~')[0])-1 == int(line[1]):
				print(line[-2], k.split('~')[0], k.split('~')[1], int(k.split('~')[1])-int(k.split('~')[0])+1, line[-1], int(line[3]), line[3], sep='\t')
		i2.seek(0)


if gap:
	for k in gap:
		for line in i2:
			line =  line.strip().split()
			if int(line[2]) <= int(line[3]):
				if int(k.split('~')[0])-1 == int(line[1]):
					gap_cor[k].append(int(line[3])+1)
				if int(k.split('~')[1])+1 == int(line[0]):
					gap_cor[k].append(int(line[2])-1)
			else:
				if int(k.split('~')[0])-1 == int(line[1]):
					gap_cor[k].append(int(line[3])+1)
				if int(k.split('~')[1])+1 == int(line[0]):
					gap_cor[k].append(int(line[2])-1)
					gap_cor[k] = gap_cor[k][::-1]
		i2.seek(0)
		if int(k.split('~')[0]) >= int(k.split('~')[1]):
			if gap_cor[k][0] < gap_cor[k][1]:
				print(line[-2], int(k.split('~')[1]),k.split('~')[1], abs(gap_cor[k][1]-gap_cor[k][0]+1-(int(k.split('~')[1])-int(k.split('~')[1]))), line[-1], gap_cor[k][0], gap_cor[k][1], sep='\t')
			else:
				print(line[-2], int(k.split('~')[1]),k.split('~')[1], abs(gap_cor[k][1]-gap_cor[k][1]+1-(int(k.split('~')[1])-int(k.split('~')[1]))), line[-1], gap_cor[k][1], gap_cor[k][1], sep='\t')
		else:
			if gap_cor[k][0] < gap_cor[k][1]:
				print(line[-2], k.split('~')[0], k.split('~')[1], abs(gap_cor[k][1]-gap_cor[k][0]+1-(int(k.split('~')[1])-int(k.split('~')[0])+1)), line[-1], gap_cor[k][0], gap_cor[k][1], sep='\t')
			else:
				print(line[-2], k.split('~')[0], k.split('~')[1], abs(gap_cor[k][1]-gap_cor[k][1]+1-(int(k.split('~')[1])-int(k.split('~')[0])+1)), line[-1], gap_cor[k][1], gap_cor[k][1], sep='\t')

	
				
	
	
	


