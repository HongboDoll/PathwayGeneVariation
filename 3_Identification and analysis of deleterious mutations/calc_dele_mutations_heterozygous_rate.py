#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # 
i2 = sys.stdin  # 

dele = {}
for line in i1:
	line = line.strip().split()
	dele[line[0]+'_'+line[1]] = ''


for line in i2:
	if '#' not in line:
		line = line.strip().split()
		if line[0]+'_'+line[1] in dele:
			alt = 0
			ref = 0
			mis = 0
			for n in range(9, len(line)):
				if line[n] == '0/0':
					ref += 1
				elif line[n] == '1/1':
					alt += 1
				elif line[n] == './.':
					mis += 1
			print('%.2f' % (alt/float(alt+ref+mis)))
			
			
			
			
			
			
			
