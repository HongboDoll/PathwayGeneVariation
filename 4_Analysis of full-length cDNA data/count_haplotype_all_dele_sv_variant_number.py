#!/usr/bin/env python3

import sys

i1 = sys.stdin # ${i}_merge_vcf/${gene}_merge4HomoChr_${nn}.vcf.snpeff
i2 = str(sys.argv[1])  # $gene
i3 = str(sys.argv[2])  # ${nn}_${hap}
i4 = open(sys.argv[3]) # starch_deleterious_prediction/Soltu.DM.01G008290.1_merge4HomoChr_Avenger_SIFTannotations.xls
i5 = open(sys.argv[4]) # starch_deleterious_prediction/Soltu.DM.01G008290.1_merge4HomoChr_Avenger_SIFTpredictions.vcf
length = int(sys.argv[5])
typ = str(sys.argv[6]) # all,dele,sv

corres = {}
al = 0
dele = 0
sv = 0
for line in i1:
	if '#CHROM' in line:
		line = line.strip().split()
		for n in range(9, len(line)):
			corres[line[n]] = n
	elif '#' not in line and i2 in line:
		line = line.strip().split()
		nline = corres[i3]
		if line[nline] != '0/0' and line[nline] != './.':
			al += 1
			if abs(len(line[3]) - len(line[4].split(',')[0])) > 50:
				sv += 1

deleterious = {}
for line in i4:
	if 'CHROM' not in line and i2 in line:
		line = line.strip().split()
		if line[16] == 'DELETERIOUS':
			deleterious[line[0]+'~'+line[1]] = ''

for line in i5:
	if '#' not in line:
		line = line.replace(' ','').strip().split()
		nline = corres[i3]
		if line[nline] != '0/0' and line[nline] != './.':
			if line[0]+'~'+line[1] in deleterious:
				dele +=	1

if typ == 'all':
	print('%.2f' % (float(al)*1000/int(length)))
elif typ == 'dele':
	print(dele)
elif typ == 'sv':
	print(sv)
	
