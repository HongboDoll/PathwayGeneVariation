#!/usr/bin/env python3

import sys

i1 = open(sys.argv[1])  # Soltu.DM.08G006240.1_merge4HomoChr_6_genomes.vcf.snpeff
o1 = open(sys.argv[2], 'w') # up
o2 = open(sys.argv[3], 'w') # down
o3 = open(sys.argv[4], 'w') # intron
o4 = open(sys.argv[5], 'w') # splice
o5 = open(sys.argv[6], 'w') # syn
o6 = open(sys.argv[7], 'w') # nonsyn
o7 = open(sys.argv[8], 'w') # all
gene = str(sys.argv[9]) # gene name

ubi=0
utri=0
utetra=0
ugreat=0
dbi=0
dtri=0
dtetra=0
dgreat=0
ibi=0
itri=0
itetra=0
igreat=0
spbi=0
sptri=0
sptetra=0
spgreat=0
sybi=0
sytri=0
sytetra=0
sygreat=0
nbi=0
ntri=0
ntetra=0
ngreat=0
abi=0
atri=0
atetra=0
agreat=0
for line in i1:
	if '#' not in line:
		line = line.strip().split()
		if line[4].count(',') == 0:
			abi += 1
		elif line[4].count(',') == 1:
			atri += 1
		elif line[4].count(',') == 2:
			atetra += 1
		else:
			agreat += 1
		
		if 'upstream' in line[7]:
			if line[4].count(',') == 0:
				ubi += 1
			elif line[4].count(',') == 1:
				utri += 1
			elif line[4].count(',') == 2:
				utetra += 1
			else:
				ugreat += 1
		if 'downstream' in line[7]:
			if line[4].count(',') == 0:
				dbi += 1
			elif line[4].count(',') == 1:
				dtri += 1
			elif line[4].count(',') == 2:
				dtetra += 1
			else:
				dgreat += 1
		if 'intron' in line[7]:
			if line[4].count(',') == 0:
				ibi += 1
			elif line[4].count(',') == 1:
				itri += 1
			elif line[4].count(',') == 2:
				itetra += 1
			else:
				igreat += 1
		if 'splice' in line[7]:
			if line[4].count(',') == 0:
				spbi += 1
			elif line[4].count(',') == 1:
				sptri += 1
			elif line[4].count(',') == 2:
				sptetra += 1
			else:
				spgreat += 1
		if 'synonymous' in line[7]:
			if line[4].count(',') == 0:
				sybi += 1
			elif line[4].count(',') == 1:
				sytri += 1
			elif line[4].count(',') == 2:
				sytetra += 1
			else:
				sygreat += 1
		if 'MODERATER' in line[7] or 'HIGH' in line[7]:
			if line[4].count(',') == 0:
				nbi += 1
			elif line[4].count(',') == 1:
				ntri += 1
			elif line[4].count(',') == 2:
				ntetra += 1
			else:
				ngreat += 1

o1.write('%s\t%s\t%s\t%s\t%s\n' % (gene, ubi, utri, utetra, ugreat))
o2.write('%s\t%s\t%s\t%s\t%s\n' % (gene, dbi, dtri, dtetra, dgreat))
o3.write('%s\t%s\t%s\t%s\t%s\n' % (gene, ibi, itri, itetra, igreat))
o4.write('%s\t%s\t%s\t%s\t%s\n' % (gene, spbi, sptri, sptetra, spgreat))
o5.write('%s\t%s\t%s\t%s\t%s\n' % (gene, sybi, sytri, sytetra, sygreat))
o6.write('%s\t%s\t%s\t%s\t%s\n' % (gene, nbi, ntri, ntetra, ngreat))
o7.write('%s\t%s\t%s\t%s\t%s\n' % (gene, abi, atri, atetra, agreat))

