#!/usr/bin/env python3

import sys

i1 = sys.stdin # starch_spunta_gene_merged.fa_info.xls
i2 = open(sys.argv[1])  # starh_SGA_gene_4_plot
o1 = open(sys.argv[2], 'w')  # hap 1
o2 = open(sys.argv[3], 'w')  # hap 2
o3 = open(sys.argv[4], 'w')  # hap 3
o4 = open(sys.argv[5], 'w')  # hap 4

gene_hap_num = {}
for line in i2:
	gene_hap_num[line.strip()] = 0

gene_hap_trans = {}
for line in i1:
	line = line.strip().split()
	g_name = line[0].split('_')[0]
	if g_name in gene_hap_num:
		gene_hap_num[g_name] += 1
		if g_name not in gene_hap_trans:
			gene_hap_trans[g_name] = [0, 0]
			if int(line[-1]) == 0:
				gene_hap_trans[g_name][0] += 1
			else:
				gene_hap_trans[g_name][1] += 1
		else:
			if int(line[-1]) == 0:
				gene_hap_trans[g_name][0] += 1
			else:
				gene_hap_trans[g_name][1] += 1
			


for g in gene_hap_num:
	if gene_hap_num[g] == 1:
		o1.write('%.2f\n' % (gene_hap_trans[g][1]*100/float(gene_hap_trans[g][0]+gene_hap_trans[g][1])))
	elif gene_hap_num[g] == 2:
		o2.write('%.2f\n' % (gene_hap_trans[g][1]*100/float(gene_hap_trans[g][0]+gene_hap_trans[g][1])))
	elif gene_hap_num[g] == 3:
		o3.write('%.2f\n' % (gene_hap_trans[g][1]*100/float(gene_hap_trans[g][0]+gene_hap_trans[g][1])))
	elif gene_hap_num[g] == 4:
		o4.write('%.2f\n' % (gene_hap_trans[g][1]*100/float(gene_hap_trans[g][0]+gene_hap_trans[g][1])))

