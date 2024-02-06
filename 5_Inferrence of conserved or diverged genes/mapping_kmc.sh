#!/bin/bash

threads=10
ref=starch_SGA_gene_seqs.fa

cd /media/bulk_01/users/li322/work/01_atlas_k-mer_starch_SGA_genes

###### mapping resequencing data to a synthetic reference genome that only incorporate starch and SGA genes

bwa index $ref

rm -rf starch_SGA_gene_mapping; mkdir starch_SGA_gene_mapping
cat part1_acc.xls part2_acc.xls | while read i
do
	echo """#!/bin/bash
	cd /media/bulk_01/users/li322/work/01_atlas_k-mer_starch_SGA_genes
	rm -rf starch_SGA_gene_mapping/${i}; mkdir -p starch_SGA_gene_mapping/${i}
	bwa mem -R "@RG\\tID:${i}\\tSM:${i}\\tPL:Illumina" -t $threads $ref <(zcat reseq_data/${i}/*_R1_*gz) <(zcat reseq_data/${i}/*_R2_*gz) | samtools view -@ $threads -F 2048 -F 12 -F 256 -F 8 - -bS > starch_SGA_gene_mapping/${i}/${i}_starch_SGA_gene_mapping.bam
	samtools sort -@ $threads -m 4G -o starch_SGA_gene_mapping/${i}/${i}_starch_SGA_gene_mapping.sort.bam starch_SGA_gene_mapping/${i}/${i}_starch_SGA_gene_mapping.bam
	gatk --java-options "-Xmx40G" MarkDuplicates -I starch_SGA_gene_mapping/${i}/${i}_starch_SGA_gene_mapping.sort.bam -O starch_SGA_gene_mapping/${i}/${i}_starch_SGA_gene_mapping.sort.markdup.bam -M starch_SGA_gene_mapping/${i}/${i}_starch_SGA_gene_mapping.markdup.metrics.txt
    samtools index starch_SGA_gene_mapping/${i}/${i}_starch_SGA_gene_mapping.sort.markdup.bam && rm starch_SGA_gene_mapping/${i}/${i}_starch_SGA_gene_mapping.bam starch_SGA_gene_mapping/${i}/${i}_starch_SGA_gene_mapping.sort.bam
	""" > mapping_${i}.sh && chmod 755 mapping_${i}.sh
done

