#!/bin/bash

ref=starch_SGA_gene_seqs.fa
threads=1

cd /media/bulk_01/users/li322/work/01_atlas_k-mer_starch_SGA_genes

####### variant calling, genotyping, and filtering using GATK

gatk CreateSequenceDictionary -R $ref -O starch_SGA_gene_seqs.dict

rm -rf starch_SGA_gatk; mkdir starch_SGA_gatk
ls starch_SGA_gene_mapping | while read i
do

	echo """#!/bin/bash
	cd /media/bulk_01/users/li322/work/01_atlas_k-mer_starch_SGA_genes
	gatk --java-options \"-Xmx6G -XX:ParallelGCThreads=1\" HaplotypeCaller -R $ref --sample-ploidy 4 --emit-ref-confidence GVCF -I starch_SGA_gene_mapping/${i}/${i}_starch_SGA_gene_mapping.sort.markdup.bam -O starch_SGA_gatk/${i}.sample.g.vcf --native-pair-hmm-threads $threads # --sample-ploidy assumes samples are polyploidy
	rm starch_SGA_gatk/${i}.sample.g.vcf.gz
	bgzip -f starch_SGA_gatk/${i}.sample.g.vcf
	tabix -f -p vcf starch_SGA_gatk/${i}.sample.g.vcf.gz
""" > ${i}_gatk.sh
done
date
chmod 755 *sh

rm sample.gvcf.xls
ls $PWD/starch_SGA_gatk/*.sample.g.vcf.gz | while read i
do
	name=`echo $i | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}'`
	echo -e "$name\t$i" >> sample.gvcf.xls 
done

rm vcf.list
grep '>' starch_SGA_gene_seqs.fa | sed 's/>//g' | while read i
do
	gatk --java-options "-Xmx200G -XX:ParallelGCThreads=70" GenomicsDBImport -R $ref -L $i --sample-name-map sample.gvcf.xls --genomicsdb-workspace-path starch_SGA_gatk/${i}.gatk.db

	gatk --java-options "-Xmx200G -XX:ParallelGCThreads=70" GenotypeGVCFs -R $ref --allow-old-rms-mapping-quality-annotation-data --sample-ploidy 4 -V gendb://starch_SGA_gatk/${i}.gatk.db -O starch_SGA_gatk/${i}.gatk.vcf # --sample-ploidy
	echo -e "$PWD/starch_SGA_gatk/${i}.gatk.vcf" >> vcf.list
done

gatk MergeVcfs -I vcf.list -O starch_SGA_gene_gatk.vcf

#
gatk SelectVariants -select-type SNP -V starch_SGA_gene_gatk.vcf -O all.snp.vcf.gz

gatk --java-options "-Xmx200G -XX:ParallelGCThreads=70" VariantFiltration -V all.snp.vcf.gz \
--filter-expression "QD < 2.0" --filter-name "LowQD" \
--filter-expression "MQ < 40.0" --filter-name "MQ40.0" \
--filter-expression "FS > 60.0" --filter-name "FS60.0" \
--filter-expression "SOR > 3.0" --filter-name "SOR3.0" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8.0" \
-G-filter "GQ < 20.0" -G-filter-name "GQ20.0" \
-O all.filter.snp.vcf.gz

gatk SelectVariants -select-type INDEL -V starch_SGA_gene_gatk.vcf -O all.indel.vcf.gz

gatk --java-options "-Xmx200G -XX:ParallelGCThreads=70" VariantFiltration -V all.indel.vcf.gz \
--filter-expression "QD < 2.0" --filter-name "LowQD" \
--filter-expression "MQ < 40.0" --filter-name "MQ40.0" \
--filter-expression "FS > 200.0" --filter-name "FS200" \
--filter-expression "SOR > 10.0" --filter-name "SOR10" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8.0" \
-G-filter "GQ < 20.0" -G-filter-name "GQ20.0" \
-O all.filter.indel.vcf.gz

gatk MergeVcfs \
-I all.filter.snp.vcf.gz \
-I all.filter.indel.vcf.gz \
-O starch_SGA_gene_gatk.filter.vcf

cat <(grep '#' starch_SGA_gene_gatk.filter.vcf) <(grep -v '#' starch_SGA_gene_gatk.filter.vcf | awk '$7=="PASS"') | bgzip -f > starch_SGA_gene_gatk.filter.vcf.gz

cat <(grep '#' starch_SGA_gene_gatk.filter.vcf) <(grep -v '#' starch_SGA_gene_gatk.filter.vcf | awk '$7=="PASS"&&$5!~/,/') | bgzip -f > starch_SGA_gene_gatk.filter.bi-allelic.vcf.gz

tabix -f -p vcf starch_SGA_gene_gatk.filter.vcf.gz
tabix -f -p vcf starch_SGA_gene_gatk.filter.bi-allelic.vcf.gz
