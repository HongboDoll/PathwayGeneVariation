#!/bin/bash

cd /media/bulk_01/users/li322/work/01_atlas_k-mer_starch_SGA_genes

######### compute population genomics indexes using the popgenome R package

vcf=starch_SGA_gene_gatk.filter.bi-allelic.vcf.gz   # argv[1]
window=500   # argv[5]
gff=starch_SGA_gene.gff3   # argv[4]

### To read in SNP data from a subset of individuals the parameter samplenames requires an character vector including the individual names. To extract the individual names from the VCF file do the following:
### vcf_handle  <- .Call("VCF_open",filename)
### ind         <- .Call("VCF_getSampleNames",vcf_handle)
### samplenames <- ind[1:10]
### pop1: starch variety
### pop2: other variety
rm -rf starch_gene_pop_ana_results; mkdir starch_gene_pop_ana_results
grep -v '#Chr' DM_v6.1_starch_pathway_gene.xls | awk '{print $1"_"$4}'| while read i
do
	end=`out_len.py starch_SGA_gene_seqs.fa | grep $i | awk '{print $2}'`  # argv[3]
	ref=$i   # argv[2]
	nuc_div=starch_gene_pop_ana_results/${i}_nucleotide_diversity.xls   # argv[6]
	fst=starch_gene_pop_ana_results/${i}_Fst.xls   # argv[7]
	tajimasd=starch_gene_pop_ana_results/${i}_TajimasD.xls   # argv[8]
	./popgenome.R $vcf $ref $end $gff $window $nuc_div $fst $tajimasd
done

### pop1: high SGA > 100
### pop2: low SGA < 50
rm -rf SGA_gene_pop_ana_results; mkdir SGA_gene_pop_ana_results
grep -v '#Chr' DM_v6.1_SGA_pathway_gene.xls | awk '{print $1"_"$4}'| while read i
do
	end=`out_len.py starch_SGA_gene_seqs.fa | grep $i | awk '{print $2}'`  # argv[3]
	ref=$i   # argv[2]
	nuc_div=SGA_gene_pop_ana_results/${i}_nucleotide_diversity.xls   # argv[6]
	fst=SGA_gene_pop_ana_results/${i}_Fst.xls   # argv[7]
	tajimasd=SGA_gene_pop_ana_results/${i}_TajimasD.xls   # argv[8]
	./popgenome.R2 $vcf $ref $end $gff $window $nuc_div $fst $tajimasd
done

rm starch_gene_pi_tajimasD_Fst_summary.xls
echo -e "gene\tPi\tTajimasD_pop1\tTajimasD_pop2\tFst" >> starch_gene_pi_tajimasD_Fst_summary.xls
grep -v '#Chr' DM_v6.1_starch_pathway_gene.xls | awk '{print $1"_"$4}' | while read i
do
    pi=`cat starch_gene_pop_ana_results/${i}_nucleotide_diversity.xls | grep -v 'pop' | grep -v 'NA' | awk '{i+=$2}END{print i/NR}'`
    taj=`cat starch_gene_pop_ana_results/${i}_TajimasD.xls | grep -v 'pop' | grep -v 'NA' | awk '{i+=$2;j+=$3}END{print i/NR"\t"j/NR}'`
    fst=`cat starch_gene_pop_ana_results/${i}_Fst.xls | grep -v 'pop' | grep -v 'NA' | awk '{i+=$2}END{print i/NR}'`
    echo -e "${i}\t${pi}\t${taj}\t${fst}" >> starch_gene_pi_tajimasD_Fst_summary.xls

done

rm SGA_gene_pi_tajimasD_Fst_summary.xls
echo -e "gene\tPi\tTajimasD_pop1\tTajimasD_pop2\tFst" >> SGA_gene_pi_tajimasD_Fst_summary.xls
grep -v '#Chr' DM_v6.1_SGA_pathway_gene.xls | awk '{print $1"_"$4}' | while read i
do
    pi=`cat SGA_gene_pop_ana_results/${i}_nucleotide_diversity.xls | grep -v 'pop' | grep -v 'NA' | awk '{i+=$2}END{print i/NR}'`
    taj=`cat SGA_gene_pop_ana_results/${i}_TajimasD.xls | grep -v 'pop' | grep -v 'NA' | awk '{i+=$2;j+=$3}END{print i/NR"\t"j/NR}'`
    fst=`cat SGA_gene_pop_ana_results/${i}_Fst.xls | grep -v 'pop' | grep -v 'NA' | awk '{i+=$2}END{print i/NR}'`
    echo -e "${i}\t${pi}\t${taj}\t${fst}" >> SGA_gene_pi_tajimasD_Fst_summary.xls

done

cat <(grep 'gene' starch_gene_pi_tajimasD_Fst_summary.xls) <(grep -v 'gene' starch_gene_pi_tajimasD_Fst_summary.xls | awk -F '_' '{print $2}' | sort -k1,1 | sed 's/-nan/NA/g') > t; mv t starch_gene_pi_tajimasD_Fst_summary.xls
cat <(grep 'gene' SGA_gene_pi_tajimasD_Fst_summary.xls) <(grep -v 'gene' SGA_gene_pi_tajimasD_Fst_summary.xls | awk -F '_' '{print $2}' | sort -k1,1 | sed 's/-nan/NA/g') > t; mv t SGA_gene_pi_tajimasD_Fst_summary.xls
