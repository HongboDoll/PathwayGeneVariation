#!/bin/bash

for i in SGA starch
do
	ls ${i}_merge_vcf/*gz | while read n
	do
		bgzip -d $n
	done
done

#### convert_gene_vcf2_genome_vcf
for i in starch SGA
do
rm -rf ${i}_vcf_genome_coordinates; mkdir ${i}_vcf_genome_coordinates
ls ${i}_merge_vcf/*.vcf | grep -v '6_genomes' | while read n
do
	name=`echo $n | awk -F '/' '{print $2}'`
    g=`echo $n | awk -F '/' '{print $2}' | awk -F '_' '{print $1}'`
    chr=`grep 'mRNA' DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3 | grep "$g" | awk '{print $1}'`
    pos=`grep 'mRNA' DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3 | grep "$g" | awk '{print $4-2000-1}'`
    awk 'BEGIN{OFS="\t"}{if($1~/#/){print $0}else{printf "'"$chr"'""\t"int("'"$pos"'")+int($2)"\t";for(i=3;i<=NF;i++){printf $i"\t"};print ""}}' $n > ${i}_vcf_genome_coordinates/${name}
done
done

for i in starch  SGA
do
	rm -rf ${i}_deleterious_prediction; mkdir ${i}_deleterious_prediction
	ls ${i}_vcf_genome_coordinates/*vcf | grep -v '6_genomes' | while read n
	do
		n_row=`cat $n | grep -v '#' | wc -l`
		name=`echo $n | awk -F '/' '{print $2}' | sed 's/\.vcf//g'`
		if [ $n_row -eq 0 ];then
			touch ${name}_SIFTannotations.xls; mv ${name}_SIFTannotations.xls ${i}_deleterious_prediction
			touch ${name}_SIFTpredictions.vcf; mv ${name}_SIFTpredictions.vcf ${i}_deleterious_prediction
		else
		java -jar /vol1/agis/huangsanwen_group/lihongbo/software/SIFT4G_Annotator_v2.3.jar -c -i ${n} -d sift_DB_DMv6.1_CL -r ${i}_deleterious_prediction
		fi
	done
done


#### number of variants in 81 genes in 6 genomes

rm starch_gene_6_genomes_deleterious_number.xls
echo -e "Gene\tAltus\tAtlantic\tAvenger\tColomba\tPOR06V12\tSpunta\tAll" >> starch_gene_6_genomes_deleterious_number.xls
ls starch_gene_seqs/* | awk -F '/' '{print $2}' | awk -F '_' '{print $1}' | while read g
do
    echo -e "$g\t\c\n" >> starch_gene_6_genomes_deleterious_number.xls
    for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
    do
        num_v=0
        num_v=`grep -v 'CHROM' starch_deleterious_prediction/${g}_merge4HomoChr_${i}_SIFTannotations.xls | grep 'DELETERIOUS' | wc -l`
        echo -e "$num_v\t\c\n" >> starch_gene_6_genomes_deleterious_number.xls
    done
	echo -e "" >> starch_gene_6_genomes_deleterious_number.xls
done

###### haplotype VS variant number correlation XY

rm starch_gene_6_genomes_haplotype_number_deleterious_number_XY.xls
for i in `seq 2 7`
do
   paste <(grep -v 'Gene' starch_gene_6_genomes_haplotype_number.xls | awk -v a=$i '{print $a}') <(grep -v 'Gene' starch_gene_6_genomes_deleterious_number.xls | awk -v a=$i '{print $a}')  >> starch_gene_6_genomes_haplotype_number_deleterious_number_XY.xls
done



##### number of variants in 38 genes in 6 genomes

rm SGA_gene_6_genomes_deleterious_number.xls
echo -e "Gene\tAltus\tAtlantic\tAvenger\tColomba\tPOR06V12\tSpunta\tAll" >> SGA_gene_6_genomes_deleterious_number.xls
ls SGA_gene_seqs/* | awk -F '/' '{print $2}' | awk -F '_' '{print $1}' | while read g
do
    echo -e "$g\t\c\n" >> SGA_gene_6_genomes_deleterious_number.xls
    for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
    do
        num_v=0
        num_v=`grep -v 'CHROM' SGA_deleterious_prediction/${g}_merge4HomoChr_${i}_SIFTannotations.xls | grep 'DELETERIOUS' | wc -l`
        echo -e "$num_v\t\c\n" >> SGA_gene_6_genomes_deleterious_number.xls
    done
	echo -e "" >> SGA_gene_6_genomes_deleterious_number.xls
done

###### haplotype VS variant number correlation XY

rm SGA_gene_6_genomes_haplotype_number_deleterious_number_XY.xls
for i in `seq 2 7`
do
   paste <(grep -v 'Gene' SGA_gene_6_genomes_haplotype_number.xls | awk -v a=$i '{print $a}') <(grep -v 'Gene' SGA_gene_6_genomes_deleterious_number.xls | awk -v a=$i '{print $a}')  >> SGA_gene_6_genomes_haplotype_number_deleterious_number_XY.xls
done


###### number of deleterious mutations in each cultivars, calculate heterozygous rates
for s in starch SGA
do
	for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
	do
		cat ${s}_deleterious_prediction/*${i}_SIFTannotations.xls |awk '$17=="DELETERIOUS"'|awk '{print $1,$2}'|sort |uniq > ${s}_${i}_all_deleterious.xls
	done
done

for s in starch SGA
do
	for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
	do
		cat ${s}_vcf_genome_coordinates/*${i}.vcf | ./calc_dele_mutations_heterozygous_rate.py ${s}_${i}_all_deleterious.xls > ${s}_${i}_all_deleterious_hetero_rate.xls
	done
done

for s in starch SGA
do
	echo $s
	cat ${s}_*hetero_rate.xls |awk '$1!=1&&$1!=0'|wc -l
	cat ${s}_*hetero_rate.xls |wc -l
done

