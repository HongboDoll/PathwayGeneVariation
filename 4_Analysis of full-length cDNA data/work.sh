#!/bin/bash

cd /media/scratchpad_02/li322/work/04_full_length_transcript
threads=20

#cp /media/scratchpad_03/lies026/PotatoMamy/stringTie/mRNA/atlantic/*fasta .
#cp /media/scratchpad_03/lies026/PotatoMamy/stringTie/mRNA/castleRusset/*fasta .
#cp /media/scratchpad_03/lies026/PotatoMamy/stringTie/mRNA/Iso-Seq_fzj/*_phased.fa .


for g in starch SGA
do
	for i in altus  atlantic  avenger  castleRusset  colomba  spunta
	do
		ls ${g}_gene_seq_6_genomes/${i}/*fa | while read s
		do
			name=`echo $s | awk -F '/' '{print $NF}' | awk -F '_' '{print $1}'`
			awk '{if($1~/>/){split($1,a,">");print ">""'"$name"'""_"a[2]}else{print $0}}' $s > ${s}.rename
		done
	done
done
	
####### map full-length cDNA sequencing data to the merged starch/SGA gene sets using gmap
for i in altus  atlantic  avenger  castleRusset  colomba  spunta
do
	for g in starch SGA
	do
		cat ${g}_gene_seq_6_genomes/${i}/*rename > ${g}_${i}_gene_merged.fa
		gmap_build -t ${threads} -d ${g}_${i}_ref -D ${g}_${i}_gmap_index ${g}_${i}_gene_merged.fa
		gmap -t ${threads} -D ${g}_${i}_gmap_index/ -d ${g}_${i}_ref --cross-species -n 0 -f samse  -z auto ${i}*.fa > ${g}_${i}_gmap.sam
		cat <(awk '$1~/^@/' ${g}_${i}_gmap.sam) <(grep "NH:i:1\b" ${g}_${i}_gmap.sam | awk '$2!=4&&$2!=2048&&$5>=1') | samtools sort -@ ${threads} > ${g}_${i}_gmap.filter.bam ### filter criteria: remove unmapped/secondarymapped, mapq <1, keep uniquely mapped reads only
		stringtie ${g}_${i}_gmap.filter.bam -p ${threads} -f 0.1 -L -o ${g}_${i}_gmap.filter.gtf
	done
done

####### proportions of transcriable haplotypes
for n in starch SGA
do
	./output_each_gene_haplotype_transcript_ratio.py ${n}_altus_gene_merged.fa_info.xls | awk '{print $1}' > ${n}_gene
	for s in altus atlantic avenger colomba castleRusset spunta
	do
		cat <(echo -e "${s}_0 ${s}_1") <(./output_each_gene_haplotype_transcript_ratio.py ${n}_${s}_gene_merged.fa_info.xls | grep -v 'gene' | awk '{print $2,$3}') > ${n}_${s}_tmp
	done
	list=`ls ${n}_*_tmp | sed ':a;N;s/\n/ /g;ta'`
	paste ${n}_gene $list > ${n}_6_genomes_haplotype_transcribe.xls
done

cat starch_6_genomes_haplotype_transcribe.xls |awk '$3==0&&$5==0&&$7==0&&$9==0&&$11==0&&$13==0'
cat starch_6_genomes_haplotype_transcribe.xls | grep -v 'Soltu.DM.03G007720.1' | grep -v 'Soltu.DM.03G007760.1' | grep -v 'Soltu.DM.07G018130.1' | grep -v 'Soltu.DM.07G018140.1' | grep -v 'Soltu.DM.08G006240.1' | grep -v 'Soltu.DM.10G007060.1' > starch_6_genomes_haplotype_transcribe.xls_4_plot

cat SGA_6_genomes_haplotype_transcribe.xls|awk '$3==0&&$5==0&&$7==0&&$9==0&&$11==0&&$13==0'
cat SGA_6_genomes_haplotype_transcribe.xls|grep -v Soltu.DM.10G016360.1 | grep -v Soltu.DM.02G026070.1 > SGA_6_genomes_haplotype_transcribe.xls_4_plot

paste starch_gene_4_plot <(grep -v 'gene' starch_6_genomes_haplotype_transcribe.xls_4_plot |awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}') > starch_6_genomes_haplotype_transcribe.xls_4_plot2

paste SGA_gene_4_plot <(grep -v 'gene' SGA_6_genomes_haplotype_transcribe.xls_4_plot |awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}') > SGA_6_genomes_haplotype_transcribe.xls_4_plot2

####### number of variants on transcriable/untranscriable haplotypes
####### all_v: variant density/1kb, dele: number of deleterious mutations, sv: number of SVs
for i in starch SGA
do
	rm ${i}_non_transcrib_gene_haplotype_all_dele_sv_number.xls ${i}_transcrib_gene_haplotype_all_dele_sv_number.xls
	for n in altus atlantic avenger colomba castleRusset spunta
	do
		if [ $n = "altus" ];then
			nn="Altus"
		elif [ $n = "atlantic" ];then
			nn="Atlantic"
		elif [ $n = "avenger" ];then
			nn="Avenger"
		elif [ $n = "colomba" ];then
			nn="Colomba"
		elif [ $n = "castleRusset" ];then
			nn="POR06V12"
		elif [ $n = "spunta" ];then
			nn="Spunta"
		fi
		awk '$2==0' ${i}_${n}_gene_merged.fa_info.xls | awk '{print $1}' | grep -v 'Soltu.DM.03G007720.1' | grep -v 'Soltu.DM.03G007760.1' | grep -v 'Soltu.DM.07G018130.1' | grep -v 'Soltu.DM.07G018140.1' | grep -v 'Soltu.DM.08G006240.1' | grep -v 'Soltu.DM.10G007060.1' | grep -v 'Soltu.DM.10G016360.1' | grep -v 'Soltu.DM.02G026070.1' | while read g
		do
			gene=`echo $g | awk -F '_' '{print $1}'`
			hap=`echo $g | awk -F ':' '{print $1}' | awk -F '_' '{print $2}'`
			start=`echo $g | awk -F ':' '{print $2}' | awk -F '-' '{print $1}'`
			end=`echo $g | awk -F ':' '{print $2}' | awk -F '-' '{print $2}'`
			length=`expr $end - $start`
			echo $gene $hap
			all_v=`cat ${i}_merge_vcf/${gene}_merge4HomoChr_${nn}.vcf.snpeff | ./count_haplotype_all_dele_sv_variant_number.py $gene ${nn}_${hap} ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTannotations.xls ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTpredictions.vcf $length all` 
			dele_v=`cat ${i}_merge_vcf/${gene}_merge4HomoChr_${nn}.vcf.snpeff | ./count_haplotype_all_dele_sv_variant_number.py $gene ${nn}_${hap} ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTannotations.xls ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTpredictions.vcf $length dele`
			sv_v=`cat ${i}_merge_vcf/${gene}_merge4HomoChr_${nn}.vcf.snpeff | ./count_haplotype_all_dele_sv_variant_number.py $gene ${nn}_${hap} ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTannotations.xls ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTpredictions.vcf $length sv`
			echo -e "${gene}_${hap}\t$all_v\t$dele_v\t$sv_v" >> ${i}_non_transcrib_gene_haplotype_all_dele_sv_number.xls
		done
		
		awk '$2==2' ${i}_${n}_gene_merged.fa_info.xls | awk '{print $1}' | grep -v 'Soltu.DM.03G007720.1' | grep -v 'Soltu.DM.03G007760.1' | grep -v 'Soltu.DM.07G018130.1' | grep -v 'Soltu.DM.07G018140.1' | grep -v 'Soltu.DM.08G006240.1' | grep -v 'Soltu.DM.10G007060.1' | grep -v 'Soltu.DM.10G016360.1' | grep -v 'Soltu.DM.02G026070.1' | while read g
		do
			gene=`echo $g | awk -F '_' '{print $1}'`
			hap=`echo $g| awk -F ':' '{print $1}' | awk -F '_' '{print $2}'`
            start=`echo $g | awk -F ':' '{print $2}' | awk -F '-' '{print $1}'`
            end=`echo $g | awk -F ':' '{print $2}' | awk -F '-' '{print $2}'`
			length=`expr $end - $start`
			all_v=`cat ${i}_merge_vcf/${gene}_merge4HomoChr_${nn}.vcf.snpeff | ./count_haplotype_all_dele_sv_variant_number.py $gene ${nn}_${hap} ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTannotations.xls ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTpredictions.vcf $length all` 
			dele_v=`cat ${i}_merge_vcf/${gene}_merge4HomoChr_${nn}.vcf.snpeff | ./count_haplotype_all_dele_sv_variant_number.py $gene ${nn}_${hap} ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTannotations.xls ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTpredictions.vcf $length dele`
			sv_v=`cat ${i}_merge_vcf/${gene}_merge4HomoChr_${nn}.vcf.snpeff | ./count_haplotype_all_dele_sv_variant_number.py $gene ${nn}_${hap} ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTannotations.xls ${i}_deleterious_prediction/${gene}_merge4HomoChr_${nn}_SIFTpredictions.vcf $length sv`
			echo -e "${gene}_${hap}\t$all_v\t$dele_v\t$sv_v" >> ${i}_transcrib_gene_haplotype_all_dele_sv_number.xls
		done	
	done
done
