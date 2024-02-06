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

