#!/bin/bash

cd /media/scratchpad_02/li322/work/01_starch_gene_allelic_variations

####### extract gene haplotype sequences from manually curated results
for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
do
	rm -rf ${i}_gene_seq; mkdir ${i}_gene_seq
	grep -v 'fasta' ${i}_results.xls | grep -v '#' | grep -v '~' | while read n
	do
		name=`echo $n | awk '{print $1}' | awk -F '_' '{print $2}'`
		chr=`echo $n | awk '{print $2}'`
		start=`echo $n | awk '{if($3<$4){print $3}else{print $4}}'`
		end=`echo $n | awk '{if($3<$4){print $4}else{print $3}}'`
		samtools faidx ${i}_Genome.fasta ${chr}:${start}-${end} > ${i}_gene_seq/${name}_${chr}.fa
	done
done

#### extract reference gene and upstream 2k sequences
rm -rf starch_gene_seqs; mkdir starch_gene_seqs
cat DM_v6.1_starch_pathway_gene.xls | grep -v '#Chr' | while read i
do
        chr=`echo $i | awk '{print $1}'`
        strand=`echo $i | awk '{print $5}'`
        if [ $strand = "+" ];then
                s=`echo $i | awk '{print $2-2000}'`
                e=`echo $i | awk '{print $3+2000}'`
        else
                s=`echo $i | awk '{print $2-2000}'`
                e=`echo $i | awk '{print $3+2000}'`
        fi
        g=`echo $i | awk '{print $4}'`
        echo $chr $strand $s $e $g
        samtools faidx $ref ${chr}:${s}-${e} > starch_gene_seqs/${g}_up_down_2k.fa
done

cd starch_gene_seqs
ls *fa | while read i
do
name=`echo $i | awk -F '_' '{print $1}'`
awk '{if($1~/>/){split($1,a,":");print a[1]"_""'"$name"'"}else{print $0}}' $i > ${i}_format.fa
done
rm *_2k.fa
cd -

##### modify GFF3 of these genes
rm DM_v6.1_starch_pathway_gene.gff3
awk '{print $4}' DM_v6.1_starch_pathway_gene.xls |while read i
do
start=`grep $i DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3 | awk '$3=="mRNA"' | awk '{print $4-1}'`
grep $i DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3 | awk 'BEGIN{OFS="\t"}{print $1"_""'"$i"'",$2,$3,$4-"'"$start"'"+2000,$5-"'"$start"'"+2000,$6,$7,$8,$9}'>> DM_v6.1_starch_pathway_gene.gff3
done

export PATH=/media/scratchpad_02/li322/software/mummer-4.0.0rc1/bin/:$PATH

##### pairwise alignment using nucmer
rm -rf delta; mkdir delta
for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
do
        ls ${i}_gene_seq/* | while read n
        do
        name=`echo $n | awk -F '/' '{print $2}' | awk -F '_' '{print $1}'`
        chr_name=`echo $n | awk -F '/' '{print $2}' | awk -F '_' '{print $2}' | sed 's/\.fa//g'`
        pre="${name}_${chr_name}_${i}"
        nucmer --maxmatch starch_gene_seqs/${name}*fa ${n} --delta delta/${pre}.delta
        delta-filter -1 delta/${pre}.delta > delta/${pre}.filter.delta
        done
done

##### variant calling
rm -rf coords/ vcf sv ; mkdir coords/ vcf sv 
ls delta | while read seq
    do
        pre=`echo $seq | awk -F '.' '{print $1"."$2"."$3"."$4}'`
        show-coords -THrcl delta/${seq} > coords/${pre}.filter.1coords
		chr_name=`echo $seq | awk -F '_' '{print $2}'`
		delta2vcf < delta/${seq} | awk 'BEGIN{OFS="\t"}{if($1~/#/&&$1!~/#CHROM/){print $0}else if($1~/#CHROM/){print $1,$2,$3,$4,$5,$6,$7,$8,$9,"'"$chr_name"'"}else{print $1,$2,$3,$4,$5,$6,$7,$8,"GT","1/1"}}' | sed 's/GT,Number=1,Type=Integer,/GT,Number=1,Type=String,/g' > vcf/${pre}.vcf
		show-diff -H delta/${seq} > sv/${pre}.sv
		./output_sv_info_from_showDiff.py sv/${pre}.sv coords/${pre}.filter.1coords > sv/${pre}.sv.coords
		f1=`head -1 delta/${seq} | awk '{print $1}'`
		f2=`head -1 delta/${seq} | awk '{print $2}'`
		awk '{if($1~/>/){print ">seq"}else{print $0}}' $f2 > tmp.fa; rm -rf tmp.fa.fai
		touch sv/${pre}.sv.vcf
		cat sv/${pre}.sv.coords | while read i
		do
			chr=`echo $i | awk '{print $1}'`
			pos1=`echo $i | awk '{print $2}'`
			pos2=`echo $i | awk '{print $3}'`
			pos3=`echo $i | awk '{print $6}'`
			pos4=`echo $i | awk '{print $7}'`
			id="."
			ref=`samtools faidx $f1 ${chr}:${pos1}-${pos2} | grep -v '>' | sed ':a;N;s/\n//g;ta'`
			echo ${pre}
			alt=`samtools faidx tmp.fa seq:${pos3}-${pos4} | grep -v '>' | sed ':a;N;s/\n//g;ta'`
			echo -e "${chr}\t${pos1}\t${id}\t${ref}\t${alt}\t40\tPASS\t.\tGT\t1/1" >> sv/${pre}.sv.vcf
		done
		./filter_vcf_N.py sv/${pre}.sv.vcf > t; mv t sv/${pre}.sv.vcf
		cat <(grep '#' vcf/${pre}.vcf) <(cat vcf/${pre}.vcf sv/${pre}.sv.vcf |grep -v '#'|sort -k1,1 -k2,2n) > t; mv t vcf/${pre}.vcf

done

#### merge variants on four homologues chromosome for each gene (missing chromosomes were ignored)
rm -rf merge_vcf; mkdir merge_vcf
for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
do
    ls vcf/*${i}.vcf | while read j
    do
        bgzip -f ${j}
        tabix -f ${j}.gz
    done
    ls vcf/*${i}.vcf.gz |awk -F '/' '{print $2}' | awk -F '_' '{print $1}'|sort |uniq | while read g
    do
        vcf_count=`ls vcf/${g}*${i}.vcf.gz | wc -l`
        if [ $vcf_count -eq 1 ];then
            zcat vcf/${g}*${i}.vcf.gz > merge_vcf/${g}_merge4HomoChr_${i}.vcf
        else
            bcftools merge -0 --merge all vcf/${g}*${i}.vcf.gz  > merge_vcf/${g}_merge4HomoChr_${i}.vcf
            cat coords/${g}*${i}*coords | ./determine_missing_genotype.py merge_vcf/${g}_merge4HomoChr_${i}.vcf > tmp; mv tmp merge_vcf/${g}_merge4HomoChr_${i}.vcf
        fi
    done
done

###### merge vcf for each gene from the 6 genomes
ls starch_gene_seqs/*fa | awk -F '/' '{print $2}' | awk -F '_' '{print $1}' | while read g
do
   for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
   do
       awk 'BEGIN{OFS="\t"}{if($1~/##/){print $0}else if($1~/#CHROM/){printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t";for(i=10;i<NF;i++){printf "'"$i"'""_"$i"\t"};print "'"$i"'""_"$NF}else{print $0}}' merge_vcf/${g}_merge4HomoChr_${i}.vcf > tmp; mv tmp merge_vcf/${g}_merge4HomoChr_${i}.vcf
       bgzip -f merge_vcf/${g}_merge4HomoChr_${i}.vcf
       tabix -f merge_vcf/${g}_merge4HomoChr_${i}.vcf.gz
   done
done

for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
do
   ls coords/*${i}* | while read n
   do
       awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"'"$i"'""_"$13}' $n > tmp; mv tmp $n
   done
done

ls starch_gene_seqs/*fa | awk -F '/' '{print $2}' | awk -F '_' '{print $1}' | while read g
do
   bcftools merge -0 --merge all merge_vcf/${g}_merge4HomoChr_*.vcf.gz > merge_vcf/${g}_merge4HomoChr_6_genomes.vcf
   cat coords/${g}*1coords | ./determine_missing_genotype.py merge_vcf/${g}_merge4HomoChr_6_genomes.vcf > tmp; mv tmp merge_vcf/${g}_merge4HomoChr_6_genomes.vcf
done

#### convert_gene_vcf2_genome_vcf
rm -rf vcf_genome_coordinates; mkdir vcf_genome_coordinates
ls merge_vcf/*_merge4HomoChr_6_genomes.vcf | while read i
do
    g=`echo $i | awk -F '/' '{print $2}' | awk -F '_' '{print $1}'`
    chr=`grep 'mRNA' DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3 | grep "$g" | awk '{print $1}'`
    pos=`grep 'mRNA' DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3 | grep "$g" | awk '{print $4-2000-1}'`
    awk 'BEGIN{OFS="\t"}{if($1~/#/){print $0}else{printf "'"$chr"'""\t"int("'"$pos"'")+int($2)"\t";for(i=3;i<=NF;i++){printf $i"\t"};print ""}}' $i > vcf_genome_coordinates/${g}_merge4HomoChr_6_genomes.vcf
done

