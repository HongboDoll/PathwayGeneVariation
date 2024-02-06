#!/bin/bash

###### expression in DM potato reference genome in seven tissues

rm DM_v6.1.starch_gene.tpm.txt
ls starch_gene_seqs/*fa | awk -F '/' '{print $2}' | awk -F '_' '{print $1}' | while read i
do
	grep $i DM_v6.1.tpm.txt >> DM_v6.1.starch_gene.tpm.txt
done
cat <(head -1 DM_v6.1.tpm.txt) DM_v6.1.starch_gene.tpm.txt > t; mv t DM_v6.1.starch_gene.tpm.txt

#### functional annotation of allelic variations

ls vcf_genome_coordinates/*vcf | while read i
do
    java -Xmx20g -jar /vol1/agis/huangsanwen_group/lihongbo/software/snpEff/snpEff.jar DM_v6.1 $i > ${i}.snpeff
done

ls merge_vcf/*merge4HomoChr*.vcf.gz | while read i
do
    uncom=`echo $i | sed 's/\.gz//g'`
    gzip -d $i
    java -Xmx20g -jar /vol1/agis/huangsanwen_group/lihongbo/software/snpEff/snpEff.jar all_starch_gene $uncom > ${uncom}.snpeff
    gzip $uncom
done

##### number of snp indel sv in each cultivar
rm starch_snp_indel_sv_6_cultivar.xls
echo -e "Cultivar\tSNP\tInDel\tSV" > starch_snp_indel_sv_6_cultivar.xls
for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
do
	zcat merge_vcf/*${i}.vcf.gz | ./count_snp_indel_sv_from_vcf.py $i >> starch_snp_indel_sv_6_cultivar.xls
done


#### number of variants in 81 genes in 6 genomes

rm starch_gene_6_genomes_variant_number.xls
echo -e "Gene\tAltus\tAtlantic\tAvenger\tColomba\tPOR06V12\tSpunta\tAll" >> starch_gene_6_genomes_variant_number.xls
ls starch_gene_seqs/* | awk -F '/' '{print $2}' | awk -F '_' '{print $1}' | while read g
do
    echo -e "$g\t\c\n" >> starch_gene_6_genomes_variant_number.xls
    for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
    do
        num_v=0
        num_v=`zcat merge_vcf/${g}_merge4HomoChr_${i}.vcf.gz | grep -v '#' | wc -l`
        echo -e "$num_v\t\c\n" >> starch_gene_6_genomes_variant_number.xls
    done
    all_v=`grep -v '#' vcf_genome_coordinates/${g}_merge4HomoChr_6_genomes.vcf | wc -l`
    echo -e "$all_v" >> starch_gene_6_genomes_variant_number.xls
done

###### number of alleles in 6 functional classes of the 81 genes

rm starch_gene_6_genomes_allele_number_AllfunctionClass.xls
echo -e "Gene\tbi-allelic\ttri-allelic\ttetra-allelic\tgreater_four_allele" > starch_gene_6_genomes_allele_number_AllfunctionClass.xls
for n in upstream downstream intron splice synonymous nonsynonymous
do
    rm starch_gene_6_genomes_allele_number_functionClass_${n}.xls
    echo -e "Gene\tbi-allelic\ttri-allelic\ttetra-allelic\tgreater_four_allele" > starch_gene_6_genomes_allele_number_functionClass_${n}.xls
done

rm -rf allele_count; mkdir allele_count
ls starch_gene_seqs/* | awk -F '/' '{print $2}' | awk -F '_' '{print $1}' | while read g
do
    ./sta_allele_number.py vcf_genome_coordinates/${g}_merge4HomoChr_6_genomes.vcf.snpeff allele_count/${g}_upstream allele_count/${g}_downstream allele_count/${g}_intron allele_count/${g}_splice allele_count/${g}_synonymous allele_count/${g}_nonsynonymous allele_count/${g}_all $g
done

for n in upstream downstream intron splice synonymous nonsynonymous
do
    cat allele_count/*_${n} >> starch_gene_6_genomes_allele_number_functionClass_${n}.xls
done

cat allele_count/*_all >> starch_gene_6_genomes_allele_number_AllfunctionClass.xls

#### functional annotation of allelic variants in 81 genes in 6 genomes

rm starch_gene_6_genomes_variant_function*
echo -e "Gene\tUpstream\tDownstream\tUTR\tIntron\tSplice_region\tSynonymous\tNon-synonymous" >> starch_gene_6_genomes_variant_function_all.xls
for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
do
    echo -e "Gene\tUpstream\tDownstream\tUTR\tIntron\tSplice_region\tSynonymous\tNon-synonymous" >> starch_gene_6_genomes_variant_function_${i}.xls
done

for n in upstream downstream intron splice synonymous nonsynonymous
do
    echo -e "Gene\tAltus\tAtlantic\tAvenger\tColomba\tPOR06V12\tSpunta" >> starch_gene_6_genomes_variant_function_${n}.xls
done

ls starch_gene_seqs/* | awk -F '/' '{print $2}' | awk -F '_' '{print $1}' | while read g
do
    for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
    do
        num_up=`grep -v '#' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | grep upstream | wc -l`
        num_down=`grep -v '#' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | grep downstream | wc -l`
        num_intron=`grep -v '#' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | grep intron | wc -l`
		num_utr=`grep -v '#' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | grep UTR | wc -l`
        num_splice=`grep -v '#' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | grep splice | wc -l`
        num_syn=`grep -v '#' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | grep synonymous_variant | wc -l`
        num_nonsyn=`grep -v '#' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | egrep 'MODERATE|HIGH' | wc -l`
        echo -e "$g\t$num_up\t$num_down\t$num_utr\t$num_intron\t$num_splice\t$num_syn\t$num_nonsyn" >> starch_gene_6_genomes_variant_function_${i}.xls
    done

    for n in upstream downstream UTR intron splice synonymous nonsynonymous
    do
        echo -e "$g\t\c\n" >> starch_gene_6_genomes_variant_function_${n}.xls
    done
	
    for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
    do
        for n in upstream downstream UTR intron splice synonymous
        do
            num=`grep -v '#' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | grep $n | wc -l`
            echo -e "$num\t\c\n" >> starch_gene_6_genomes_variant_function_${n}.xls
        done
        num=`grep -v '#' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | egrep 'MODERATE|HIGH' | wc -l`
        echo -e "$num\t\c\n" >> starch_gene_6_genomes_variant_function_nonsynonymous.xls
    done

    for n in upstream downstream UTR intron splice synonymous nonsynonymous
    do
        echo -e "" >> starch_gene_6_genomes_variant_function_${n}.xls
    done
	
	num_up=`cat vcf_genome_coordinates/${g}_merge4HomoChr_6_genomes.vcf.snpeff | grep -v '#' | grep upstream | wc -l`
    num_down=`cat vcf_genome_coordinates/${g}_merge4HomoChr_6_genomes.vcf.snpeff | grep -v '#' | grep downstream | wc -l`
    num_intron=`cat vcf_genome_coordinates/${g}_merge4HomoChr_6_genomes.vcf.snpeff | grep -v '#' | grep intron | wc -l`
    num_utr=`cat vcf_genome_coordinates/${g}_merge4HomoChr_6_genomes.vcf.snpeff | grep -v '#' | grep UTR | wc -l`
    num_splice=`cat vcf_genome_coordinates/${g}_merge4HomoChr_6_genomes.vcf.snpeff | grep -v '#' | grep splice | wc -l`
    num_syn=`cat vcf_genome_coordinates/${g}_merge4HomoChr_6_genomes.vcf.snpeff | grep -v '#' | grep synonymous_variant | wc -l`
    num_nonsyn=`cat vcf_genome_coordinates/${g}_merge4HomoChr_6_genomes.vcf.snpeff | grep -v '#' | egrep 'MODERATE|HIGH' | wc -l`
    echo -e "$g\t$num_up\t$num_down\t$num_utr\t$num_intron\t$num_splice\t$num_syn\t$num_nonsyn" >> starch_gene_6_genomes_variant_function_all.xls

done

##### number of haplotypes in 81 genes in 6 genomes

rm starch_gene_6_genomes_haplotype_number.xls
echo -e "Gene\tAltus\tAtlantic\tAvenger\tColomba\tPOR06V12\tSpunta\tMean" >> starch_gene_6_genomes_haplotype_number.xls
ls starch_gene_seqs/* | awk -F '/' '{print $2}' | awk -F '_' '{print $1}' | while read g
do
    echo -e "$g\t\c\n" >> starch_gene_6_genomes_haplotype_number.xls
    for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
    do
        if [ -f merge_vcf/${g}_merge4HomoChr_${i}.vcf.gz ];then
			num_hap=`zcat merge_vcf/${g}_merge4HomoChr_${i}.vcf.gz | grep -v '##' | awk '{print NF-9}' | head -1`
        else
            num_hap=0
        fi
        echo -e "$num_hap\t\c\n" >> starch_gene_6_genomes_haplotype_number.xls

    done
    echo -e "Mean" >> starch_gene_6_genomes_haplotype_number.xls
done

awk 'BEGIN{OFS="\t"}{if($1~/Gene/){print $0}else{print $1,$2,$3,$4,$5,$6,$7,($2+$3+$4+$5+$6+$7)/6}}' starch_gene_6_genomes_haplotype_number.xls > tmp

mv tmp starch_gene_6_genomes_haplotype_number.xls

##### haplotype VS variant number correlation XY

rm starch_gene_6_genomes_haplotype_number_variant_number_XY.xls
for i in `seq 2 7`
do
   paste <(grep -v 'Gene' starch_gene_6_genomes_haplotype_number.xls | awk -v a=$i '{print $a}') <(grep -v 'Gene' starch_gene_6_genomes_variant_number.xls | awk -v a=$i '{print $a}')  >> starch_gene_6_genomes_haplotype_number_variant_number_XY.xls
done


##### number of variants in each gene haplotype of the 6 genoems
for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
do
    rm starch_gene_6_genomes_number_variant_each_haplotype_${i}.xls
    echo -e "Gene\tHaplotype1\tHaplotype2\tHaplotype3\tHaplotype4" > starch_gene_6_genomes_number_variant_each_haplotype_${i}.xls
done

ls starch_gene_seqs/* | awk -F '/' '{print $2}' | awk -F '_' '{print $1}' | while read g
do
    for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
    do
        echo -e "$g\t\c\n" >> starch_gene_6_genomes_number_variant_each_haplotype_${i}.xls
        if [ -f "merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff" ];then
            nf=`grep -v '##' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | awk '{print NF}' | head -1`
                for n in `seq 10 $nf`
                do
                    count=`grep -v '#' merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff | awk -v a=$n '{print $a}' | awk '$1!="./."&&$1!="0/0"' | wc -l`
                    echo -e "$count\t\c\n" >> starch_gene_6_genomes_number_variant_each_haplotype_${i}.xls
                done
				echo  $nf merge_vcf/${g}_merge4HomoChr_${i}.vcf.snpeff
				if [ $nf -lt 13 ];then
					for n in `seq $nf 12`
					do
						echo -e "NA\t\c\n" >> starch_gene_6_genomes_number_variant_each_haplotype_${i}.xls
					done
				fi
        else
		            nf=10
		            echo -e "0\t\c\n" >> starch_gene_6_genomes_number_variant_each_haplotype_${i}.xls
		        fi
		        echo -e "" >> starch_gene_6_genomes_number_variant_each_haplotype_${i}.xls
		    done
		done

###### transpose the matrix for plotting
for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
do
	./t.R starch_gene_6_genomes_number_variant_each_haplotype_${i}.xls starch_gene_6_genomes_number_variant_each_haplotype_${i}_t.xls
done

for i in Altus Atlantic Avenger Colomba POR06V12 Spunta
do
	./add_geneID_2_matrix.py starch_gene_6_genomes_number_variant_each_haplotype_${i}_t.xls starch_gene_ID.xls > t; mv t starch_gene_6_genomes_number_variant_each_haplotype_${i}_t.xls
done


##### number of three types of haplotypes

rm starch_gene_three_haplotypes.xls
ls vcf_genome_coordinates/*snpeff | while read i
do
	gene=`echo $i | awk -F '/' '{print $2}' | awk -F '_' '{print $1}'`
	paste <(grep -v '#' $i | cut -f '1-7') <(grep -v '#' $i | awk '{print $8}'|awk -F ',' '{print $1}') <(grep -v '#' $i |cut -f '9-999') > tmp
	amino_acid_hap=`cat <(grep '#CHROM' $i) <(cat tmp | grep $gene | grep -v '#' | egrep 'MODERATE|HIGH') | ./output_haplotype_from_vcf.py |awk '{print $2}'|sort|uniq| wc -l`
	transcript_hap=`cat <(grep '#CHROM' $i) <(cat tmp | grep $gene | grep -v '#' | grep -v 'upstream' | grep -v 'downstream') |./output_haplotype_from_vcf.py |awk '{print $2}'|sort|uniq|wc -l`
	regulation_hap=`cat <(grep '#CHROM' $i) <(cat tmp | grep $gene | grep -v '#') |./output_haplotype_from_vcf.py |awk '{print $2}'|sort|uniq|wc -l`
	echo -e "$gene\t${amino_acid_hap}\t${transcript_hap}\t${regulation_hap}" >> starch_gene_three_haplotypes.xls
done

###### domestication and haplotype numbers

./overlap_gene_with_sweeps.py sweep_region.v6.txt DM_v6.1_starch_pathway_gene.xls > DM_v6.1_starch_pathway_gene_sweep.xls
./overlap_gene_with_sweeps.py sweep_region.v6.txt DM_v6.1_SGA_pathway_gene.xls > DM_v6.1_SGA_pathway_gene_sweep.xls

rm starch_SGA_domestication_gene_three_haplotypes.xls
cat DM_v6.1_starch_pathway_gene_sweep.xls ../022_SGA_gene_allelic_variations/DM_v6.1_SGA_pathway_gene_sweep.xls | grep 'domestication' | awk '{print $1}' | while read i
do
	grep $i starch_SGA_gene_three_haplotypes.xls >> starch_SGA_domestication_gene_three_haplotypes.xls
done




