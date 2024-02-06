#!/usr/bin/env Rscript

library('PopGenome')
argv<-commandArgs(TRUE)

GENOME.class <- readVCF(argv[1],numcols=1000, tid=argv[2], from=1, to=argv[3], approx=FALSE, include.unknown=TRUE, parallel=FALSE, gffpath=argv[4])

#### population assignment, each accession each line
p1 <- as.character(read.table("pop1_starch.txt")[[1]])
p2 <- as.character(read.table("pop2_starch.txt")[[1]])
#p3 <- as.character(read.table("pop3.txt")[[1]])


slide <- sliding.window.transform(GENOME.class,as.numeric(argv[5]),100, type=2)

#length(slide@region.names)
slide <- diversity.stats(slide)
nucdiv <- slide@nuc.diversity.within
nucdiv <- nucdiv/(as.numeric(argv[5]))

#slide <- neutrality.stats(slide, FAST=FALSE)
#tajimaD <- slide@Tajima.D

GENOME.class <- set.populations(GENOME.class,list(p1, p2), diploid=FALSE)
slide <- sliding.window.transform(GENOME.class,as.numeric(argv[5]),100, type=2)

slide <- F_ST.stats(slide, mode="nucleotide")
pairwise.FST <- t(slide@nuc.F_ST.pairwise)

slide <- neutrality.stats(slide, FAST=FALSE)
tajimaD <- slide@Tajima.D

write.table(nucdiv,argv[6],row.names = T,col.names = T,quote = FALSE,sep='\t')
write.table(pairwise.FST,argv[7],row.names = T,col.names = T,quote = FALSE,sep='\t')
write.table(tajimaD,argv[8],row.names = T,col.names = T,quote = FALSE,sep='\t')

