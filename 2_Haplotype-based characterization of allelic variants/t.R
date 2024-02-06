#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)
a <- read.table(argv[1], header=T, row.names=1)
#row.names(a) <- a$Gene
#a <- a[,-1]
b <- data.frame(t(a))
write.table(b,argv[2],row.names = T,col.names = T,quote = FALSE,sep='\t')
