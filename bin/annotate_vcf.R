#!/usr/bin/env Rscript
library(VariantAnnotation)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

args = commandArgs(trailingOnly=TRUE)

vcf <- readVcf(args[1], "sacCer3", param=ScanVcfParam(samples=NA))
message('VCF Read')
seqlevels(vcf) <- sapply(seqlevels(vcf), function(x){paste0('chr',as.roman(gsub('(C|c)hromosome|(C|c)hr','',x)))})
  
sc <- BSgenome.Scerevisiae.UCSC.sacCer3
#seqnames(sc) <- paste0('chromosome',c(1:16,'M'))
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
#seqlevels(txdb) <- paste0('chromosome',c(1:16,'M'))

coding <- data.frame(predictCoding(vcf, txdb, seqSource=sc))
coding$ALT <- sapply(coding$ALT, function(x){paste(as.character(x), collapse = ',')})
coding$PROTEINLOC <- sapply(coding$PROTEINLOC, function(x){paste(x, collapse = ',')})
write.table(coding, stdout(), sep='\t', row.names=F, quote=F)
