#!/usr/bin/env Rscript
library(VariantAnnotation)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

args = commandArgs(trailingOnly=TRUE)

vcf <- readVcf(args[1], "sacCer3")
seqlevels(vcf) <- sapply(seqlevels(vcf), function(x){paste0('chr',as.roman(gsub('(C|c)hromosome|(C|c)hr','',x)))})
  
sc <- BSgenome.Scerevisiae.UCSC.sacCer3
#seqnames(sc) <- paste0('chromosome',c(1:16,'M'))
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
#seqlevels(txdb) <- paste0('chromosome',c(1:16,'M'))

coding <- data.frame(predictCoding(vcf, txdb, seqSource=sc))

## Function to get the correct alt allele based on a coding mutant row
getALT <- function(x){
  if (x$strand == '-'){
    return(as.character(reverseComplement(DNAString(x$varAllele))))
  } else {
    return(x$varAllele)
  }
}

coding$ALT <- apply(coding, 1, getALT)

## Fix coding for mutfunc
coding <- coding[coding$width == 1 & !is.na(coding$ALT),]
coding$mut_id <- paste0(coding$seqnames,':',coding$start,'_',coding$REF,'/',coding$ALT)
coding <- coding[c('mut_id', 'REFCODON', 'VARCODON', 'PROTEINLOC',
                   'REFAA', 'VARAA', 'GENEID', 'CONSEQUENCE')]

colnames(coding) <- c('mut_id','ref_codon','alt_codon','pos_aa',
                      'ref_aa','alt_aa','gene','type')

write.table(coding, args[2], sep='\t', row.names=F, quote=F)
