#!/usr/bin/env Rscript
# Initial analysis script

source('src/config.R')

proteomic <- read_tsv('data/raw/1k_quant_wide.tsv') %>%
  rename(gene = X1) %>%
  pivot_longer(-gene, names_to = 'strain', values_to = 'abundance')

transcriptomic <- read_csv('data/raw/tpm_FinalSet_969strains.csv') %>%
  rename(uniprot = systematic_name) %>%
  pivot_longer(-gene, names_to = 'strain', values_to = 'abundance')
