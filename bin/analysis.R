#!/usr/bin/env Rscript
# Initial analysis script
source('src/config.R')

proteomic <- read_csv('data/raw/1k_quant_wide_systematic_name.csv', ) %>%
  select(-X1) %>%
  rename(gene=symbol, uniprot=systematic_name) %>%
  pivot_longer(!one_of(c('gene', 'uniprot')), names_to = 'strain', values_to = 'abundance')

transcriptomic <- read_csv('data/raw/tpm_FinalSet_969strains.csv') %>%
  rename(uniprot = systematic_name) %>%
  pivot_longer(-uniprot, names_to = 'strain', values_to = 'abundance')


