#!/usr/bin/env Rscript
# Identify Early Stop and Frameshift locations
source('src/config.R')

# Rscript bin/paff.R data/variants/chr1.tsv data/variants/chr1.genotypes
args = commandArgs(trailingOnly=TRUE)

# Import Data
gene_lengths <- read_tsv('meta/yeast_gene_lengths') %>%
  rename_all(~str_to_lower(.) %>% str_replace_all(' ', '_')) %>%
  select(uniprot = entry, length)

annotation <- read_tsv(args[1], col_types = cols(PROTEINLOC=col_character(), CDSID=col_character())) %>%
  rename_all(str_to_lower) %>%
  mutate(uniprot = unname(NAME_TO_UNIPROT[geneid])) %>%
  filter(consequence %in% c('nonsense', 'frameshift'))

genotypes <- read_tsv(args[2]) %>%
  pivot_longer(!chromosome:alt, names_to = 'strain', values_to = 'genotype') %>%
  filter(!genotype == 0)

comb <- left_join(select(annotation, systematic=geneid, uniprot, chromosome = seqnames, position = start, proteinloc, consequence, ref, alt, refaa, varaa),
                        select(genotypes, -genotype), by=c('chromosome', 'position', 'ref', 'alt')) %>%
  left_join(gene_lengths, by = "uniprot") %>%
  mutate(protein_pos = as.numeric(str_remove(proteinloc, ',.*')),
         prop = protein_pos / length) %>%
  select(strain, systematic, uniprot, consequence, everything()) %>%
  drop_na(strain) %>%
  arrange(strain, chromosome, position)

write(format_tsv(comb), file = stdout())
