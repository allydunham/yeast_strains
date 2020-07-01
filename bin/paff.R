#!/usr/bin/env Rscript
# Generate P(Aff) type scores
source('src/config.R')

sift <- read_tsv('data/raw/sift_filtered.tab')

genotypes <- read_tsv('data/variants/chr1.genotypes') %>%
  pivot_longer(!one_of(c('chromosome', 'position', 'ref', 'alt')), names_to = 'strain', values_to = 'genotype') %>%
  filter(!genotype == 0)

annotation <- read_tsv('data/variants/chr1.tsv') %>%
  rename_all(str_to_lower) %>%
  mutate(uniprot = unname(NAME_TO_UNIPROT[geneid]))

sift <- filter(sift, acc %in% annotation$uniprot) %>%
  mutate(pos = as.character(pos))

annotation <- left_join(annotation, select(sift, uniprot=acc, proteinloc=pos, refaa=ref, varaa=alt, score),
                        by = c("uniprot", "proteinloc", "refaa", "varaa"))

View(filter(annotation, !consequence == 'synonymous'))
