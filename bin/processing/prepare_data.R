#!/usr/bin/env Rscript
# Process data for analyses
source('src/config.R')

# Strain information
strains <- read_tsv('meta/strain_information.tsv')

# Paff Scores
paff <- read_tsv('data/paff_scores.tsv') %>%
  rename(systematic=geneid)
# TODO Strains CEN.PK and Reference? called XTRA_DXL and FY4-6 in their data
write_rds(paff, 'data/rdata/paff.rds')

# Identify genes with poor Paff distributions
low_paff_range_genes <- group_by(paff, systematic) %>%
  filter(max(paff) - min(paff) < 0.5 | max(paff) < 0.33 | min(paff) > 0.66) %>%
  ungroup() %>%
  pull(systematic) %>%
  unique()
write_rds(low_paff_range_genes, 'data/rdata/low_paff_range.rds')

# Genes with Early Stop 
early_stops <- read_tsv('data/early_stops.tsv') %>%
  filter(systematic %in% paff$systematic) %>%
  mutate(low_paff_range = systematic %in% low_paff_range_genes)
write_rds(early_stops, 'data/rdata/early_stops.rds')

# Omics
proteomic <- read_csv('data/raw/1k_quant_wide_systematic_name.csv', ) %>%
  select(-X1) %>%
  rename(gene=symbol, systematic=systematic_name) %>%
  pivot_longer(!one_of(c('gene', 'systematic')), names_to = 'strain', values_to = 'abundance') %>% 
  group_by(gene, systematic) %>%
  mutate(fc = abundance / median(abundance),
         fc = log2(fc + min(fc[fc > 0]))) %>%
  ungroup() %>%
  filter(systematic %in% paff$systematic)

transcriptomic <- read_csv('data/raw/tpm_FinalSet_969strains.csv') %>%
  rename(systematic = systematic_name) %>%
  pivot_longer(-systematic, names_to = 'strain', values_to = 'abundance') %>% 
  group_by(systematic) %>%
  mutate(fc = abundance / median(abundance),
         fc = log2(fc + min(fc[fc > 0]))) %>%
  ungroup() %>%
  filter(systematic %in% paff$systematic)

comb <- full_join(select(proteomic, -gene) %>% rename(proteomic_raw=abundance, proteomic=fc),
                  rename(transcriptomic, transcriptomic_raw=abundance, transcriptomic=fc),
                  by = c('systematic', 'strain')) %>%
  left_join(paff, by = c('systematic', 'strain')) %>%
  mutate(low_paff_range = systematic %in% low_paff_range_genes)
write_rds(comb, 'data/rdata/omics.rds')

# Gene Annotation
annotation <- map(1:16, ~read_tsv(str_c('data/variants/chr', ., '.tsv'), col_types = cols(PROTEINLOC=col_character(), CDSID=col_character()))) %>%
  bind_rows() %>%
  rename_all(str_to_lower) %>%
  rename(systematic = geneid) %>%
  mutate(uniprot = unname(NAME_TO_UNIPROT[systematic]),
         low_paff_range = systematic %in% low_paff_range_genes)
write_rds(annotation, 'data/rdata/annotation.rds')

# Phenotypes
liti_phenotypes <- read_tsv('data/raw/phenoMatrix_35ConditionsNormalizedByYPD.tab') %>%
  rename(strain = X1) %>%
  pivot_longer(-strain, names_to = 'condition', values_to = 'score')
write_rds(liti_phenotypes, 'data/rdata/liti_phenotypes.rds')

bede_phenotypes <- read_tsv('data/raw/yeasts_natural.tsv', col_types = cols(strain=col_character())) %>%
  filter(!info == 'undefined') %>%
  select(strain=info, os_strain=strain, everything())
write_rds(bede_phenotypes, 'data/rdata/bede_phenotypes.rds')
