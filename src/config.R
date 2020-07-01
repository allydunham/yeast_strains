#!/usr/bin/env Rscript
# Setup namespace for Yeast Strains project

library(tidyverse)

# Import Uniprot to Gene map
NAME_TO_UNIPROT <- read_tsv('meta/uniprot_map', col_names = c('uniprot', 'database', 'id')) %>%
    filter(database == 'EnsemblGenome') %>%
    select(uniprot, id) %>%
    drop_na()
NAME_TO_UNIPROT <- structure(NAME_TO_UNIPROT$uniprot, names = NAME_TO_UNIPROT$id)
