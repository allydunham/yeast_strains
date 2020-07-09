#!/usr/bin/env Rscript
# Setup namespace for Yeast Strains project
library(tidyverse)
library(ggpubr)
library(ggtext)
library(broom)
library(plotlistr)

# Import Uniprot to Gene map
NAME_TO_UNIPROT <- read_tsv('meta/uniprot_map', col_names = c('uniprot', 'database', 'id')) %>%
    filter(database == 'EnsemblGenome') %>%
    select(uniprot, id) %>%
    drop_na()
NAME_TO_UNIPROT <- structure(NAME_TO_UNIPROT$uniprot, names = NAME_TO_UNIPROT$id)

# GGPlot Theme
theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))