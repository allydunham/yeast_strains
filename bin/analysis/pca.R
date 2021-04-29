#!/usr/bin/env Rscript
# Determine PCAs for each data type
library(topGO)
library(org.Sc.sgd.db)
source('src/config.R')

### Import Data ###
phenotypes <- read_rds('data/rdata/bede_phenotypes.rds')
omics <- read_rds('data/rdata/omics.rds') %>%
  filter(!low_paff_range)
genetic_distance <- read_rds('data/rdata/genetic_distance_matrix.rds')

plots <- list()
### Calculate Dim Reduction ###
proteomic_pca <- select(omics, systematic, strain, proteomic) %>%
  drop_na() %>%
  pivot_wider(names_from = 'systematic', values_from = 'proteomic') %>%
  tibble_to_matrix(-strain, row_names = 'strain')
proteomic_pca[is.na(proteomic_pca)] <- 0
proteomic_pca <- prcomp(proteomic_pca, rank. = 100)

transcriptomic_pca <- select(omics, systematic, strain, transcriptomic) %>%
  drop_na() %>%
  pivot_wider(names_from = 'systematic', values_from = 'transcriptomic') %>%
  tibble_to_matrix(-strain, row_names = 'strain')
transcriptomic_pca[is.na(transcriptomic_pca)] <- 0
transcriptomic_pca <- prcomp(transcriptomic_pca, rank. = 100)

paff_pca <- select(omics, systematic, strain, paff) %>%
  drop_na() %>%
  complete(systematic, strain) %>%
  group_by(systematic) %>%
  mutate(paff = ifelse(is.na(paff), mean(paff, na.rm = TRUE), paff)) %>%
  ungroup() %>%
  pivot_wider(names_from = 'systematic', values_from = 'paff') %>%
  tibble_to_matrix(-strain, row_names = 'strain')
paff_pca <- prcomp(paff_pca, rank. = 100)

dist_pca <- prcomp(genetic_distance, rank. = 100)
  
omic_pcas <- as_tibble(proteomic_pca$x, rownames = 'strain') %>% 
  rename_with(~str_c('proteomic_', .), -strain) %>%
  full_join(phenotypes, ., by = 'strain') %>%
  full_join(as_tibble(transcriptomic_pca$x, rownames = 'strain') %>% rename_with(~str_c('transcriptomic_', .), -strain), by = 'strain') %>%
  full_join(as_tibble(paff_pca$x, rownames = 'strain') %>% rename_with(~str_c('paff_', .), -strain), by = 'strain') %>%
  full_join(as_tibble(dist_pca$x, rownames = 'strain') %>% rename_with(~str_c('dist_', .), -strain), by = 'strain') %>%
  select(strain, condition, score, qvalue, starts_with('proteomic_'), starts_with('transcriptomic_'), starts_with('paff_'))
write_rds(omic_pcas, 'data/rdata/omic_pcas.rds')

### Evaluate PCA Loadings ###
trans_clust <- hclust(dist(transcriptomic_pca$rotation))
trans_loadings <- as_tibble(transcriptomic_pca$rotation, rownames='systematic') %>%
  pivot_longer(-systematic, names_to = 'pc', values_to = 'loading') %>%
  mutate(systematic = factor(systematic, levels = trans_clust$labels[trans_clust$order]),
         pc = factor(pc, levels = str_c('PC', 1:50)))

plots$trans_pca_loadings <- ggplot(trans_loadings, aes(x = pc, y = systematic, fill = loading)) +
  geom_raster() +
  scale_fill_distiller(type = 'div', palette = 'RdBu') +
  guides(fill = guide_colourbar(title = 'Loading')) +
  theme(axis.title = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())
plots$trans_pca_loadings <- labeled_plot(plots$trans_pca_loadings, units='cm', height=100, width=50, file_format='jpg')

prot_clust <- hclust(dist(proteomic_pca$rotation))
prot_loadings <- as_tibble(proteomic_pca$rotation, rownames='systematic') %>%
  pivot_longer(-systematic, names_to = 'pc', values_to = 'loading') %>%
  mutate(systematic = factor(systematic, levels = prot_clust$labels[prot_clust$order]),
         pc = factor(pc, levels = str_c('PC', 1:50)))

plots$prot_pca_loadings <- ggplot(prot_loadings, aes(x = pc, y = systematic, file_format='jpg')) +
  geom_raster() +
  scale_fill_distiller(type = 'div', palette = 'RdBu') +
  guides(fill = guide_colourbar(title = 'Loading')) +
  theme(axis.title = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())
plots$prot_pca_loadings <- labeled_plot(plots$trans_pca_loadings, units='cm', height=100, width=50, file_format='jpg')

paff_clust <- hclust(dist(paff_pca$rotation))
paff_loadings <- as_tibble(paff_pca$rotation, rownames='systematic') %>%
  pivot_longer(-systematic, names_to = 'pc', values_to = 'loading') %>%
  mutate(systematic = factor(systematic, levels = paff_clust$labels[paff_clust$order]),
         pc = factor(pc, levels = str_c('PC', 1:50)))

plots$paff_pca_loadings <- ggplot(paff_loadings, aes(x = pc, y = systematic, fill = loading)) +
  geom_raster() +
  scale_fill_distiller(type = 'div', palette = 'RdBu') +
  guides(fill = guide_colourbar(title = 'Loading')) +
  theme(axis.title = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())
plots$paff_pca_loadings <- labeled_plot(plots$paff_pca_loadings, units='cm', height=100, width=50, file_format='jpg')

### GO Analysis ###
# go_map <- toTable(org.Sc.sgdGO) %>%
#   as_tibble() %>%
#   filter(Ontology == 'BP') %>%
#   group_by(systematic_name)
# go_map <- group_map(go_map, ~c(.x$go_id)) %>%
#   set_names(group_keys(go_map)$systematic_name)
# 
# go_analysis <- function(tbl, type, pc){
#   message(pc)
#   all_genes <- structure(abs(tbl$loading), names=tbl$systematic)
#   go_data <- new("topGOdata",
#                  ontology = "BP",
#                  allGenes = all_genes,
#                  geneSel = function(x){x[order(x, decreasing = TRUE)][1:100]},
#                  description = str_c("GO analysis of", type, pc, sep = ' '),
#                  annot = annFUN.gene2GO,
#                  gene2GO=go_map)
# 
#   fisher <- runTest(go_data, algorithm = "classic", statistic = "fisher")
#   ks <- runTest(go_data, algorithm = "classic", statistic = "ks")
#   ks_elim <- runTest(go_data, algorithm = "elim", statistic = "ks")
#   return(GenTable(go_data, fisher = fisher, ks = ks, ks_elim = ks_elim,
#                   orderBy = "ks_elim", ranksOf = "fisher", topNodes = 10))
# }
# 
# map_over_pcs <- function(tbl, type){
#   message(type)
#   tbl <- group_by(tbl, pc)
#   group_map(tbl, ~go_analysis(.x, type, .y)) %>%
#     set_names(group_keys(tbl)$pc)
# }
# 
# top_gene_go <- bind_rows(Transcriptomic=trans_loadings, Proteomic=prot_loadings, `P(Aff)`=paff_loadings, .id = 'type') %>%
#   mutate(systematic = as.character(systematic)) %>%
#   filter(as.integer(str_remove(pc, 'PC')) <= 10) %>%
#   group_by(type)
# top_gene_go <- group_map(top_gene_go, ~map_over_pcs(.x, .y)) %>%
#   set_names(group_keys(top_gene_go)$type)
# write_rds(top_gene_go, 'data/rdata/top_gene_go.rds')
top_gene_go <- read_rds('data/rdata/top_gene_go.rds')

### Save Plots ###
save_plotlist(plots, 'figures/pcas/', overwrite = "all")
