#!/usr/bin/env Rscript
# Assess impact of each gene on each condition
library(biomaRt)
source('src/config.R')

### Import Data ###
phenotypes <- read_rds('data/rdata/bede_phenotypes.rds') %>%
  select(strain, condition, score, qvalue) %>%
  mutate(del = qvalue < 0.05 & score < 0)

omics <- read_rds('data/rdata/omics.rds') %>%
  filter(!low_paff_range) %>%
  select(systematic, strain, proteomic, transcriptomic, paff)

kos <- read_tsv('data/raw/ko_scores.txt') %>%
  select(-position) %>%
  filter(!gene == 'WT')

uniprot_map <- read_tsv('meta/uniprot_map', col_names = c('uniprot', 'type', 'value')) %>%
  filter(type == 'Gene_OrderedLocusName') %>%
  select(uniprot, systematic=value)

### Calculate LMs ###
analyse_gene_con <- function(tbl, key){
  if (nrow(tbl) < 3){
    return(tibble(rsquared = NA, adj_rsquared = NA, fstatistic = NA, pvalue = NA, estimate_intercept = NA,
                  estimate_proteomic = NA, estimate_transcriptomic = NA, estimate_paff = NA, std_error_intercept = NA,
                  std_error_proteomic = NA, std_error_transcriptomic = NA, std_error_paff = NA, statistic_intercept = NA,
                  statistic_proteomic = NA, statistic_transcriptomic = NA, statistic_paff = NA, p_value_intercept = NA,
                  p_value_proteomic = NA, p_value_transcriptomic = NA, p_value_paff = NA))
  }
  
  m <- lm(score ~ proteomic + transcriptomic + paff, data = tbl)
  s <- summary(m)
  t <- tidy(m)
  t$term[t$term == '(Intercept)'] <- 'intercept'
  pivot_wider(t, names_from = term, values_from = c(-term)) %>%
    rename_all(~str_replace_all(., '\\.', '_')) %>%
    mutate(rsquared = s$r.squared, adj_rsquared = s$adj.r.squared, fstatistic = s$fstatistic[1], pvalue = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)) %>%
    select(rsquared, adj_rsquared, fstatistic, pvalue, everything()) %>%
    return()
}

analyse_gene <- function(tbl, key){
  left_join(tbl, phenotypes, by = 'strain') %>%
    drop_na() %>%
    group_by(condition) %>%
    group_modify(analyse_gene_con) %>%
    return()
}

# genes <- group_by(omics, systematic) %>%
#   group_modify(analyse_gene)
# write_rds(genes, 'data/rdata/per_gene_lms.rds')
genes <- read_rds('data/rdata/per_gene_lms.rds')

genes_processed <- drop_na(genes) %>%
  mutate(padj = p.adjust(pvalue, 'fdr'))

### Plots ###
plots <- list()

plots$overview <- (ggplot(genes_processed, aes(x = rsquared, y = padj, colour = padj < 0.001, label = systematic)) +
  facet_wrap(~condition, labeller = label_wrap_gen(15)) +
  geom_point(shape = 20) +
  geom_text_repel(data = filter(genes_processed, rsquared > 0.2), colour = 'black', nudge_x = 0.2, nudge_y = 0.5, size = 2) +
  scale_colour_manual(name = 'p<sub>adj</sub>', values = c(`TRUE`='red', `FALSE`='black'), labels = c(`TRUE`='< 0.001', `FALSE`='> 0.001')) +
  labs(x = expression(R^2), y = expression(p[adj])) + 
  theme(text = element_text(size=8),
        legend.title = element_markdown())) %>%
  labeled_plot(units = 'cm', height = 29, width = 21)

p_value_cutoff <- 0.01
factor_summary <- mutate(genes_processed,
       trans = p.adjust(p_value_transcriptomic, 'fdr') < p_value_cutoff,
       prot = p.adjust(p_value_proteomic, 'fdr') < p_value_cutoff,
       paff = p.adjust(p_value_paff, 'fdr') < p_value_cutoff) %>%
  group_by(condition) %>%
  summarise(Transcriptomic = sum(trans), Proteomic = sum(prot), `P(Aff)` = sum(paff),
            `Transcriptomic & Proteomic` = sum(trans & prot), `Transcriptomic & P(Aff)` = sum(trans & paff), `Proteomic & P(Aff)` = sum(prot & trans),
            All = sum(trans & prot & paff), .groups = 'drop') %>%
  pivot_longer(-condition, names_to = 'factors', values_to = 'count') %>%
  mutate(n_factors = c(Transcriptomic = 1, Proteomic = 1, `P(Aff)` = 1, `Transcriptomic & Proteomic` = 2,
                       `Transcriptomic & P(Aff)` = 2, `Proteomic & P(Aff)` = 2, All = 3)[factors],
         factors = factor(factors, levels = c('All', 'Proteomic & P(Aff)', 'Transcriptomic & P(Aff)', 'Transcriptomic & Proteomic',
                                              'P(Aff)', 'Proteomic', 'Transcriptomic')))
plots$factor_combs <- (ggplot(factor_summary, aes(x = factors, y = count, fill = factor(n_factors))) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~condition) +
  scale_fill_brewer(type = 'qual', palette = 'Set2') +
  coord_flip() +
  labs(x = '', y = str_c('Count of Genes where p<sub>adj</sub> < ', p_value_cutoff)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted'),
        axis.title.x = element_markdown())) %>%
  labeled_plot(units = 'cm', height = 30, width = 30)

omics_high_r_squared <- filter(omics, systematic %in% (filter(genes_processed, rsquared > 0.25) %>% pull(systematic) %>% unique())) %>%
  left_join(phenotypes, by = 'strain') %>% 
  left_join(genes_processed, by = c('systematic', 'condition')) %>%
  filter(rsquared > 0.25) %>%
  select(systematic, strain, condition, score, qvalue, del, transcriptomic, proteomic, paff) %>%
  drop_na() %>%
  mutate(del = ifelse(del, 'q < 0.05', 'q > 0.05'))

plots$high_r_paff <- (ggplot(omics_high_r_squared, aes(x = paff, y = score, colour = del)) +
                        geom_point() +
                        geom_smooth(method = 'lm', formula = y ~ x, mapping = aes(x = transcriptomic, y = score), inherit.aes = FALSE) +
                        facet_wrap(~ systematic + condition, ncol = 5) +
                        scale_colour_manual(name = '', values = c(`q < 0.05` = 'red', `q > 0.05` = 'black')) +
                        xlim(c(0, 1)) +
                        labs(x = 'P(Aff)', y = 'S-Score')) %>%
  labeled_plot(units = 'cm', width = 120, height = 120)

plots$high_r_prot <- (ggplot(omics_high_r_squared, aes(x = proteomic, y = score, colour = del)) +
                        geom_point() +
                        geom_smooth(method = 'lm', formula = y ~ x, mapping = aes(x = transcriptomic, y = score), inherit.aes = FALSE) +
                        facet_wrap(~ systematic + condition, ncol = 5) +
                        scale_colour_manual(name = '', values = c(`q < 0.05` = 'red', `q > 0.05` = 'black')) +
                        labs(x = 'Proteomic FC', y = 'S-Score')) %>%
  labeled_plot(units = 'cm', width = 120, height = 120)

plots$high_r_trans <- (ggplot(omics_high_r_squared, aes(x = transcriptomic, y = score, colour = del)) +
                        geom_point() +
                        geom_smooth(method = 'lm', formula = y ~ x, mapping = aes(x = transcriptomic, y = score), inherit.aes = FALSE) +
                        facet_wrap(~ systematic + condition, ncol = 5) +
                        scale_colour_manual(name = '', values = c(`q < 0.05` = 'red', `q > 0.05` = 'black')) +
                        labs(x = 'Transcriptomic FC', y = 'S-Score')) %>%
  labeled_plot(units = 'cm', width = 120, height = 120)

## Analyse genes ##
common_genes <- filter(genes_processed, adj_rsquared > 0.1, padj < 0.01) %>% 
  count(systematic) %>%
  arrange(desc(n))
# Enrichmed for basic metabolism type processes at r > 0.05:
# Enrichment for r > 0.1, n > 5
#  - NADH oxidation
#  - alpha-amino acid biosynthetic process
#  - energy derivation by oxidation of organic compounds
#  - nucleobase-containing small molecule metabolic process
#  - mitochondrion organization
# No enrichment for r > 0.25

plots$gene_count_hist <- ggplot(common_genes, aes(x = n)) +
  geom_histogram(binwidth = 1, fill = 'cornflowerblue') +
  labs(x = 'Number of significant conditions (Adj. r<sup>2</sup> > 0.1, p<sub>adj</sub> < 0.01)', y = 'Count') +
  theme(axis.title.x = element_markdown())

per_con_top_genes <- filter(genes_processed, adj_rsquared > 0.1, padj < 0.01) %>% 
  group_by(condition) %>%
  {
    n <- group_keys(.)$condition
    group_map(., ~.$systematic) %>%
    set_names(n)
  }
# Enrichments generally are broad biosynthesis/metabolism type stuff
# Possible interesting cases:
# * 39ºC (72H) - protein refoldxing
# * Anerobic growth - sulphate assimilation?
# 

# Compare to KOs
ko_cors <- select(kos, strain, condition, systematic = gene, score) %>%
  distinct(strain, condition, systematic, .keep_all = TRUE) %>%
  left_join(genes_processed, by = c("condition", "systematic")) %>%
  drop_na(rsquared) %>%
  group_by(strain, condition) %>%
  group_modify(~tidy(cor.test(.$score, .$rsquared))) %>%
  mutate(padj = p.adjust(p.value, method = 'fdr'),
         cat = ifelse(padj <= 0.05, 'p<sub>adj</sub> &le; 0.05', 'p<sub>adj</sub> &gt; 0.05'))

plots$ko_score_cors <- (ggplot(ko_cors, aes(x = estimate, y = -log10(padj), colour = cat, label = condition)) +
  facet_wrap(~strain) +
  geom_point() +
  geom_text_repel(data = filter(ko_cors, padj <= 0.05), colour='black') +
  scale_colour_brewer(type = 'qual', palette = 'Set1', direction = -1, name = '') +
  labs(x = 'Correlation Coefficient between KO S-Score and R-Squared', y = '-log<sub>10</sub> p<sub>adj</sub>') +
  theme(legend.text = element_markdown(),
        axis.title.y = element_markdown())) %>%
  labeled_plot(units = 'cm', height = 20, width = 20)

# KO Jaccard
calc_jaccard <- function(set1, set2){
  length(intersect(set1, set2))/length(union(set1, set2))
}

calc_prop <- function(lm, ko){
  length(intersect(lm, ko))/length(lm)
}

lm_sig_genes <- filter(genes_processed, padj < 0.01) %>%
  group_by(condition) %>%
  {
    n <- group_keys(.)$condition
    group_map(., ~.$systematic) %>%
      set_names(n)
  }
  

jaccard <- filter(kos, qvalue < 0.01) %>% 
  select(strain, condition, gene) %>%
  group_by(strain, condition) %>%
  summarise(jaccard = calc_jaccard(gene, lm_sig_genes[[condition[1]]]),
            prop_lm = calc_prop(lm_sig_genes[[condition[1]]], gene),
            .groups='drop')

plots$ko_jaccard <- ggplot(jaccard, aes(x = condition, y = jaccard, fill = strain)) +
  geom_col(position = 'dodge') +
  coord_flip() +
  labs(x = '', y = 'Jaccard Index') +
  scale_fill_brewer(type='qual', palette = 'Dark2', name='Strain') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour='grey', linetype='dotted'))

plots$ko_overlap_per_strain <- ggplot(jaccard, aes(x = condition, y = prop_lm, fill = strain)) +
  geom_col(position = 'dodge') +
  coord_flip() +
  labs(x = '', y = 'Proportion of significant genes with significant KO phenotype') +
  scale_fill_brewer(type='qual', palette = 'Dark2', name='Strain') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour='grey', linetype='dotted'))

total_prop <- filter(kos, qvalue < 0.01) %>%
  select(condition, gene) %>%
  distinct() %>%
  group_by(condition) %>%
  summarise(prop_lm = calc_prop(lm_sig_genes[[condition[1]]], gene), .groups='drop')

plots$ko_overlap <- ggplot(total_prop, aes(x = condition, y = prop_lm)) +
  geom_col(position = 'dodge', fill = 'cornflowerblue') +
  coord_flip() +
  labs(x = '', y = 'Proportion of significant genes with significant KO phenotype') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour='grey', linetype='dotted'))

# Top genes per condition
top_genes <- genes_processed %>%
  filter(padj < 0.01) %>%
  group_by(condition) %>%
  slice_max(adj_rsquared, n = 10) %>%
  left_join(uniprot_map, by = 'systematic') %>%
  select(uniprot, systematic, everything())

ensembl = useEnsembl(biomart="ensembl", dataset="scerevisiae_gene_ensembl")

descriptions <- getBM(attributes=c('uniprotswissprot','description'), filters = 'uniprotswissprot',
                      values = unique(top_genes$uniprot), mart=ensembl) %>%
  as_tibble() %>%
  rename(uniprot=uniprotswissprot)

top_genes <- select(top_genes, condition, uniprot, systematic, adj_rsquared, padj) %>%
  left_join(descriptions, by = 'uniprot')
write_tsv(top_genes, 'data/top_lm_genes.tsv')

# Common KO genes
common_ko_genes <- filter(kos, score < 0, qvalue < 0.01) %>%
  distinct(strain, condition, gene) %>%
  count(condition, gene) %>%
  filter(n == 4) %>%
  left_join(uniprot_map, by = c(gene='systematic'))

ko_descriptions <- getBM(attributes=c('uniprotswissprot','description'), filters = 'uniprotswissprot',
                         values = unique(common_ko_genes$uniprot), mart=ensembl) %>%
  as_tibble() %>%
  rename(uniprot=uniprotswissprot)

common_ko_genes <- left_join(common_ko_genes, ko_descriptions, by = 'uniprot')
write_tsv(common_ko_genes, 'data/common_ko_genes.tsv')

### Save plots ###
save_plotlist(plots, 'figures/per_gene_lms', overwrite = "all")

