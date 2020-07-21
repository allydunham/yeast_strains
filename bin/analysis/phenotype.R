#!/usr/bin/env Rscript
# Analyse Phenotypes
source('src/config.R')

### Import Data ###
phenotypes <- read_rds('data/rdata/phenotypes.rds')
omics <- read_rds('data/rdata/omics.rds') %>%
  filter(!low_paff_range)

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

omic_pcas <- as_tibble(proteomic_pca$x, rownames = 'strain') %>% 
  rename_with(~str_c('proteomic_', .), -strain) %>%
  full_join(phenotypes, ., by = 'strain') %>%
  full_join(as_tibble(transcriptomic_pca$x, rownames = 'strain') %>% rename_with(~str_c('transcriptomic_', .), -strain), by = 'strain') %>%
  full_join(as_tibble(paff_pca$x, rownames = 'strain') %>% rename_with(~str_c('paff_', .), -strain), by = 'strain')

### Generate Linear Models ###
tidy_lm <- function(x, type){
  s <- summary(x)
  return(tibble(type = type, adj_r_squared = s$adj.r.squared, f_Statistic = s$fstatistic[1],
                p_value = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)))
}

calc_lms <- function(tbl, ...){
  tbl <- drop_na(tbl)
  
  if (nrow(tbl) < 3){
    return(tibble(type=NA, adj_r_squared=NA, f_Statistic=NA, p_value=NA))
  }
  
  prot_lm <- lm(sscore ~ ., data = select(tbl, sscore, starts_with('proteomic')))
  trans_lm <- lm(sscore ~ ., data = select(tbl, sscore, starts_with('transcriptomic')))
  paff_lm <- lm(sscore ~ ., data = select(tbl, sscore, starts_with('paff')))
  
  prot_trans_lm <- lm(sscore ~ ., data = select(tbl, sscore, starts_with('proteomic'), starts_with('transcriptomic')))
  prot_paff_lm <- lm(sscore ~ ., data = select(tbl, sscore, starts_with('proteomic'), starts_with('paff')))
  trans_paff_lm <- lm(sscore ~ ., data = select(tbl, sscore, starts_with('transcriptomic'), starts_with('paff')))
  
  comb_lm <- lm(sscore ~ ., data = select(tbl, -strain))
  
  return(
    bind_rows(
      tidy_lm(prot_lm, type = 'Proteomic'),
      tidy_lm(trans_lm, type = 'Transcriptomic'),
      tidy_lm(paff_lm, type = 'P(Aff)'),
      tidy_lm(prot_trans_lm, type = 'Proteomic/Transcriptomic'),
      tidy_lm(prot_paff_lm, type = 'Proteomic/P(Aff)'),
      tidy_lm(trans_paff_lm, type = 'Transcriptomic/P(Aff)'),
      tidy_lm(comb_lm, type = 'All'),
    )
  )
}

phenotype_lms <- group_by(omic_pcas, condition) %>%
  group_modify(calc_lms) %>%
  drop_na() %>%
  mutate(type = factor(type, levels = c('P(Aff)', 'Transcriptomic', 'Proteomic', 'Transcriptomic/P(Aff)',
                                        'Proteomic/P(Aff)', 'Proteomic/Transcriptomic', 'All')))
  
  
### Analyse ###
lm_colours <- c(`P(Aff)`='yellow', `Transcriptomic`='magenta', `Proteomic`='cyan',
                `Transcriptomic/P(Aff)`='red', `Proteomic/P(Aff)`='green', `Proteomic/Transcriptomic`='blue',
                `All` = 'black')
plots$lm_summary <- ggplot(phenotype_lms, aes(x = type, y = adj_r_squared, fill = type)) +
  geom_col(show.legend = FALSE) +
  coord_flip() + 
  scale_fill_manual(values = lm_colours) +
  facet_wrap(~condition, nrow = 6, ncol = 6) +
  theme(panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted'),
        panel.grid.major.y = element_blank())
plots$lm_summary <- labeled_plot(plots$lm_summary, units = 'cm', height = 25, width = 25)

### Save Plots ###
save_plotlist(plots, 'figures/phenotypes/')
