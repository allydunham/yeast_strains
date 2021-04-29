#!/usr/bin/env Rscript
# Regress Out RNA Level
source('src/config.R')

### Import Data ###
omics <- read_rds('data/rdata/omics.rds') %>%
  filter(!low_paff_range)

### Analyse ###
plots <- list()

tidy_lm <- function(x, type){
  s <- summary(x)
  tidy(x) %>%
    mutate(term = str_to_lower(term) %>% str_remove_all('[()]')) %>%
    rename_all(~str_replace_all(., '\\.', '_')) %>%
    pivot_wider(names_from = term, values_from = c(estimate, std_error, statistic, p_value), names_glue = "{term}_{.value}") %>%
    select(sort(tidyselect::peek_vars())) %>%
    bind_cols(tibble(type = type, r = sqrt(s$r.squared), r_squared = s$r.squared,
                     adj_r_squared = s$adj.r.squared, f_Statistic = s$fstatistic[1],
                     p_value = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)), .) %>%
    return()
}

calc_gene_lms <- function(tbl, ...){
  tbl <- drop_na(tbl)
  
  if (nrow(tbl) < 3){
    return(tibble(type=NA, r=NA, r_squared = NA, adj_r_squared=NA, f_Statistic=NA, p_value=NA,
                  intercept_estimate=NA, intercept_p_value=NA, intercept_statistic=NA, intercept_std_error=NA,
                  transcriptomic_estimate=NA, transcriptomic_p_value=NA, transcriptomic_statistic=NA,
                  transcriptomic_std_error=NA))
  }
  
  trans_lm <- lm(proteomic ~ transcriptomic, data = tbl)
  both_lm <- lm(proteomic ~ transcriptomic + paff, data = tbl)
  return(bind_rows(tidy_lm(trans_lm, type = 'Transcriptome'),
                   tidy_lm(both_lm, type = 'Transcriptome + P(aff)')))
}

proteome_paff_lm <- group_by(omics, systematic) %>%
  group_modify(calc_gene_lms) %>%
  filter(!is.na(type))

plots$paff_proteome_lm_r_squared <- select(proteome_paff_lm, systematic, type, adj_r_squared) %>%
  pivot_wider(names_from = type, values_from = adj_r_squared) %>%
  mutate(diff = `Transcriptome + P(aff)` - Transcriptome,
         big = ifelse(diff > 0.05, '> 0.05', '< 0.05')) %>%
  ggplot(aes(x = Transcriptome, y = `Transcriptome + P(aff)`, colour = big)) + 
  geom_point(shape = 20) +
  guides(colour=guide_legend(title = 'R<sup>2</sup> Difference')) +
  scale_colour_manual(values = c(`> 0.05`='firebrick2', `< 0.05`='cornflowerblue')) +
  labs(subtitle = 'Adj. R-Squared for LMs using different factors to predict Proteome FC') +
  theme(legend.title = element_markdown())

plots$correlation_distributions_r_squared <- select(proteome_paff_lm, systematic, type, adj_r_squared) %>%
  ggplot(aes(x = adj_r_squared, y = ..scaled.., colour = type)) +
  stat_density(geom = 'line', position = 'identity') +
  geom_vline(xintercept = 0, colour = 'grey', linetype = 'dotted') +
  scale_colour_manual(values = c(`Transcriptome + P(aff)`='firebrick2', Transcriptome='cornflowerblue')) +
  labs(x = 'Adj. R Squared', y = 'Scaled Density') +
  guides(colour = guide_legend(title = ''))

plots$paff_proteome_lm_r <- select(proteome_paff_lm, systematic, type, r) %>%
  pivot_wider(names_from = type, values_from = r) %>%
  mutate(diff = `Transcriptome + P(aff)` - Transcriptome,
         big = ifelse(diff > 0.05, '> 0.05', '< 0.05')) %>%
  ggplot(aes(x = Transcriptome, y = `Transcriptome + P(aff)`, colour = big)) + 
  geom_point(shape = 20) +
  guides(colour=guide_legend(title = 'R Difference')) +
  scale_colour_manual(values = c(`> 0.05`='firebrick2', `< 0.05`='cornflowerblue')) +
  labs(subtitle = 'Adj. R for LMs using different factors to predict Proteome FC') +
  theme(legend.title = element_markdown())

plots$correlation_distributions_r <- select(proteome_paff_lm, systematic, type, r) %>%
  ggplot(aes(x = r, y = ..scaled.., colour = type)) +
  stat_density(geom = 'line', position = 'identity') +
  geom_vline(xintercept = 0, colour = 'grey', linetype = 'dotted') +
  scale_colour_manual(values = c(`Transcriptome + P(aff)`='firebrick2', Transcriptome='cornflowerblue')) +
  labs(x = 'R', y = 'Scaled Density') +
  guides(colour = guide_legend(title = ''))

### Examples ###
paff_diff_genes <- select(proteome_paff_lm, systematic, type, adj_r_squared) %>%
  pivot_wider(names_from = type, values_from = adj_r_squared) %>%
  filter(`Transcriptome + P(aff)` - Transcriptome > 0.05) %>% 
  pull(systematic)

plots$high_diff_gene_cors <- (filter(omics, systematic %in% paff_diff_genes) %>%
                                select(systematic, strain, proteomic, transcriptomic, paff) %>%
                                pivot_longer(c(transcriptomic, paff), names_to = 'type', values_to = 'value') %>%
                                mutate(type = c(transcriptomic='Transcriptomic', paff='P(Aff)')[type]) %>%
                                ggplot(aes(x = value, y = proteomic, colour = type)) +
                                facet_grid(rows = vars(systematic), cols = vars(type), switch = 'x', scales = 'free') +
                                geom_point(shape = 20, show.legend = FALSE) +
                                scale_colour_brewer(type = 'qual') +
                                theme(strip.placement = 'outside')) %>%
  labeled_plot(units = 'cm', width = 20, height = length(paff_diff_genes) * 3, limitsize = FALSE)

### Save Plots ###
save_plotlist(plots, 'figures/regressions/', overwrite = "all")
