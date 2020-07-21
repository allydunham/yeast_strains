#!/usr/bin/env Rscript
# Initial analysis script
source('src/config.R')

### Import Data ###
paff <- read_tsv('data/paff_scores.tsv') %>%
  rename(systematic=geneid)
# TODO Strains CEN.PK and Reference? called XTRA_DXL and FY4-6 in their data

early_stops <- read_tsv('data/early_stops.tsv') %>%
  filter(systematic %in% paff$systematic)

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
  left_join(paff, by = c('systematic', 'strain'))

annotation <- map(1:16, ~read_tsv(str_c('data/variants/chr', ., '.tsv'), col_types = cols(PROTEINLOC=col_character(), CDSID=col_character()))) %>%
  bind_rows() %>%
  rename_all(str_to_lower) %>%
  mutate(uniprot = unname(NAME_TO_UNIPROT[geneid])) %>%
  rename(systematic = geneid)

plots <- list()
### Paff Distribution ###
paff_distribution <- group_by(paff, systematic) %>%
  summarise(Max = max(paff),
            Median = median(paff),
            Min = min(paff),
            .groups='drop') %>%
  mutate(Range = Max - Min) %>%
  pivot_longer(-systematic, names_to = 'metric', values_to = 'value')

plots$paff_distribution <- ggplot(paff_distribution, aes(x = metric, y = value)) +
  geom_boxplot(shape=20) +
  coord_flip() +
  theme(axis.title.y = element_blank())

plots$paff_extreme_distributions <- (filter(paff, systematic %in% (filter(paff_distribution, metric == 'Range') %>%
                                             filter(value < quantile(value, probs = 0.01)) %>%
                                             pull(systematic))) %>%
  ggplot(aes(x = systematic, y = paff)) +
  geom_boxplot(shape=20) +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dotted', colour = 'grey'))) %>%
  labeled_plot(width = 20, height = 20, units = 'cm')

###  Regress out RNA vs Proteome ###
tidy_lm <- function(x, type){
  s <- summary(x)
  tidy(x) %>%
    mutate(term = str_to_lower(term) %>% str_remove_all('[()]')) %>%
    rename_all(~str_replace_all(., '\\.', '_')) %>%
    pivot_wider(names_from = term, values_from = c(estimate, std_error, statistic, p_value), names_glue = "{term}_{.value}") %>%
    select(sort(tidyselect::peek_vars())) %>%
    bind_cols(tibble(type = type, adj_r_squared = s$adj.r.squared, f_Statistic = s$fstatistic[1],
                     p_value = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)), .) %>%
    return()
}

calc_gene_lms <- function(tbl, ...){
  tbl <- drop_na(tbl)
  
  if (nrow(tbl) < 3){
    return(tibble(type=NA, adj_r_squared=NA, f_Statistic=NA, p_value=NA, intercept_estimate=NA,
                  intercept_p_value=NA, intercept_statistic=NA, intercept_std_error=NA,
                  transcriptomic_estimate=NA, transcriptomic_p_value=NA, transcriptomic_statistic=NA,
                  transcriptomic_std_error=NA))
  }
  
  trans_lm <- lm(proteomic ~ transcriptomic, data = tbl)
  both_lm <- lm(proteomic ~ transcriptomic + paff, data = tbl)
  return(bind_rows(tidy_lm(trans_lm, type = 'Transcriptome'),
                   tidy_lm(both_lm, type = 'Transcriptome + P(aff)')))
}

proteome_paff_lm <- group_by(comb, systematic) %>%
  group_modify(calc_gene_lms) %>%
  filter(!is.na(type))

plots$paff_proteome_lm <- select(proteome_paff_lm, systematic, type, adj_r_squared) %>%
  pivot_wider(names_from = type, values_from = adj_r_squared) %>%
  mutate(diff = `Transcriptome + P(aff)` - Transcriptome,
         big = ifelse(diff > 0.05, '> 0.05', '< 0.05')) %>%
  ggplot(aes(x = Transcriptome, y = `Transcriptome + P(aff)`, colour = big)) + 
  geom_point(shape = 20) +
  guides(colour=guide_legend(title = 'R<sup>2</sup> Difference')) +
  scale_colour_manual(values = c(`> 0.05`='firebrick2', `< 0.05`='cornflowerblue')) +
  labs(subtitle = 'Adj. R-Squared for LMs using different factors to predict Proteome FC') +
  theme(legend.title = element_markdown())

### Save Plots ###
save_plotlist(plots, 'figures/')
