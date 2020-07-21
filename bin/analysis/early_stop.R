#!/usr/bin/env Rscript
# Analyse Early Stops
source('src/config.R')

### Import Data ###
early_stops <- read_rds('data/rdata/early_stops.rds') %>%
  filter(!low_paff_range)
omics <- read_rds('data/rdata/omics.rds') %>%
  filter(!low_paff_range)

### Analyse ###
plots <- list()

stops <- group_by(early_stops, strain, systematic) %>%
  filter(protein_pos <= min(protein_pos)) %>%
  ungroup() %>%
  select(strain, systematic, consequence, prop, pos=protein_pos) %>%
  left_join(omics, ., by = c('systematic', 'strain')) %>%
  mutate(consequence = ifelse(is.na(consequence), 'Missense', str_to_title(consequence)))

# Against Prop
plots$transcriptomic_fc_prop <- (ggplot(drop_na(stops, prop), aes(x = prop, group = cut_width(prop, 0.1), y = transcriptomic)) +
                                facet_wrap(~consequence, nrow = 2, labeller = as_labeller(str_to_title)) +
                                geom_boxplot() +
                                labs(x = 'Proportion through protein', y = 'Transcriptomic FC')) %>% 
  labeled_plot(unit = 'cm', height = 20, width = 20)

plots$proteomic_fc_prop <- (ggplot(drop_na(stops, prop), aes(x = prop, group = cut_width(prop, 0.1), y = proteomic)) +
                               facet_wrap(~consequence, nrow = 2, labeller = as_labeller(str_to_title)) +
                               geom_boxplot() +
                               labs(x = 'Proportion through protein', y = 'Proteomic FC')) %>% 
  labeled_plot(unit = 'cm', height = 20, width = 20)

# Against Absolute Position
classify_position <- function(x){
  c <- rep('> 500', length(x))
  c[x < 500] <- '< 500'
  c[x < 250] <- '< 250'
  c[x < 100] <- '< 100'
  c[x < 50] <- '< 50'
  c[x < 25] <- '< 25'
  c <- factor(c, levels = c('< 25', '< 50', '< 100', '< 250', '< 500', '> 500'))
  return(c)
}

plots$trancriptomic_fc_pos <- (drop_na(stops, pos) %>%
  mutate(pos_class = classify_position(pos)) %>%
  ggplot(aes(x = pos_class, y = transcriptomic)) +
   facet_wrap(~consequence, nrow = 2, labeller = as_labeller(str_to_title)) +
   geom_boxplot() +
   labs(x = 'Protein Position', y = 'Transcriptomic FC')) %>% 
  labeled_plot(unit = 'cm', height = 20, width = 20)

plots$proteomic_fc_pos <- (drop_na(stops, pos) %>%
  mutate(pos_class = classify_position(pos)) %>%
  ggplot(aes(x = pos_class, y = proteomic)) +
  facet_wrap(~consequence, nrow = 2, labeller = as_labeller(str_to_title)) +
  geom_boxplot() +
  labs(x = 'Protein Position', y = 'Proteomic FC')) %>% 
  labeled_plot(unit = 'cm', height = 20, width = 20)

# Overall distribution
plots$overall <- ggplot(stops, aes(x = consequence, y = proteomic)) +
  geom_boxplot() +
  labs(x = '', y = 'Proteomic FC') +
  stat_compare_means(comparisons = list(c('Frameshift', 'Missense'), c('Missense', 'Nonsense'), c('Frameshift', 'Nonsense')))

### Save Plots ###
save_plotlist(plots, 'figures/early_stops')
