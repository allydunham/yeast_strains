#!/usr/bin/env Rscript
# Analyse Paff distribution
source('src/config.R')

### Import Data ###
paff <- read_rds('data/rdata/paff.rds')

plots <- list()
### Analyse ###
paff_distribution <- group_by(paff, systematic) %>%
  summarise(Max = max(paff),
            Median = median(paff),
            Min = min(paff),
            .groups='drop') %>%
  mutate(Range = Max - Min) %>%
  pivot_longer(-systematic, names_to = 'metric', values_to = 'value')

plots$distribution <- ggplot(paff_distribution, aes(x = metric, y = value)) +
  geom_boxplot(shape=20) +
  coord_flip() +
  theme(axis.title.y = element_blank())

plots$extreme_distributions <- (filter(paff, systematic %in% (filter(paff_distribution, metric == 'Range') %>%
                                             filter(value < quantile(value, probs = 0.01)) %>%
                                             pull(systematic))) %>%
  ggplot(aes(x = systematic, y = paff)) +
  geom_boxplot(shape=20) +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dotted', colour = 'grey'))) %>%
  labeled_plot(width = 20, height = 20, units = 'cm')

### Save Plots ###
save_plotlist(plots, 'figures/paff/', overwrite = "all")
