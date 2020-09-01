#!/usr/bin/env Rscript
# Analyse Bede's S-Score Phenotype Data
source('src/config.R')
library(caret)

### Import Data ###
phenotypes <- read_rds('data/rdata/bede_phenotypes.rds') %>%
  drop_na()
omics <- read_rds('data/rdata/omics.rds') %>%
  filter(!low_paff_range)
omic_pcas <- read_rds('data/rdata/omic_pcas.rds') %>%
  select(strain, condition, score, qvalue, num_range('proteomic_PC', range = 1:50),
         num_range('transcriptomic_PC', range = 1:50), num_range('paff_PC', range = 1:50))

vae <- bind_rows(`VAE - All` = read_tsv('data/vae/profiles.tsv'),
                 `VAE - Omics` = read_tsv('data/vae_omics/profiles.tsv'),
                 .id = 'model')

# Identify conditions with reasonable number of negative examples
good_cons <- group_by(phenotypes, condition) %>%
  drop_na(score, qvalue) %>%
  summarise(n = n(),
            n_del = sum(score < 0 & qvalue < 0.05),
            n_neut = sum(qvalue >= 0.05),
            n_pos = sum(score > 0 & qvalue < 0.05),
            .groups = 'drop') %>%
  filter(n_del > 25)

plots <- list()
### S-Score distributions and selection of conditions with good ranges ###
plots$sscores <- (
  mutate(phenotypes, sig = ifelse(qvalue < 0.05, 'q < 0.05', 'q > 0.05')) %>%
  ggplot(aes(x = score, fill = sig)) +
    facet_wrap(~condition, scales = 'free_y') +
    geom_histogram(bins = 30) +
    labs(x = 'S-Score', y = 'Count') +
    scale_fill_brewer(type='qual', palette = 'Set1', direction = 1) +
    guides(fill = guide_legend(title = ''))
) %>% labeled_plot(units = 'cm', width = 40, height = 30)


### Generate Linear Models ###
tidy_lm <- function(x, type){
  s <- summary(x)
  return(tibble(type = type, adj_r_squared = s$adj.r.squared, f_Statistic = s$fstatistic[1],
                p_value = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)))
}

vae_lms <- select(phenotypes, strain, condition, score) %>%
  left_join(vae, by = 'strain') %>%
  drop_na(model) %>%
  select(-strain) %>%
  group_by(condition, model) %>%
  group_modify(~tidy_lm(lm(score ~ ., data = .), type = 'VAE')) %>%
  select(-type) %>%
  rename(type = model)

calc_lms <- function(tbl, ...){
  tbl <- drop_na(tbl)
  
  if (nrow(tbl) < 3){
    return(tibble(type=NA, adj_r_squared=NA, f_Statistic=NA, p_value=NA))
  }
  
  prot_lm <- lm(score ~ ., data = select(tbl, score, starts_with('proteomic')))
  trans_lm <- lm(score ~ ., data = select(tbl, score, starts_with('transcriptomic')))
  paff_lm <- lm(score ~ ., data = select(tbl, score, starts_with('paff')))
  
  prot_trans_lm <- lm(score ~ ., data = select(tbl, score, starts_with('proteomic'), starts_with('transcriptomic')))
  prot_paff_lm <- lm(score ~ ., data = select(tbl, score, starts_with('proteomic'), starts_with('paff')))
  trans_paff_lm <- lm(score ~ ., data = select(tbl, score, starts_with('transcriptomic'), starts_with('paff')))
  
  comb_lm <- lm(score ~ ., data = select(tbl, -strain))
  
  return(
    bind_rows(
      tidy_lm(prot_lm, type = 'Proteomic'),
      tidy_lm(trans_lm, type = 'Transcriptomic'),
      tidy_lm(paff_lm, type = 'P(Aff)'),
      tidy_lm(prot_trans_lm, type = 'Proteomic/Transcriptomic'),
      tidy_lm(prot_paff_lm, type = 'Proteomic/P(Aff)'),
      tidy_lm(trans_paff_lm, type = 'Transcriptomic/P(Aff)'),
      tidy_lm(comb_lm, type = 'All')
    )
  )
}


calc_logistic_models <- function(tbl, ...){
  tbl <- drop_na(tbl) %>%
    mutate(del = as.factor(ifelse(score < 0 & qvalue < 0.05, 'deleterious', 'neutral'))) %>%
    select(del, starts_with('proteomic_'), starts_with('transcriptomic_'), starts_with('paff_'))
  
  if (nrow(tbl) < 3){
    return(tibble(type=NA, accuracy=NA, kappa=NA))
  }
  
  prot_lm <- train(del ~ ., data = select(tbl, del, starts_with('proteomic')), method = "glm", family=binomial())
  trans_lm <- train(del ~ ., data = select(tbl, del, starts_with('transcriptomic')), method = "glm", family=binomial())
  paff_lm <- train(del ~ ., data = select(tbl, del, starts_with('paff')), method = "glm", family=binomial())
  
  prot_trans_lm <- train(del ~ ., data = select(tbl, del, starts_with('proteomic'), starts_with('transcriptomic')), method = "glm", family=binomial())
  prot_paff_lm <- train(del ~ ., data = select(tbl, del, starts_with('proteomic'), starts_with('paff')), method = "glm", family=binomial())
  trans_paff_lm <- train(del ~ ., data = select(tbl, del, starts_with('transcriptomic'), starts_with('paff')), method = "glm", family=binomial())
  
  comb_lm <- train(del ~ ., data = tbl, method = "glm", family=binomial())
  
  return(
    bind_rows(
      tibble(type='Proteomic', accuracy=prot_lm$results$Accuracy, kappa=prot_lm$results$Kappa),
      tibble(type='Transcriptomic', accuracy=trans_lm$results$Accuracy, kappa=trans_lm$results$Kappa),
      tibble(type='P(Aff)', accuracy=paff_lm$results$Accuracy, kappa=paff_lm$results$Kappa),
      tibble(type='Proteomic/Transcriptomic', accuracy=prot_trans_lm$results$Accuracy, kappa=prot_trans_lm$results$Kappa),
      tibble(type='Proteomic/P(Aff)', accuracy=prot_paff_lm$results$Accuracy, kappa=prot_paff_lm$results$Kappa),
      tibble(type='Transcriptomic/P(Aff)', accuracy=trans_paff_lm$results$Accuracy, kappa=trans_paff_lm$results$Kappa),
      tibble(type='All', accuracy=comb_lm$results$Accuracy, kappa=comb_lm$results$Kappa)
    )
  )
}

phenotype_lms <- select(omic_pcas, -qvalue) %>%
  group_by(condition) %>%
  group_modify(calc_lms) %>%
  drop_na() %>%
  bind_rows(vae_lms) %>%
  mutate(type = factor(type, levels = c('P(Aff)', 'Transcriptomic', 'Proteomic', 'Transcriptomic/P(Aff)',
                                        'Proteomic/P(Aff)', 'Proteomic/Transcriptomic', 'All', 'VAE - Omics', 'VAE - All')))

# phenotype_logisitics <- group_by(omic_pcas, condition) %>%
#   group_modify(calc_logistic_models) %>%
#   drop_na() %>%
#   mutate(type = factor(type, levels = c('P(Aff)', 'Transcriptomic', 'Proteomic', 'Transcriptomic/P(Aff)',
#                                         'Proteomic/P(Aff)', 'Proteomic/Transcriptomic', 'All')))
phenotype_logisitics <- read_rds('data/rdata/phenotype_logistics.rds')

### Analyse  Models ###
lm_colours <- c(`P(Aff)`='yellow', `Transcriptomic`='magenta', `Proteomic`='cyan',
                `Transcriptomic/P(Aff)`='red', `Proteomic/P(Aff)`='green', `Proteomic/Transcriptomic`='blue',
                `All` = 'black', 'VAE - All'='orange', 'VAE - Omics'='orange')
plots$lm_factor_r_squared <- ggplot(phenotype_lms, aes(x = type, y = adj_r_squared, fill = type)) +
  facet_wrap(~condition, nrow = 5, ncol = 10) +
  geom_col(show.legend = FALSE) +
  coord_flip() + 
  scale_fill_manual(values = lm_colours) +
  labs(x = 'Adj. R Squared', y = '') +
  theme(panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted'),
        panel.grid.major.y = element_blank())
plots$lm_factor_r_squared <- labeled_plot(plots$lm_factor_r_squared, units = 'cm', height = 30, width = 50)

nfactor_map <- c(`P(Aff)`=1, `Transcriptomic`=1, `Proteomic`=1, `Transcriptomic/P(Aff)`=2,
                 `Proteomic/P(Aff)`=2, `Proteomic/Transcriptomic`=2, `All`=3,
                 `VAE - All`=3, `VAE - Omics`=2)
plots$lm_factor_r_squared_summary <- mutate(phenotype_lms, nfactors = as.character(nfactor_map[as.character(type)])) %>%
  ggplot(aes(x = type, y = adj_r_squared, fill = nfactors)) +
  geom_boxplot() +
  ylim(c(0, NA)) +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  guides(fill = guide_legend(title = 'Number of Factors')) +
  labs(x = '', y = 'Adj. R Squared')
plots$lm_factor_r_squared_summary <- labeled_plot(plots$lm_factor_r_squared_summary, units = 'cm', height = 15, width = 40)

plots$log_factors_kappa <- ggplot(phenotype_logisitics, aes(x = type, y = kappa, fill = type)) +
  facet_wrap(~condition, ncol = 6) +
  geom_col(show.legend = FALSE) +
  coord_flip() + 
  scale_fill_manual(values = lm_colours) +
  labs(x = "Cohen's Kappa", y = '') +
  theme(panel.grid.major.x = element_line(colour = 'grey', linetype = 'dotted'),
        panel.grid.major.y = element_blank())
plots$log_factors_kappa <- labeled_plot(plots$log_factors_kappa, units = 'cm', height = 25, width = 35)

plots$log_factors_kappa_summary <- mutate(phenotype_logisitics, nfactors = as.character(nfactor_map[type])) %>%
  ggplot(aes(x = type, y = kappa, fill = nfactors)) +
  geom_boxplot() +
  ylim(c(0, NA)) +
  scale_fill_brewer(type = 'qual', palette = 'Set1') +
  guides(fill = guide_legend(title = 'Number of Factors')) +
  labs(x = '', y = "Cohen's Kappa")
plots$log_factors_kappa_summary <- labeled_plot(plots$log_factors_kappa_summary, units = 'cm', height = 15, width = 35)

plots$log_kappa_accuracy <- ggplot() +
  geom_point(data = phenotype_logisitics, mapping = aes(x = accuracy, y = kappa, colour = type)) +
  geom_vline(data = mutate(good_cons, base_acc = n_neut / n), mapping = aes(xintercept = base_acc), linetype='dotted') +
  facet_wrap(~condition, ncol = 3) +
  scale_colour_manual(values = lm_colours) +
  guides(colour = guide_legend(title = '')) +
  labs(x = 'Accuracy', y = "Cohen's Kappa")
plots$log_kappa_accuracy <- labeled_plot(plots$log_kappa_accuracy, units = 'cm', height = 30, width = 30)

### Cyclohexamide model ###
cyclo <- filter(omic_pcas, condition == 'Cyclohexamide (72H)') %>%
  mutate(del = as.factor(ifelse(score < 0 & qvalue < 0.05, 'deleterious', 'neutral'))) %>%
  select(del, starts_with('proteomic_'), starts_with('transcriptomic_'), starts_with('paff_')) %>%
  drop_na()

trainInd <- createDataPartition(cyclo$del, times = 1, p = 0.8, list = FALSE)
cyclo_train <- cyclo[trainInd[,1],]
cyclo_test <- cyclo[-trainInd[,1],]

tc <- trainControl(method = "repeatedcv", repeats = 3, savePredictions = TRUE, classProbs = TRUE, summaryFunction = twoClassSummary)
cyclo_paff <- train(del ~ ., data = select(cyclo_train, del, starts_with('paff')), method = "glm", family=binomial(), trControl=tc)

preds <- predict(cyclo_paff, cyclo_test)
cm <- confusionMatrix(preds, cyclo_test$del)

### Save Plots ###
save_plotlist(plots, 'figures/bede_phenotypes/')
