#!/us/bin/env Rscript
# Generate figure proteomics phenotypes for thesis
source("src/config.R")
library(multipanelfigure)

omics <- read_rds('data/rdata/omics.rds') %>%
  filter(!low_paff_range)

genes <- read_tsv("meta/sacc_gene_loci")

#### Panel - Data Description ####
data_summary <- select(omics, systematic, strain, Proteome=proteomic, Transcriptome=transcriptomic, `P(Aff)`=paff) %>%
  pivot_longer(c(-systematic, -strain), names_to = "measure", values_to = "value") %>%
  drop_na() %>%
  group_by(measure) %>%
  summarise(Genes = n_distinct(systematic), Strains = n_distinct(strain)) %>%
  pivot_longer(-measure, names_to = "item", values_to = "count") %>%
  mutate(measure = factor(measure, levels = c("Proteome", "Transcriptome", "P(Aff)")))

p_data <- ggplot(data_summary, aes(x = measure, y = count)) +
  facet_wrap(~item, scales = "free_x", strip.position = "bottom") +
  geom_col(fill = "cornflowerblue", width = 0.55) +
  stat_summary(geom = "text", fun.data = function(i) {data.frame(y = i, label = i)}, hjust = -0.2, size = 2.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        panel.spacing = unit(3, "mm"),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        strip.placement = "outside")

#### Panel - Distributions ####
omics_long <- select(omics, systematic, strain, `"Proteomic log"[2]~"FC"`=proteomic,
                     `"Transcriptomic log"[2]~"FC"`=transcriptomic, `P(Aff)`=paff) %>%
  pivot_longer(c(-systematic, -strain), names_to = "measure", values_to = "value")

p_data_dist <- ggplot(omics_long, aes(x = value, y = ..scaled.., colour = measure)) +
  facet_wrap(~measure, scales = "free", strip.position = "bottom", labeller = labeller(measure=label_parsed)) +
  stat_density(geom = "line", show.legend = FALSE) +
  labs(y = "Scaled\nDensity") +
  scale_x_continuous(name = "", labels = function(x) str_remove(x, "\\.?0*^")) +
  theme(strip.placement = "outside")

#### Panel - Correlation ####
group_cor_test <- function(tbl, var1, var2){
  var1 <- enquo(var1)
  var2 <- enquo(var2)
  tbl <- drop_na(tbl, !!var1, !!var2)
  
  if (nrow(tbl) < 3){
    return(tibble(estimate=NA, statistic=NA, p.value=NA, parameter=NA, method=NA, alternative=NA))
  }
  return(tidy(cor.test(pull(tbl, !!var1), pull(tbl, !!var2))))
}

per_gene_cor <- select(omics, systematic, strain, proteomic, transcriptomic, paff) %>%
  pivot_longer(c(proteomic, transcriptomic), names_to = "name", values_to = "value") %>%
  group_by(name, systematic) %>%
  group_modify(~group_cor_test(.x, value, paff)) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(p.adj = p.adjust(p.value))

p_correlation <- ggplot(per_gene_cor, aes(x = estimate, colour = str_to_sentence(name))) +
  stat_density(geom = 'line', position = 'identity') +
  labs(x = 'Correlation Coefficient', y = 'Density') +
  scale_colour_brewer(name = '', type = 'qual', palette = 'Dark2', direction = -1) +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,0,-10,0))

#### Panel - Regression ####
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

proteome_lm <- group_by(omics, systematic) %>%
  group_modify(calc_gene_lms) %>%
  ungroup() %>%
  filter(!is.na(type)) %>%
  select(systematic, type, r) %>%
  pivot_wider(names_from = type, values_from = r) %>%
  mutate(diff = `Transcriptome + P(aff)` - Transcriptome,
         big = diff > 0.05,
         big_str = ifelse(big,
                          str_c('> 0.05 (n = ', sum(big, na.rm = TRUE), ')'),
                          str_c('≤ 0.05 (n = ', sum(!big, na.rm = TRUE), ')')))

p_regression <- ggplot(proteome_lm, aes(x = Transcriptome, y = `Transcriptome + P(aff)`, colour = big_str)) + 
  geom_point(shape = 20, size = 0.5) +
  scale_colour_manual(name = "R difference", values = c(`> 0.05 (n = 223)`='firebrick2', `≤ 0.05 (n = 838)`='cornflowerblue')) +
  guides(colour = guide_legend(override.aes = list(size = 1.8))) +
  labs(x = "Transcriptome", y = "Transcriptome + P(aff)") +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,0,-10,0),
        legend.key.size = unit(1.8, "mm"))

#### Panel - Linear Models ####
phenotype_lms <- read_rds('data/rdata/phenotype_lms.rds')

nfactor_map <- c(`Clade`=1, `Genetic Distance`=1, `P(Aff)`=1, `Transcriptomic`=1, `Proteomic`=1, `Transcriptomic/P(Aff)`=2,
                 `Proteomic/P(Aff)`=2, `Proteomic/Transcriptomic`=2, `All`=3,
                 `VAE - All`=3, `VAE - Omics`=2)
p_linear_models <- mutate(phenotype_lms, nfactors = as.character(nfactor_map[as.character(type)])) %>%
  ggplot(aes(x = type, y = adj_r_squared, fill = nfactors)) +
  geom_boxplot(outlier.shape = 20, outlier.size = 0.5) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  ylim(c(0, 0.6)) +
  scale_fill_brewer(name = "Factors", type = 'qual', palette = 'Set1') +
  labs(x = '', y = 'Adj. R Squared') +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,0,-10,0),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

#### Panel - Per Gene Model Summary ####
gene_lms <- read_rds('data/rdata/per_gene_lms.rds') %>%
  drop_na() %>%
  mutate(padj = p.adjust(pvalue, 'fdr')) %>%
  ungroup()

p_value_cutoff <- 0.01
factor_summary <- mutate(gene_lms,
                         trans = p.adjust(p_value_transcriptomic, 'fdr') < p_value_cutoff,
                         prot = p.adjust(p_value_proteomic, 'fdr') < p_value_cutoff,
                         paff = p.adjust(p_value_paff, 'fdr') < p_value_cutoff) %>%
  summarise(Transcriptomic = sum(trans & !prot & !paff),
            Proteomic = sum(prot & !trans & !paff),
            `P(Aff)` = sum(paff & !prot & !trans),
            `Transcriptomic & Proteomic` = sum(trans & prot & !paff),
            `Transcriptomic & P(Aff)` = sum(trans & paff & !prot),
            `Proteomic & P(Aff)` = sum(prot & paff & !trans),
            All = sum(trans & prot & paff), .groups = 'drop') %>%
  pivot_longer(everything(), names_to = 'factors', values_to = 'count') %>%
  mutate(n_factors = c(None = 0, Transcriptomic = 1, Proteomic = 1, `P(Aff)` = 1, `Transcriptomic & Proteomic` = 2,
                       `Transcriptomic & P(Aff)` = 2, `Proteomic & P(Aff)` = 2, All = 3)[factors],
         factors = factor(factors, levels = c('All', 'Proteomic & P(Aff)', 'Transcriptomic & P(Aff)', 'Transcriptomic & Proteomic',
                                              'P(Aff)', 'Proteomic', 'Transcriptomic', 'None')))

p_per_gene_summary <- ggplot(factor_summary, aes(x = factors, y = count, fill = as.character(n_factors))) +
  geom_col() +
  stat_summary(geom = "text", fun.data = function(i) {data.frame(y = i, label = i)}, hjust = -0.2, size = 2.5) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  scale_fill_brewer(name = "Factors", type = "qual", palette = "Set1") +
  labs(y = "Associations", x = "") +
  theme(legend.position = "top",
        legend.title = element_text(vjust = 1),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,-5,0,-10),
        legend.key.size = unit(4, units = "mm"),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

#### Panel - Per Gene Model Examples ####
examples_genes <- filter(gene_lms, condition %in% c("NaCl 0.6M (48H)", "42ºC (48H)", "Caffeine 20mM (48H)")) %>%
  left_join(select(genes, systematic = id, name), by = "systematic") %>%
  mutate(condition = str_remove(condition, " \\(48H\\)"))

p_per_gene_examples <- ggplot(examples_genes, aes(x = rsquared, y = padj)) +
  facet_wrap(~condition, nrow = 1) +
  geom_point(aes(colour = padj < 0.05), shape = 20, size = 0.75) +
  geom_text_repel(data = filter(examples_genes, rsquared > 0.25), mapping = aes(label = name),
                  ylim = c(0, Inf), force = 40, size = 2.5) +
  labs(y = expression(p[adj])) +
  scale_x_continuous(name = expression(R^2), labels = function(x) str_remove(x, "\\.?0*^")) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_colour_manual(name = 'p<sub>adj</sub>',
                      values = c(`TRUE`='red', `FALSE`='black'),
                      labels = c(`TRUE`='< 0.05', `FALSE`='> 0.05')) + 
  theme(legend.title = element_markdown(),
        legend.position = "bottom",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,0,-5,0),
        axis.ticks.y = element_blank())

#### Figure - Data Properties ####
size <- theme(text = element_text(size = 11))
p1_prop <- p_data + labs(tag = 'A') + size
p2_prop <- p_data_dist + labs(tag = 'B') + size
p3_prop <- p_correlation + labs(tag = 'C') + size
p4_prop <- p_regression + labs(tag = 'D') + size
  
figure_prop <- multi_panel_figure(width = c(90,90), height = c(40, 40, 60),
                                  panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1_prop, row = 1, column = 1:2) %>%
  fill_panel(p2_prop, row = 2, column = 1:2) %>%
  fill_panel(p3_prop, row = 3, column = 1) %>%
  fill_panel(p4_prop, row = 3, column = 2)

#### Figure - Models ####
size <- theme(text = element_text(size = 11))
p1_mod <- p_linear_models + guides(fill = "none") + labs(tag = 'A') + size
p2_mod <- p_per_gene_summary + labs(tag = 'B') + size
p3_mod <- p_per_gene_examples + labs(tag = "C") + size

figure_mod <- multi_panel_figure(width = 180, height = 140, rows = 2, columns = 2,
                                 panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1_mod, row = 1:2, column = 1) %>%
  fill_panel(p2_mod, row = 1, column = 2) %>%
  fill_panel(p3_mod, row = 2, column = 2)

### Save ###
ggsave('figures/thesis_figure_properties.pdf', figure_prop, width = figure_width(figure_prop), height = figure_height(figure_prop), units = 'mm',
       device = cairo_pdf)
ggsave('figures/thesis_figure_properties.png', figure_prop, width = figure_width(figure_prop), height = figure_height(figure_prop), units = 'mm')

ggsave('figures/thesis_figure_model.pdf', figure_mod, width = figure_width(figure_mod), height = figure_height(figure_mod), units = 'mm',
       device = cairo_pdf)
ggsave('figures/thesis_figure_model.png', figure_mod, width = figure_width(figure_mod), height = figure_height(figure_mod), units = 'mm')
