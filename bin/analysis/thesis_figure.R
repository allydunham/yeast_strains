#!/us/bin/env Rscript
# Generate figure proteomics phenotypes for thesis
source("src/config.R")
library(multipanelfigure)

blank_plot <- function(text = ''){
  ggplot(tibble(x=c(0, 1)), aes(x=x, y=x)) +
    geom_blank() +
    annotate(geom = 'text', x = 0.5, y = 0.5, label = text) +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
}

omics <- read_rds('data/rdata/omics.rds') %>%
  filter(!low_paff_range)

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

classify_p <- function(trans, prot){
  out <- rep('Neither', length(trans))
  out[trans < 0.05] <- 'Transcriptomic'
  out[prot < 0.05] <- 'Proteomic'
  out[trans < 0.05 & prot < 0.05] <- 'Both'
  return(factor(out, levels = c('Neither', 'Transcriptomic', 'Proteomic', 'Both')))
}

cor_cats <- select(per_gene_cor, name, systematic, p.adj) %>%
  pivot_wider(names_from = name, values_from = p.adj) %>%
  mutate(cat = classify_p(transcriptomic, proteomic)) %>%
  count(cat)

p_correlation_cat <- ggplot(cor_cats, aes(x = reorder(cat, -n), y = n)) +
  geom_col(width = 0.5, fill = "cornflowerblue", show.legend = FALSE) +
  coord_flip() +
  stat_summary(geom = "text", fun.data = function(x) {data.frame(y = x, label = x)}, hjust = -0.2, size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "", y = "Number of genes where p<sub>adj</sub> < 0.05") +
  theme(axis.title.x = element_markdown(),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

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
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(x = "Transcriptome", y = "Transcriptome + P(aff)") +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-5,0,-10,0))

#### Panel - Linear Models ####

p_linear_models <- blank_plot("Linear Models")

#### Panel - Per Gene Models ####

p_per_gene <- blank_plot("Per Gene")

#### Figure Assembly ####
size <- theme(text = element_text(size = 8))
p1 <- p_correlation + labs(tag = 'A') + size
p2 <- p_correlation_cat + labs(tag = "B") + size
p3 <- p_regression + labs(tag = 'C') + size
p4 <- p_linear_models + labs(tag = 'D') + size
p5 <- p_per_gene + labs(tag = 'E') + size

figure <- multi_panel_figure(width = 180, height = 120, columns = 2, rows = 4,
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1) %>%
  fill_panel(p2, row = 2, column = 1) %>%
  fill_panel(p3, row = 1:2, column = 2) %>%
  fill_panel(p4, row = 3:4, column = 1) %>%
  fill_panel(p5, row = 3:4, column = 2)

ggsave('figures/thesis_figure.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
