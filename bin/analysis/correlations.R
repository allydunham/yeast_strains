#!/usr/bin/env Rscript
# Analyse Basic Correlations
source('src/config.R')

### Import Data ###
omics <- read_rds('data/rdata/omics.rds') %>%
  filter(!low_paff_range)

### Analyse ###
plots <- list()
plots$proteome_vs_transcriptome <- ggplot(omics, aes(x = transcriptomic, y = proteomic)) + 
  geom_hex()

plots$proteome_vs_paff <- ggplot(omics, aes(x = paff, y = proteomic, group = cut_width(paff, 0.05))) + 
  geom_violin()

plots$transcriptome_vs_paff <- ggplot(omics, aes(x = paff, y = transcriptomic, group = cut_width(paff, 0.05))) + 
  geom_violin()

### Per gene Correlations ###
group_cor_test <- function(tbl, var1, var2){
  var1 <- enquo(var1)
  var2 <- enquo(var2)
  tbl <- drop_na(tbl, !!var1, !!var2)
  
  if (nrow(tbl) < 3){
    return(tibble(estimate=NA, statistic=NA, p.value=NA, parameter=NA, method=NA, alternative=NA))
  }
  return(tidy(cor.test(pull(tbl, !!var1), pull(tbl, !!var2))))
}

per_gene_cross <- group_by(omics, systematic) %>%
  group_modify(~group_cor_test(.x, proteomic, transcriptomic)) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(p.adj = p.adjust(p.value))

per_gene_transcriptomic <- group_by(omics, systematic) %>%
  group_modify(~group_cor_test(.x, transcriptomic, paff)) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(p.adj = p.adjust(p.value))

per_gene_proteomic <- group_by(omics, systematic) %>%
  group_modify(~group_cor_test(.x, proteomic, paff)) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(p.adj = p.adjust(p.value))

gene_prot_trans <- full_join(select(per_gene_cross, systematic, p_prot_trans=p.adj),
                             select(per_gene_transcriptomic, systematic, p_trans_paff=p.adj),
                             by = c('systematic')) %>%
  full_join(select(per_gene_proteomic, systematic, p_prot_paff=p.adj),
            by = c('systematic'))

plots$gene_cors <- ggplot() +
  geom_point(data = filter(gene_prot_trans, p_trans_paff < 0.05, p_prot_paff < 0.05),
             mapping = aes(x = p_trans_paff, y = p_prot_paff), colour = 'red') + 
  geom_point(data = filter(gene_prot_trans, p_trans_paff < 0.05, p_prot_paff > 0.05),
             mapping = aes(x = p_trans_paff, y = p_prot_paff), colour = 'green', shape = 20) + 
  geom_point(data = filter(gene_prot_trans, p_trans_paff > 0.05, p_prot_paff < 0.05),
             mapping = aes(x = p_trans_paff, y = p_prot_paff), colour = 'green', shape = 20) + 
  geom_point(data = filter(gene_prot_trans, p_trans_paff > 0.05, p_prot_paff > 0.05),
             mapping = aes(x = p_trans_paff, y = p_prot_paff), colour = 'black', shape = 20) + 
  coord_cartesian(clip = 'off')

classify_p <- function(trans, prot){
  out <- rep('Neither', length(trans))
  out[trans < 0.05] <- 'Transcriptomic'
  out[prot < 0.05] <- 'Proteomic'
  out[trans < 0.05 & prot < 0.05] <- 'Both'
  return(factor(out, levels = c('Neither', 'Transcriptomic', 'Proteomic', 'Both')))
}

plots$gene_cors_bars <- mutate(gene_prot_trans, type = classify_p(p_trans_paff, p_prot_paff)) %>%
  ggplot(aes(x = type)) +
  geom_bar() + 
  scale_y_log10()

### Direction of Correlation
cors <- bind_rows(Proteomic = per_gene_proteomic, Transcriptomic = per_gene_transcriptomic, .id = 'omic') %>%
  mutate(cat = ifelse(p.adj <= 0.05, 'p<sub>adj</sub> &le; 0.05', 'p<sub>adj</sub> &gt; 0.05'))

plots$cor_volcano <- ggplot(cors, aes(x = estimate, y = -log10(p.value), colour = cat)) +
  facet_wrap(~omic) +
  geom_point() +
  labs(x = 'Correlation Coefficient', y = '-log<sub>10</sub> p') +
  scale_colour_brewer(name = '', type = 'qual', palette = 'Set1', direction = -1) +
  theme(legend.text = element_markdown(), axis.title.y = element_markdown())

plots$cor_distribuition <- ggplot(cors, aes(x = estimate, colour = omic)) +
  stat_density(geom = 'line', position = 'identity') +
  labs(x = 'Correlation Coefficient', y = 'Density') +
  scale_colour_brewer(name = '', type = 'qual', palette = 'Dark2', direction = -1)

### Save Plots ###
save_plotlist(plots, 'figures/correlation/')
