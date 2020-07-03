#!/usr/bin/env Rscript
# Initial analysis script
source('src/config.R')

### Import Data ###
proteomic <- read_csv('data/raw/1k_quant_wide_systematic_name.csv', ) %>%
  select(-X1) %>%
  rename(gene=symbol, systematic=systematic_name) %>%
  pivot_longer(!one_of(c('gene', 'systematic')), names_to = 'strain', values_to = 'abundance')

transcriptomic <- read_csv('data/raw/tpm_FinalSet_969strains.csv') %>%
  rename(systematic = systematic_name) %>%
  pivot_longer(-systematic, names_to = 'strain', values_to = 'abundance')

paff <- read_tsv('data/paff_scores.tsv') %>%
  rename(systematic=geneid)
# TODO Strains CEN.PK and Reference? called XTRA_DXL and FY4-6 in their data

comb <- full_join(select(proteomic, -gene) %>% rename(proteomic=abundance),
                  rename(transcriptomic, transcriptomic=abundance),
                  by = c('systematic', 'strain')) %>%
  left_join(paff, by = c('systematic', 'strain'))

### Analyse ###
plots <- list()

# Basic Correlations
plots$proteome_vs_transcriptome <- ggplot(comb, aes(x = transcriptomic, y = proteomic)) + 
  scale_x_log10() +
  scale_y_log10() +
  geom_hex()

plots$proteome_vs_paff <- ggplot(comb, aes(x = paff, y = proteomic, group = cut_width(paff, 0.05))) + 
  scale_y_log10() +
  geom_violin()

plots$transcriptome_vs_paff <- ggplot(comb, aes(x = paff, y = transcriptomic, group = cut_width(paff, 0.05))) + 
  scale_y_log10() +
  geom_violin()

# Per gene
group_cor_test <- function(tbl, var1, var2){
  var1 <- enquo(var1)
  var2 <- enquo(var2)
  tbl <- drop_na(tbl, !!var1, !!var2)
  
  if (nrow(tbl) < 3){
    return(tibble(estimate=NA, statistic=NA, p.value=NA, parameter=NA, method=NA, alternative=NA))
  }
  return(tidy(cor.test(pull(tbl, !!var1), pull(tbl, !!var2))))
}

per_gene_cross <- group_by(comb, systematic) %>%
  group_modify(~group_cor_test(.x, proteomic, transcriptomic)) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(p.adj = p.adjust(p.value))

per_gene_transcriptomic <- group_by(comb, systematic) %>%
  group_modify(~group_cor_test(.x, transcriptomic, paff)) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(p.adj = p.adjust(p.value))

per_gene_proteomic <- group_by(comb, systematic) %>%
  group_modify(~group_cor_test(.x, proteomic, paff)) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(p.adj = p.adjust(p.value))
