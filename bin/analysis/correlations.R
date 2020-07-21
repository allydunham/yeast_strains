#!/usr/bin/env Rscript
# Analyse Basic Correlations
source('src/config.R')

### Import Data ###
paff <- read_tsv('data/paff_scores.tsv') %>%
  rename(systematic=geneid)
# TODO Strains CEN.PK and Reference? called XTRA_DXL and FY4-6 in their data

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

### Analyse ###
plots <- list()
plots$proteome_vs_transcriptome <- ggplot(comb, aes(x = transcriptomic, y = proteomic)) + 
  geom_hex()

plots$proteome_vs_paff <- ggplot(comb, aes(x = paff, y = proteomic, group = cut_width(paff, 0.05))) + 
  geom_violin()

plots$transcriptome_vs_paff <- ggplot(comb, aes(x = paff, y = transcriptomic, group = cut_width(paff, 0.05))) + 
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

### Save Plots ###
save_plotlist(plots, 'figures/')
