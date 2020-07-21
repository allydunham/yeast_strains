#!/usr/bin/env Rscript
# Analyse Early Stops
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

### Analyse ###
plots <- list()

stops <- group_by(early_stops, strain, systematic) %>%
  filter(protein_pos <= min(protein_pos)) %>%
  ungroup() %>%
  select(strain, systematic, consequence, prop) %>%
  left_join(comb, ., by = c('systematic', 'strain')) %>%
  mutate(consequence = ifelse(is.na(consequence), 'Missense', str_to_title(consequence)))

plots$early_stop_trans_fc <- (ggplot(drop_na(stops, prop), aes(x = prop, group = cut_width(prop, 0.1), y = transcriptomic)) +
                                facet_wrap(~consequence, nrow = 2, labeller = as_labeller(str_to_title)) +
                                geom_boxplot() +
                                labs(x = 'Proportion through protein', y = 'Proteomic FC')) %>% 
  labeled_plot(unit = 'cm', height = 20, width = 20)

plots$early_stop_prot_fc <- (ggplot(drop_na(stops, prop), aes(x = prop, group = cut_width(prop, 0.1), y = proteomic)) +
                               facet_wrap(~consequence, nrow = 2, labeller = as_labeller(str_to_title)) +
                               geom_boxplot() +
                               labs(x = 'Proportion through protein', y = 'Proteomic FC')) %>% 
  labeled_plot(unit = 'cm', height = 20, width = 20)

plots$early_stop_prot_overall <- ggplot(stops, aes(x = consequence, y = proteomic)) +
  geom_boxplot() +
  labs(x = '', y = 'Proteomic FC') +
  stat_compare_means(comparisons = list(c('Frameshift', 'Missense'), c('Missense', 'Nonsense'), c('Frameshift', 'Nonsense')))

### Save Plots ###
save_plotlist(plots, 'figures/')
