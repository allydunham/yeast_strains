#!/us/bin/env Rscript
# Test phenotype models based on the number of lm associated genes knocked out
source("src/config.R")

filtered_strains <- read_tsv('meta/strain_information.tsv') %>%
  filter(Ploidy <= 2, Aneuploidies == 'euploid') %>%
  pull(`Standardized name`) %>%
  setdiff(. , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))

conditions <- c(
  "2,4-Dichlorophenoxyacetic acid (48H)", "20ºC (72H)", "39ºC (72H)", "42ºC (72H)", "5-FU (48H)",
  "6-AU (48H)", "6-AU + 39ºC (72H)", "aa starvation (48H)", "Acetic acid (48H)",
  "Amphotericin B (48H)", "Amphotericin B + anaerobic (48H)", "Anaerobic growth (48H)", "Cadmium chloride (48H)",
  "Caffeine 15mM (72H)", "Caffeine 20mM (72H)", "Caspofungin (48H)", "Clioquinol (72H)",
  "Clozapine (48H)", "Cyclohexamide (72H)", "DMSO 1%  (48H)", "Glucose 20% (48H)",
  "Glycerol 2%  (72H)", "NaCl 0.4M (72H)", "NaCl 0.4M + 39ºC (72H)", "NaCl 0.6M (72H)",
  "NaCl 0.6M + 39ºC (72H)", "NiSO4 (48H)", "Nitrogen starvation (48H)", "Nystatin (48H)", "Paraquat (48H)",
  "SC (48H)", "SC + hepes (48H)", "Sorbitol 1M (48H)", "YPD (24H)"
)

growth <- read_rds("data/rdata/bede_phenotypes.rds") %>%
  filter(strain %in% filtered_strains, condition %in% conditions) %>%
  filter(subset == 'liti')

paff <- read_rds("data/rdata/paff.rds") %>%
  filter(strain %in% filtered_strains) %>%
  mutate(ko = paff > 0.5)

gene_lms <- read_rds('data/rdata/per_gene_lms.rds') %>%
  drop_na() %>%
  mutate(padj = p.adjust(pvalue, 'fdr')) %>%
  ungroup()

ko_counts <- filter(gene_lms, rsquared > 0.1, padj < 0.01) %>%
  select(condition, systematic) %>%
  left_join(select(paff, systematic, strain, ko), by = "systematic") %>%
  group_by(condition, strain) %>%
  summarise(ko_count = sum(ko), .groups = "drop")

growth_lms <- select(growth, condition, strain, score) %>%
  left_join(ko_counts, by = c("strain", "condition")) %>%
  group_by(condition) %>%
  drop_na() %>%
  group_modify(~broom::glance(lm(score ~ ko_count, data = .)))
