#!/usr/bin/env Rscript
# Assign mutfunc data to variants from two strains
source("src/config.R")
mutfunc_root <- "~/phd/mutfunc/"

# IDs
mutfunc_ids <- read_tsv("meta/mutfunc_ids.tsv") %>% select(-type)
uniprot_orf <- read_tsv(str_c(mutfunc_root, "id_mapping/uniprot2orf.tab")) %>% rename(uniprot = alias)

# SIFT
sift <- read_tsv(str_c(mutfunc_root, "conservation/sift_parsed_all.tab"), comment = "#") %>%
  select(uniprot = acc, pos_aa = pos, ref_aa = ref, alt_aa = alt, sift_score = score, sift_median_ic = median_ic)

# FoldX
foldx <- read_tsv(str_c(mutfunc_root, "structure/exp_ddg1.tab")) %>%
  select(uniprot = uniprot_id, pos_aa = uniprot_pos, ref_aa = aa_wt, alt_aa = aa_mt, foldx_ddg = ddG)

# Int
ints <- read_tsv("structure/interfaces_fixed.tab")) %>%
  select(uniprot = PROTEIN, pos_aa = POS_UNIPROT, ref_aa = AA, int_partner = PARTNER, int_acc_diff_percent = ACCDIFFPERCENT)

# ELM
elm <- read_tsv(str_c(mutfunc_root, "elm/elm_mut.tsv")) %>%
  select(gene = orf, pos_aa = pos, ref_aa = wt, alt_aa = mt, elm, elm_wt = seq_wt,
         elm_alt = seq_mt, elm_lost = lost, elm_start = start, elm_end = end)

# PTM
dbptm <- as_tibble(readRDS(str_c(mutfunc_root, "ptm/dbptm_yeast_parsed.rds"))) %>% 
    select(gene = orf, pos_aa = position, ref_aa = modified_residue, ptm = modification, ptm_seq = sequence,
           ptm_type = residue_type, ptm_pubmed = pubmed)

phosphogrid <- as_tibble(readRDS(str_c(mutfunc_root, "ptm/Phosphogrid_df_parsed.rds"))) %>%
  select(gene = Orf, pos_aa = Phosphosite, ref_aa = Residue, ptm_seq = seq, ptm_kinase = Kinase_orfs, ptm_pubmed = Site_evidence_pubmed) %>%
  mutate(ptm = "phosphorylation")

ptms <- bind_rows(DBPTM = dbptm, Phosphogrid = phosphogrid, .id = "ptm_source")

mimp <- as_tibble(readRDS(str_c(mutfunc_root, "ptm/mimp_results_yeast_lr1_all_possible_muts.rds"))) %>%
  select(gene, pos_aa = mut_pos, ref_aa, alt_aa, mimp_wt = score_wt,
         mimp_mut = score_mt, mimp_ratio = log_ratio, mimp_prob = prob, mimp_effect = effect)

# TFBS
tfbs <- as_tibble(readRDS("~/phd/mutfunc/tfbs/muts_chip_tfko_nofilter.rds")) %>%
  select(mut_id, tf, tf_ko_pval = ko_pval)

# Join
variants <- read_tsv("data/raw/variants_in_sets.txt") %>%
  select(mut_id = new_id, chr, pos, pos_ed, ref_allele, alt_allele, type) %>%
  left_join(mutfunc_ids, by = "mut_id") %>%
  left_join(uniprot_orf, by = c("gene" = "orf")) %>%
  left_join(sift, by = c("uniprot", "pos_aa", "ref_aa", "alt_aa")) %>%
  left_join(foldx, by = c("uniprot", "pos_aa", "ref_aa", "alt_aa")) %>%
  left_join(ints, by = c("uniprot", "pos_aa", "ref_aa")) %>%
  left_join(elm, by = c("gene", "pos_aa", "ref_aa", "alt_aa")) %>%
  left_join(ptms, by = c("gene", "pos_aa", "ref_aa")) %>%
  left_join(mimp, by = c("gene", "pos_aa", "ref_aa", "alt_aa")) %>%
  left_join(tfbs, by = "mut_id")
write_tsv(variants, "data/rm11_sk1_variants_mutfunc.tsv")
