#!/usr/bin/env Rscript
# Generate P(Aff) type scores
source('src/config.R')

# Rscript bin/paff.R data/variants/chr1.tsv data/variants/chr1.genotypes
args = commandArgs(trailingOnly=TRUE)

# Import data
annotation <- read_tsv(args[1], col_types = cols(PROTEINLOC=col_character())) %>%
  rename_all(str_to_lower) %>%
  mutate(uniprot = unname(NAME_TO_UNIPROT[geneid]))

sift <- read_rds('data/raw/sift_all_yeast.rds') %>%
  filter(acc %in% annotation$uniprot) %>%
  mutate(pos = as.character(pos))

# Add sift scores to annotation and drop:
# - Substitutions with no score
# - Genes with no sift results (these are generally pseudogenes) 
annotation <- left_join(annotation, select(sift, uniprot=acc, proteinloc=pos, refaa=ref, varaa=alt, score),
                        by = c("uniprot", "proteinloc", "refaa", "varaa")) %>%
  filter(!(consequence %in% c('nonsynonymous', 'synonymous') & is.na(score)), uniprot %in% sift$acc)

# Add consequences to genotypes and filter to those with consequences
genotypes <- read_tsv(args[2]) %>%
  pivot_longer(!one_of(c('chromosome', 'position', 'ref', 'alt')), names_to = 'strain', values_to = 'genotype') %>%
  filter(!genotype == 0) %>%
  left_join(select(annotation, chromosome = seqnames, position = start, ref, alt, uniprot, geneid, proteinloc, refaa, altaa=varaa, consequence, score),
            by = c('chromosome', 'position', 'ref', 'alt')) %>%
  drop_na(uniprot)

# Calulcate P(Aff) scores
# P(Aff) = 1 - Prod(P(Neut))
# P(Neut) rules:
# - From SIFT: 1/(1 + exp(-1.312424 * ln(sift_score + 1.598027e-05) - 4.103955))
# - Frameshift: 0.05
# - Nonsense: 0.05
calc_p_neut <- function(type, sift){
  pneut <- 1/(1 + exp(-1.312424 * log(sift + 1.598027e-05) - 4.103955))
  pneut[type %in% c('frameshift', 'nonsense')] <- 0.05
  return(pneut)
}

scores <- mutate(genotypes, pneut = calc_p_neut(consequence, score)) %>%
  group_by(strain, uniprot, geneid) %>%
  summarise(paff = 1 - prod(pneut), .groups = 'drop')

write(format_tsv(scores), file = stdout())
