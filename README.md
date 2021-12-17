# Phenotype prediction of 1,011 *S. cerevisiae* strains using genotype, transcriptome and proteome data

This anlysis extends some of my previous [work](https://github.com/allydunham/hog_pathway), developing genotype to phenotype prediction methods based on omics data from 1,011 *S. cerevisiae* strains (see [Peter et al. 2018](https://www.nature.com/articles/s41586-018-0030-5)).
The transcriptomics data was provided by the [Schacherer lab](https://www.usias.fr/en/fellows/2017-fellows/joseph-schacherer/) and the proteomics data by the [Ralser lab](https://www.crick.ac.uk/research/labs/markus-ralser).
Neither dataset is currently available.
The growth phenotypes were measured by members of the Beltrao lab, where I performed this work, as detailed in [Galardini et al. (2019)](https://onlinelibrary.wiley.com/doi/abs/10.15252/msb.20198831).

This repo contains several phenotype prediction analyses, based on genotype (modelled with [P(Aff) scores](https://www.nature.com/articles/ng.1007)), gene expression and abundance scores, expressed as fold changes compared to that genes median expression:

* Associations and correlations between P(Aff), expression and abundance, showing generally weak relationships.
* Phenotype prediction based on linear models using the first 50 PCs of the P(Aff), abundance and expression scores.
* Variational Auto-Encoder based linear phenotype prediction models, using a custom VAE implementation.
* Gene based linear models, assessing the strength of association between each gene and phenotype, based on genotype, expression and abundance.

The project is split as follows:

* `bin`:
  * `analysis` - Scripts performing data analysis and figure generation
  * `processing` - Scripts to parse, format and normalise the raw data
  * `util` - Additional utility scripts
* `docs` - Two lists of genes identified
* `src` - shared modules and R config
