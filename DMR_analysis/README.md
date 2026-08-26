A collection of scripts that are meant for DMR analysis 3 of the scripts perform a group pairwise analysis between the BRCA2 and control samples with each of them focusing on different parts of the genome, BRCA2_genomewide_methylation_permanova_promoters.py focuses on promoter regions, BRCA2_genomewide_methylation_permanova_filler.py focuses on filler regions and BRCA2_genomewide_methylation_permanova_genebodies.py focuses on genebodies.
A seperate version focuses on doing a global comparison of each genomic position between the two groups and determining in which of these positions any difference in genomic methylation appears at all.

## Dependencies
- python (v3.9+)
  - pandas (v2.2.3)
  - numpy (v2.0.2)
  - scipy (v1.13.1)
  - matplotlib (v3.9.1)
  - requests (v2.32.5)
  - seaborn (v0.13.2)
