Simple scripts that extract qwidth data from bam files into 5 bp bins and stores them in a csv file. This way it's possible to quickly redo statistics using the csv file rather than having to reextract qwidth data for each run. First run Extract_qwidths.py then use the output from that to run cfDNA_fragment_size-PCA.py and cfDNA_group_ks_mwu_analysis_BRCA2.py.

## Dependencies

Python - version 3.9+
pandas (v2.2.3)
numpy (v2.0.2)
scipy (v1.13.1)
sckit-learn (v1.6.1)
pysam (v0.23.0)
seaborn (v0.13.2)
