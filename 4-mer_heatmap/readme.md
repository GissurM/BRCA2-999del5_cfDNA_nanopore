This script is based on FragmentomicsGenomBiol. Download the directory and set up the run_fragmentomic_analysis.sh and then run it using bash. Doing so should fully run a pipeline that also includes nanopore_fragmentomics.py and count_motif.R. This should create motif.csv files which contain information on the percentage of fragments with each 4-mer end motif these files should also be usable for downstream statistics.

## Dependencies
- Python 3.9+
  - pandas (v2.2.3)
  - numpy (v2.0.2)
  - pysam (v0.23.0)
- R (v4.3.3)
  - data.table (v1.18.2.1)

Python packages used for statistics scripts only
- scipy (v1.13.1)
- sklearn (v1.6.1)
- seaborn (v0.13.2)
- matplotlib (v3.9.1)
