This is a script that converts Nanopore bam files to the 450k methylation array format that is compatible with meth atlas.
The generated file should be directly compatible with the MethAtlas tool: https://github.com/nloyfer/meth_atlas. 
In order to run the script a reference manifest will be necessary that can be downloaded form here: https://github.com/methylgrammarlab/cfdna-ont/tree/main/deconvolution_code.

## Dependencies
- Python (3.0+)
  - Deconvolve.py (this is a program that should be included in the meth_atlas directory you git cloned)
  - "Infinium_HumanMethylation450k_manifests" (either hg38, hg19 or custom)
  - Python packages:
    - pandas (v2.2.3)
    - numpy (v2.0.2)
    - scikit-learn (v1.6.1)
    - seaborn (v0.13.2)
    - matplotlib (v3.9.1)
    - scipy (1.13.1)
