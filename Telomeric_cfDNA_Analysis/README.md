This analysis relies on TelSeq which can be approached here: https://github.com/zd1/telseq. The TelSeq output can then be processed with the telseq_parser.py script to process the raw output and then run telseq_metadata_integration.py and then telseq_visualization.py
These scripts will only work if you have a metadata file with age, gender and BMI if any of these are missing it should be possible to remove certain aspects from these scripts to enable them to work regardless but it won't be possible to replicate the study directly.

## Dependencies
- TelSeq (v0.0.1)
- GCC (≥v4.8)
- bamtools (v2.5.3)
- python (v3.9+)
  - pandas (v2.2.3)
  - numpy (v2.0.2)
  - scipy (v1.13.1)
  - sckit-learn (v1.6.1)
  - seaborn (v0.13.2)
