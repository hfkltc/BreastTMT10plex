# BreastTMT10plex (Working in progress)
This repository includes codes used in "Deep Coverage of Global Protein Expression and Phosphorylation in Breast Tumor Cell Lines Using TMT 10-plex Isobaric Labeling" - http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00374

This readme.md explains how all of the scripts work and how they are connected.

The dataset includes the following folders and files:
- README.md 
- SILAC_Analysis_Functions.R
- phenotypicinfo.csv 
- ms2ProteinGroup
- ms2pSTY
- ms2pIP

SILAC_Analysis_Functions.R includes customrized functions

ms2ProteinGroup includes input files and R codes to analyze proteinGroup data acquired using QE-HF

ms2pSTY includes input files and R codes to analyze TiO2-enriched phosphopeptide data using QE-HF

ms2pYIP includes input files and R codes to analyze pY-IP-enriched phosphopeptide data using QE-HF
