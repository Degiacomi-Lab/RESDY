# CARBAMYLATION FINDER

## Introduction

This software scans collections of protein structures, looking for lysines that may undergo a carbamlyation post-translational modification.
An example of usage is available in the Jupyter notebook `carbamylation.ipynb`

 The code will:
* identify PDB or AlphaFold files from UNIPROT codes either:
  - associated with an organism
  - provided manually
  - contained in an input .csv file
  
* download and curate each identified PDB file. Curation operations are:
  - mutation of MSE to MET
  - removal all HETATM, ions excluded
  - removal carboxylations from lysines
  - saving alternate conformations exist (e.g. NMR ensemble) in individual files
  - saving alternate side chain rotamers in individual files
  - addition of missing regions with Modeller. This operation is only allowed if size of gaps in sequence is smaller <8 amino acids, if larger the protein is disregarded.
  
* calculate pKa and solvent accessible area for every lysine in every curated structure
 
* Plot scatter plots aggregating all data
 

## Dependencies

The following Python packages are required:
* numpy
* pandas
* matplotlib
* bs4
* modeller
* propka
* plotly
* nglview


## Notes

* planned refactoring and new features:
  - Convert postprocessing operations into a class (Analysis), and refactor
  - code works by chdir into working directory takes place, instead of referring to different working folders as subfolders of results
  - pickle instances of Uniprot and PDB classes in results folder
  - visualize protein structures using nglview within the Jupyter notebook

