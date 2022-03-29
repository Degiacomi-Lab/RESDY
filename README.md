# CARBAMYLATION FINDER

## Introduction

This software scans collections of protein structures, looking for lysines that may undergo a carbamlyation post-translational modification.
An example of usage is available in the Jupter notebook `carbamylation.ipynb`

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
 
Further insight on produced data can be obtained by observing the data in the Jupyter notebook `carbamylation.ipynb`.


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

* known bugs:
  - a known bug in the inquirer package makes it non functional in Python >3.7 under Windows. Current workaround involves using the numerical pad instead of arrows.
  - missing implementation for options in postprocessing.py when launched alone

* planned refactoring:
  - Convert measuring and postprocessing into classes
  - edit code so that a chdir into working directory takes place, instead of referring to different working folders as subfolders of results
  
* planned new features
  - save classes instances of Uniprot and PDB in results data
  - visualize protein structures using nglview within the Jupyter notebook
  - simplify addition of additional lysine scoring metrics (e.g. depth, sequence conservation).


