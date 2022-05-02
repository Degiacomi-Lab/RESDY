# CARBAMYLATION FINDER

## Introduction

This software scans collections of protein structures, looking for lysines that may undergo a carbamlyation post-translational modification.

It is implemented in a set of Python classes, assembled as a pipeline in the Jupyter notebook `carbamylation.ipynb`.
A full description of the operations carried out by the pipeline is provided in the notebook. In short, the code will:

* Identify PDB or AlphaFold files from UNIPROT codes (see `Uniprot` class). Uniprot codes are either:
  - associated with an organism
  - provided manually
  - contained in an input .csv file
  
* download and curate each identified PDB file (see `PDB` class). Results are saved in the CSV file `result\proteins.csv`. Curation operations are:
  - mutation of MSE to MET
  - removal all HETATM, ions excluded
  - removal carboxylations from lysines
  - saving alternate conformations exist (e.g. NMR ensemble) in individual files
  - saving alternate side chain rotamers in individual files
  - addition of missing regions with Modeller. This operation is only allowed if size of gaps in sequence is smaller <8 amino acids, if larger the protein is disregarded.
  
* calculate pKa and solvent accessible area for every lysine in every curated structure (see `Measure` class). Results are saved in the CSV file `result\measures.csv`.
 
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
* biobox
* cython
* inquirer
* (jupyter)
## Notes

* The `Measure` class has been implemented to facilitate the addition of new measurable features. This is done by:
  - implementing a method taking a filename as input and returning a pandas DataFrame with three columns [resid, chain, feature].
  - adding the function name and its label in `self.measures` within `__init__`.

* The `Analysis` class is not yet implemented (partially implemented now).
  - `Analysis` class takes M.df rather than M.df.dropna() as the input in order to make both the aggregation and interactive plot work.
  - method `aggregate` will aggregate data by [Uniprot Entry, Resid] and return a dataframe with statistics produced for each [Uniprot Entry, Resid]
  - method `plot_graph` takes the type of plot ('histogram' or 'boxplot'), Uniprot Entry, Resid, and the feature ('pKa' or 'sasa') as inputs and produce the corresponding plot.
  - method `get_outliers` takes Uniprot Entry, Resid, the feature ('pKa' or 'sasa'), and a number 'whis' that controls the extremeness of the outliers (usually is 1.5) as inputs and return a dataframe with outlier data included.
  - method `interactive_plot` seems to be working (!!! I wish)
* planned feature for jupyter notebook: visualize protein structures using nglview.

