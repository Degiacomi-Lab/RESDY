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
 
* Plot aggregated data


## Dependencies

The following Python packages are required:
* numpy
* pandas
* matplotlib
* cython
* biobox
* bs4
* modeller
* propka
* plotly
* nglview
* jupyter (optional, required to execute the notebook `carbamylation.ipynb`)

## Technical Notes

* The `Measure` class has been implemented to facilitate the addition of new measurable features. This is done by:
  - implementing a method taking a filename as input and returning a pandas DataFrame with three columns [resid, chain, feature].
  - adding the function name and its label in `self.measures` within `__init__`.

* The `Analysis` class is responsible for aggregating the data (calculating mean, std, range for the feature of each lysine), scrape GO Terms from the Uniprot Database for each protein.
  - `Analysis` class takes M.df as the input.
  - method `aggregate` will aggregate data and calculate descriptive statistics for each lysine. The resulting dataframe will be stored in `self.df_aggregated`. 
  - method `subset` will aggregate data in two ways using either 'average' or 'south_east': (1) only include the average values of the two features for each lysine (2) only include the measure of the lysine with relatively lower pKa and higher sasa. The resulting dataframe will be stored in `self.df_sub`.
  - method `Go_Get_Data` will extract GO Terms associated with each distinct lysine from the Uniprot Database. This operation is sped up by applying multi-threading. The data will be stored in a dictionary (GO ID: a list of uniprot codes).

* The `Viewer` class can integrate different functionalities all together in an interactive plot.
  - `Viewer` class takes an instance of the `Analysis` class as the input.
  - method `enrichment_analysis` will do identify if there is any GO Term gets enriched within a group of proteins. The process of the enrichment analysis is: (1) for each GO Term, calculate an contingency table (2) compute the p-value for that contigency table and obtain a list of p-values in the end (3) adjust these p values using BH method to correct possible false postives.
  - method `advanced_plot` can make the interactive plot that integrates the scatter plot, enrichment analysis, and 3D visualisation of protein structures all together.
  - note the 3D visualisation requires downloaded PDB files.

# Installation notes for Mac environment:
- Installing modeller:

`conda install -c salilab modeller=10.2`
(There are issues with version 10.3)

You will also need to get a modeller licence key: https://salilab.org/modeller/registration.html

- Additional dependencies for Jupyter notebook: 
   - seaborn
