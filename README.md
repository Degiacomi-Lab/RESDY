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

* The `Analysis` class is responsible for aggregating the data (calculating mean, std, range for the feature of each lysine), scrape GO Terms from the Uniprot Database for each unique protein.
  - `Analysis` class takes M.df as the input.
  - method `aggregate` will aggregate data and calculate descriptive statistics for each lysine. The resulting dataframe will be stored in `self.aggregated_df`. 
  - method `concise` will aggregate data in two ways using either 'average' or 'south_east': (1) only include the average values of the two features for each lysine (2) only include the measure of the lysine with relatively lower pKa and higher sasa, the trade-off between the two features can be controlled by a parameter 'weight'. The resulting dataframe will be stored in `self.df_concise`.
  - method `GO_Get_Data` will identify the unique proteins in the dataset and search for all the GO Terms associated with them from the Uniprot Database. The data will then be stored in a dictionary with the item in the strcture of (GO ID: a list of uniprot codes).

* The `CoolPlots` class can integrate different functionalities all together in an interactive plot.
  - `CoolPlots` class takes an instance of the `Analysis` class as the input.
  - method `enrichment_analysis` will do identify if there is any GO Term gets enriched within a group of proteins. The process of the enrichment analysis is: (1) for each GO Term, calculate an contingency table (2) compute the p value for that contigency table and obtain a list of p values in the end (3) adjust these p values using BH method to correct possible false-postives.
  - method `advanced_plot` can make the interactive plot that integrates the scatter plot, enrichment analysis, and 3D visualisation of protein structures all together.
