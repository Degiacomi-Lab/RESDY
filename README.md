# CARBAMYLATION FINDER

## Introduction

This software scans collections of protein structures, looking for lysines that may undergo a carbamlyation post-translational modification (PTM).

It is implemented in a set of Python classes, assembled as a pipeline in the Jupyter notebook `carbamylation.ipynb`.
A full description of the operations carried out by the pipeline is provided in the notebook.
In order to run the notebook without errors, please download the entire repository.
In short, the code will:

* Identify PDB or AlphaFold files from UNIPROT codes (see `Uniprot` class). UNIPROT codes are either:
  - associated with an organism
  - provided manually
  - contained in an input CSV file
  
* download and curate each identified PDB file (see `PDB` class). Results are saved in the CSV file `result\proteins.csv`. Curation operations are:
  - mutation of MSE to MET
  - removal of all HETATM, ions excluded
  - removal of carbamylations from lysines (revert to lysine)
  - saving alternate conformations (e.g. NMR ensemble) in individual files
  - saving alternate side chain rotamers in individual files
  - addition of missing regions with Modeller. This operation is only allowed if size of gaps in sequence is smaller <8 amino acids, if larger the protein is disregarded.
  
* calculate pKa and solvent accessible area for every lysine in every curated structure (see `Measure` class). Results are saved in the CSV file `result\measures.csv`.
 
* Plot aggregated data


## Dependencies

### Required

The following Python packages are required to run the overall pipeline:
* <a href="https://anaconda.org/conda-forge/numpy">numpy</a>
* <a href="https://anaconda.org/conda-forge/pandas">pandas</a>
* <a href="https://anaconda.org/conda-forge/biobox">biobox</a>
* <a href="https://anaconda.org/conda-forge/matplotlib">matplotlib</a>
* <a href="https://anaconda.org/conda-forge/seaborn">seaborn</a>
* <a href="https://anaconda.org/salilab/modeller">Modeller</a>

> **Note**
> Installation of Modeller requires a license key: https://salilab.org/modeller/registration.html.
> Note: there are issues with Modeller version 10.3, but any more recent version works correctly.

### Optional

There are different requirements for the methods for calculating different features for the lysines. The table below documents the dependencies for each.

<div class="table_component" role="region" tabindex="0">
<table>
    <caption><br></caption>
    <thead>
        <tr>
            <th>Feature<br></th>
            <th>Methods<br></th>
            <th>Requirements</th>
            <th>Where to find</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td rowspan='2'>pKa<br></td>
            <td>PROPKA3</td>
            <td>propka</td>
            <td><a href="https://anaconda.org/conda-forge/propka">anaconda</a></td>
        </tr>
        <tr>
            <td>pKaANI</td>
            <td>pKaANI</td>
            <td><a href="https://github.com/isayevlab/pKa-ANI">pKaANI GitHub Repository</a></td>
        </tr>
        <tr>
            <td rowspan='2'>Depth<br></td>
            <td rowspan='2'>Biopython</td>
            <td>Biopython</td>
            <td><a href="https://anaconda.org/conda-forge/biopython">anaconda</a></td>
        </tr>
        <tr>
            <td>MSMS</td>
            <td><a href="https://ccsb.scripps.edu/mgltools/">from this website</a></td>
        </tr>
        <tr>
            <td>Solvent Accessible Surface Area<br></td>
            <td>Biobox</td>
            <td>no extra requirements</td>
            <td></td>
        </tr>
        <tr>
            <td rowspan='4'>Atomic Environment Vectors<br></td>
            <td rowspan='4'>ANI-2x</td>
            <td>ase</td>
            <td>anaconda</td>
        </tr>
        <tr>
            <td>torch</td>
            <td><a href="https://anaconda.org/pytorch/pytorch">anaconda</a></td>
        </tr>
        <tr>
            <td><a href="https://github.com/aiqm/torchani/">torchani</a></td>
            <td><a href="https://anaconda.org/conda-forge/torchani">anaconda</a></td>
        </tr>
        <tr>
            <td>cuaev (optional: for gpu accelerated calculation)</td>
            <td><a href="https://github.com/aiqm/torchani/tree/master/torchani/cuaev">GitHub Repository</a></td>
        </tr>
        <tr>
            <td>Dynamically Accessible Surface<br></td>
            <td>Biobox</td>
            <td>no extra requirements</td>
            <td></td>
        </tr>
        <tr>
            <td>Sequence Charge<br></td>
            <td>Biobox</td>
            <td>Biopython</td>
            <td><a href="https://anaconda.org/conda-forge/biopython">anaconda</a></td>
        </tr>
        <tr>
            <td>Curvature, Writhing, Torsion, arc_length, phi, psi<br></td>
            <td>Melodia</td>
            <td><a href="https://github.com/rwmontalvao/Melodia_py">Melodia-py</a></td>
            <td><a href="https://github.com/rwmontalvao/Melodia_py">GitHub Repository</a> or pip</td>
        </tr>
        <tr>
            <td rowspan='3'>Frustration, density<br></td>
            <td rowspan='3'>Frustratometer</td>
            <td><a href="https://github.com/HanaJaafari/Frustratometer?tab=readme-ov-file">Frustratometer</a></td>
            <td><a href="https://github.com/rwmontalvao/Melodia_py">GitHub Repository</a></td>
        </tr>
        <tr>
            <td>openmm</td>
            <td><a href="https://anaconda.org/conda-forge/openmm">anaconda</a></td>
        </tr>
        <tr>
            <td>pdbfixer</td>
            <td><a href="https://anaconda.org/conda-forge/pdbfixer">anaconda</a></td>
        </tr>
    </tbody>
</table>
</div>

If running the Jupyter notebook `carbamylation.ipynb`:
* plotly
* nglview
* jupyter

## Technical Notes

* The `Uniprot` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/uniprot.py">uniprot.py</a>) is used to handle collecting the structural information for the proteins required by mining the UNIPROT database. There are multiple ways to use this:
  - from_csv_file() function allow you to pass a list of UNIPROT codes as a csv to mine. Respective headers for the file should be from: 'Uniprot_Entry', PDB_Code', 'Resid'. 'Uniprot_Entry' is the only required column. Add in data in the 'PDB_Code' if you only require a specific PDB file from this protein. 'Resid' allows you to note a specific residue of interest within this protein, can be used for later analysis.
  - get_protein_data() function allows you to pass a UNIPROT code directly as a function input
  - count_organism_proteins() function allows you to pass a code for a whole proteome and extract information about all the proteins within this.

* The `Protein` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/protein.py">protein.py</a>) is used to handle extracting PDB files for proteins specified within the Uniprot class. Structures are downloaded and patched to ensure good quality structures as used for calculations.

* THe `Alphafold` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/alphafold.py">alphafold.py</a>) is used to handle extracting AlphaFold data for proteins specified within the Uniprot class. Structures are downloaded and patched to ensure good quality structures as used for calculations.

* The `Measure` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/measure.py">measure.py</a>) has been implemented to facilitate the addition of new measurable features. This is done by:
  - implementing a method taking a filename as input and returning a pandas DataFrame with three columns [resid, chain, feature].
  - adding the function name and its label in `self.measures` within `__init__`.

* The `Analysis` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/analysis.py">analysis.py</a>) is responsible for aggregating the data (calculating mean, std, range for the feature of each lysine), scrape GO Terms from the Uniprot Database for each protein.
  - `Analysis` class takes `Measure.df` as the input.
  - method `aggregate` will aggregate data and calculate descriptive statistics for each lysine. The resulting dataframe will be stored in `self.df_aggregated`. 
  - method `subset` will aggregate data in two ways using either 'average' or 'south_east': (1) only include the average values of the two features for each lysine (2) only include the measure of the lysine with relatively lower pKa and higher sasa. The resulting dataframe will be stored in `self.df_sub`.
  - method `Go_Get_Data` will extract GO Terms associated with each distinct lysine from the Uniprot Database. This operation is sped up by applying multi-threading. The data will be stored in a dictionary (GO ID: a list of uniprot codes).

* The `Viewer` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/viewer.py">viewer.py</a>) can integrate different functionalities all together in an interactive plot.
  - `Viewer` class takes an instance of the `Analysis` class as the input.
  - method `enrichment_analysis` will do identify if there is any GO Term gets enriched within a group of proteins. The process of the enrichment analysis is: (1) for each GO Term, calculate an contingency table (2) compute the p-value for that contigency table and obtain a list of p-values in the end (3) adjust these p values using BH method to correct possible false postives.
  - method `advanced_plot` can make the interactive plot that integrates the scatter plot, enrichment analysis, and 3D visualisation of protein structures all together.
  - Note: the 3D visualisation requires downloaded PDB files.
