# RESDY

[![Tests](https://github.com/Degiacomi-Lab/carbamylation/actions/workflows/tests.yml/badge.svg)](https://github.com/Degiacomi-Lab/carbamylation/actions/workflows/tests.yml)

## Introduction

This toolkit is subdivided in a set of classes that together operate as a pipeline enabling the rapid featurisation of aminoacids from collections of protein structures. Features and associated protein metadata can be explored with dedicated analysis and visualisation tools. In short, the pipeline will:

* Identify PDB or AlphaFold files from UNIPROT codes (see `Uniprot` class). UNIPROT codes are either:
  - associated with an organism
  - provided manually
  - contained in an input CSV file
  
* download and curate each identified PDB file (see `PDB` class). Results are saved in the CSV file `result\proteins.csv`. Curation operations are:
  - mutation of MSE to MET
  - removal of all HETATM, ions excluded
  - reversion to modified aminoacid to their wild type counterpart.
  - saving alternate conformations (e.g., NMR ensemble) in individual files
  - saving alternate side chain rotamers in individual files
  - addition of missing regions, if their size falls within a user-defined length (if larger the protein is disregarded).
  
* calculate a set of features for every amino acid of interest in every curated structure (see `Measure` class). Results are saved in the CSV file `result\measures.csv`.
 
* Plot aggregated data


## Installation

Clone the repository and install it from its root:

```
git clone https://github.com/Degiacomi-Lab/carbamylation.git
cd carbamylation
pip install -e .
```

The `-e` flag means an edit to a file under `src/resdy/` takes effect on the next import, with no reinstall. Drop it for a normal install.
This pulls in the required dependencies listed below, and makes the package importable from any working directory.

Several of the features rely on external programs rather than on Python packages, and are listed in the Optional table below. 
To check the installation:

```
python -c "import resdy; print(resdy.__version__)"
```

## Dependencies

### Required

The required dependencies are installed by `pip install -e .`, and are needed to run the overall pipeline:
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

Features denoted in the table below with specific package requirements or programme installations may not be available unless manually added to the working environment. A feature whose dependency is absent is skipped, with a message on the console, rather than stopping the run.

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
            <td><a href="https://ccsb.scripps.edu/mgltools/">from this website</a>
                <a href="https://anaconda.org/channels/conda-forge/packages/msms/overview">anaconda</a></td>
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
            <td><a href="https://github.com/HanaJaafari/Frustratometer">GitHub Repository</a></td>
        </tr>
        <tr>
            <td>openmm</td>
            <td><a href="https://anaconda.org/conda-forge/openmm">anaconda</a></td>
        </tr>
        <tr>
            <td>pdbfixer</td>
            <td><a href="https://anaconda.org/conda-forge/pdbfixer">anaconda</a></td>
        </tr>
        <tr>
            <td>15N NMR shift<br></td>
            <td>Legolas</td>
            <td><a href="https://github.com/roitberg-group/legolas">Legolas</a></td>
            <td>Follow the installation instructions at <a href="https://github.com/roitberg-group/legolas">GitHub repo</a>.</td>
        </tr>
        <tr>
            <td>Flexibility<br></td>
            <td>Biobox</td>
            <td>no extra requirements</td>
            <td></td>
        </tr>
        <tr>
            <td>Root Mean Square Fluctuation<br></td>
            <td>Biobox</td>
            <td>no extra requirements</td>
            <td></td>
        </tr>
        <tr>
            <td>Evolution<br></td>
            <td>ESM</td>
            <td>Fair-ESM</td>
            <td><a href="https://anaconda.org/channels/conda-forge/packages/fair-esm/overview">fair-esm</a></td>
        </tr>
        <tr>
            <td rowspan='2'>Secondary Structure<br></td>
            <td rowspan='2'>Biobox</td>
            <td>Biopython</td>
            <td><a href="https://anaconda.org/conda-forge/biopython">anaconda</a></td>
        </tr>
        <tr>
            <td>DSSP</td>
            <td><a href="https://anaconda.org/channels/conda-forge/packages/dssp/overview">DSSP</a></td>
        </tr>
        <tr>
            <td>Legolas (15N nmr prediction)<br></td>
            <td>Legolas</td>
            <td>Legolas</td>
            <td><a href="https://github.com/roitberg-group/legolas">Install through their GitHub</a></td>
        </tr>
    </tbody>
</table>
</div>

If running the Jupyter notebook `resdy.ipynb` (`pip install -e ".[notebook]"`), additional dependencies are:
* <a href="https://anaconda.org/conda-forge/seaborn">plotly</a>
* <a href="https://anaconda.org/conda-forge/seaborn">nglview</a>
* <a href="https://anaconda.org/conda-forge/jupyter">jupyter</a>

## Usage

After installation, the package is usable as follows.
```python
import resdy as RD

UP = RD.Uniprot()
```

`Uniprot`, `PDB`, `Measure`, `Analysis`, `Preprocessing`, `Aggregation` and `Viewer` are reachable straight from the package. 


* A tutorial demonstrating the pipeline main functionalities is provided in the notebook `resdy.ipynb`. 
* The API is available on readthedocs (LINK SOON).

## Technical Notes

* The `Uniprot` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/src/resdy/uniprot.py">src/resdy/uniprot.py</a>) is used to handle the collection of structural information for required proteins by mining the UNIPROT database. There are multiple ways to use this:
  - `from_csv_file()` function allow you to pass a list of UNIPROT codes as a csv to mine. Respective headers for the file should be from: 'Uniprot_Entry', PDB_Code', 'Resid'. 'Uniprot_Entry' is the only required column. Add data in 'PDB_Code' if you only require a specific PDB file from this protein. 'Resid' allows you to note a specific residue of interest within this protein, can be used for later analysis.
  - `get_protein_data()` function allows you to pass a UNIPROT code directly as a function input
  - `count_organism_proteins()` function allows you to pass a code for a whole proteome and extract information about all the proteins within this.

* The `Protein` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/src/resdy/protein.py">src/resdy/protein.py</a>) is used to handle extracting PDB and AlphaFold files for proteins specified within the Uniprot class. Structures are downloaded and patched to ensure good quality structures as used for calculations.


* The `Measure` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/src/resdy/measure.py">src/resdy/measure.py</a>) has been implemented to facilitate the addition of new measurable features. This is done by implementing a new measuring class in a new file, saved in the `src/resdy/features` folder.

* The `Analysis` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/src/resdy/analysis.py">src/resdy/analysis.py</a>) is responsible for aggregating the data (calculating mean, std, range for the feature of each lysine), scrape GO Terms from the Uniprot Database for each protein.
  - `Analysis` class takes the `Measure.df` pandas dataframe as the input.
  - method `aggregate` will aggregate data and calculate descriptive statistics for each lysine. The resulting dataframe will be stored in `self.df_aggregated`. 
  - method `subset` will aggregate data in two ways using either 'average' or 'south_east': (1) only include the average values of the two features for each lysine (2) only include the measure of the lysine with relatively lower pKa and higher sasa. The resulting dataframe will be stored in `self.df_sub`.
  - method `Go_Get_Data` will extract GO Terms associated with each distinct lysine from the Uniprot Database. This operation is sped up by applying multi-threading. The data will be stored in a dictionary (GO ID: a list of uniprot codes).

* The `Viewer` class (<a href="https://github.com/Degiacomi-Lab/carbamylation/blob/main/src/resdy/viewer.py">src/resdy/viewer.py</a>) can integrate different functionalities all together in an interactive plot.
  - `Viewer` class takes an instance of the `Analysis` class as the input.
  - method `enrichment_analysis` will do identify if there is any GO Term gets enriched within a group of proteins. The process of the enrichment analysis is: (1) for each GO Term, calculate an contingency table (2) compute the p-value for that contigency table and obtain a list of p-values in the end (3) adjust these p values using BH method to correct possible false postives.
  - method `advanced_plot` can make the interactive plot that integrates the scatter plot, enrichment analysis, and 3D visualisation of protein structures all together.
  - Note: the 3D visualisation requires downloaded PDB files.
