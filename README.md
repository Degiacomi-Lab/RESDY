# CARBAMYLATION FINDER

This software scans protein structure, looking for lysines that may undergo a carbamlyation post-translational modification.
To run it, type the following in a terminal:

python carbamylation.py
 
and respond to the question the code asks you.
The code will identify PDB or AlphaFold files from UNIPROT codes, and calculate pKa a Solvent accessible area for every lysine.
If multiple atomic structures exist for the same UNIPROT code, or a lysine exists in many alternative conformations, calculations are run for each of them.

Plots analysing data are produced at the end of the run. Further insight can be obtained by observing the data in the Jupyter notebook carbamylation.ipynb.

*** this software is experimental ***

## TO DO

* call PDB loader from carbamylation.py, and not from uniprot_loader, i.e. uniprot loader only helps creating the database of codes
* AlphaFold download, and PDB download + cleaning + patching should happen at the same time
* removing hydrogens before splitting would speed things up
* missing implementation for options in postprocessing.py when launched alone
* fix separator issue in Windows for CSV file user input (add path autocomplete)
 
