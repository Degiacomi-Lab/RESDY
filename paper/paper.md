---
title: 'RESDY: a featurisation toolkit of protein sidechain structure and properties for classification tasks'
tags:
  - Python
  - proteins
  - amino acids
  - featurization
  
authors:
  - name: George Weston
    orcid: 0009-0000-1924-7183
    equal-contrib: false
    affiliation: 1
  - name: Matteo T. Degiacomi
    orcid: 0000-0003-4672-471X
    equal-contrib: false
    corresponding: true
    affiliation: 2, 3
affiliations:
 - name: Department of Biosciences, Durham University, United Kingdom
   index: 1
 - name: School of Informatics, University of Edinburgh, United Kingdom
   index: 2
 - name: EaStCHEM School of Chemistry, University of Edinburgh, United Kingdom
   index: 3
date: 1 September 2026
bibliography: paper.bib

---

# Summary

We present `resdy` [TODO].

# Statement of need

Assigning numerical descriptors to individual amino acids within a protein structure is a task carried out repeatedly, and largely independently, across structural bioinformatics.
It underpins the prediction of post-translational modification sites [@huang2013], the identification of catalytic residues, the ranking of residues for covalent ligand design, and the assessment of chemical liabilities in biologics.
In each case the practitioner assembles a broadly similar codebase, involving a means of gathering relevant atomic structures, tools for cleaning them, a set of calls to third-party descriptor programs, and a means of reducing the result to a table suitable for a classifier.

As a first, practical hurdle, building a datasets of interest often involves processing thousands of structures, so that computational tractability requires the whole procedure to be streamlined. More profoundly though, several properties of the underlying data make this task substantially harder than it appears.
First, many deposited protein structures should not be considered as directly usable, as they may contain mutations, unresolved regions, alternate side chain conformations, non-standard aminoacids, cofactors, or other macromolecular binding partners. Those best placed to exploit the resulting datasets are not always those most familiar with the conventions of structural data, leading to potentially unsubstantiated curation decisions.
Second, proteins are flexible, whereby any protein structure is only a representative of a broader conformational ensemble. As a result, the value of a descriptor computed for an aminoacid will vary depending on the specific conformational representative the user retrieved.
Third, structural coverage is supplied by sources that may not be fully compatible. The protein databank features protein models obtained from x-ray crystallography, Nuclear Magnetic Resonance, and cryo-electron microscopy measurements, whereas AI engines such as AlphaFold [@alphafold] generate, and subsequently energy minimise in vacuo, models based on a training set primarily constituted by x-ray data.

Existing software addresses parts of the problems highted above.
ProCaliper [@procaliper] retrieves annotations and structures from UniProt [@uniprot] and the AlphaFold Protein Structure Database [@afdb], and computes residue-level charge, solvent accessibility and protonation state, but operates on a single structure per protein and does not address curation or the comparability of experimental and predicted sources.
structuremap [@structuremap] places post-translational modifications in structural context at proteome scale, yet considers predicted models exclusively.
Graphein [@graphein] and PyUUL [@pyuul] convert structures into graphs, voxels and point clouds for geometric deep learning, and therefore target learned rather than interpretable per-residue quantities.
ProtDCal [@protdcal] computes large batteries of residue-level indices, though from a single theoretical framework and through a Java graphical interface, while iFeatureOmega [@ifeatureomega] and ProFeatX [@profeatx] provide extensive featurisation that is predominantly sequence-derived.
Orthogonal to these, individual programs are explicitly developed to produce a single descriptor.
Among them we find PROPKA [@propka] and pKa-ANI [@pkaani] for protonation, MSMS [@msms] for molecular surfaces, frustratometer [@frustratometer] for local energetic frustration, Melodia [@melodia] for backbone differential geometry, LEGOLAS [@legolas] for NMR chemical shifts, TorchANI [@torchani] for atomic environment vectors, and ESM-2 [@esm2] for language model embeddings.
None of them, however, is concerned with where the structures came from, with what state they were in on arrival, or with how many of them describe the same residue.

`resdy` is designed to address these concerns.
Protein structure curation is aimed at ensuring the data gathered is a close representative of the ensemble of conformation a protein occupies in its native environment.
As such, for each protein of interest every deposited entry, NMR candidate model, and side chain alternate conformation is curated and featurised separately, so that the dispersion of a descriptor across structures can become itself informative.
Curation involves modelling missing regions within a user-defined gap length, retaining ions and discarding heteroatoms, and modifying mutated residues to their wild type counterpart, which is a prerequisite for asking whether they would be modified at all. Furthermore, an optional energy minimisation in implicit solvent is available to ensure input structures are relaxed in a suitable proxy for their native environment.
Finally, to ensure the user can assess any residual dependence on provenance, metadata including determination method and resolution are saved alongside every measurement.

To ensure that datasets can be assembled at scale, within `resdy` curation and featurisation run in parallel, per-structure failures are recorded rather than allowed to interrupt a run, and the aggregated output is prepared for downstream learning through dedicated filtering and data aggregation tools.
As the residue type(s) of interest and the atoms defining its reactive centre are supplied as parameters, this modelling pipeline may be applied to any residue, be it standard of modified.

# Package Description

`resdy` is subdivided in a set of classes that together operate as a pipeline enabling the rapid featurisation of aminoacids from collections of protein structures. Features and associated protein metadata can be explored with dedicated analysis and visualisation tools (see Figure \autoref{fig:gui}).

![here we could maybe have that artistic representation of the pipeline, minus the classifier? \label{fig:greatfigure}.](greatfigure.png)

The main components of the pipeline are as follows.

*	*Uniprot*. Given a UNIPROT proteome identifier or a list of UNIPROT accession codes [@uniprot], this class fetches all associated Protein Data Bank [@pdb] and AlphaFold Protein Structure Database [@afdb] codes.
*	*Protein*. Given a list of raw protein structures, this class downloads and curates them.
Curation includes modelling missing regions (using Modeller [@modeller]), removing cofactors, and ensuring only standard aminoacids are present.
If alternative conformations or side chain rotamers are available in the file, these are split in individual structures.
For AlphaFold models, an optional energy minimisation step in implicit solvent using OpenMM [@openmm] is also available.
If multiple cores are available on the computer running `resdy`, all these operations can be carried out in parallel.
*	*Measure*. Given a list of curated protein structures, this class featurises all the aminoacids of interest, optionally in parallel.
While a range of features are already made available (see below), the code architecture has been designed to facilitate the addition of custom features.
This class yields a pandas [@pandas] dataframe coupling the features extracted with relevant metadata (e.g., structure determination technique and resolution).
*	*Analysis*. This class enables exploring and manipulating the data gathered as a pandas dataframe by the Measure class.
TALK ABOUT DATA AGGREGATION.
*	*Viewer*. TBD

The features currently available within `resdy`, with associated origin packages given in parentheses, are:

* Solvent accessible surface area (biobox [@biobox], after Shrake and Rupley [@shrake])
* Dynamically accessible surface (biobox [@biobox])
* Depth (Biopython [@biopython] and MSMS [@msms])
* pKa (PROPKA3 [@propka] or pKa-ANI [@pkaani])
* Atomic environment vectors (ANI-2x [@ani2x] via TorchANI [@torchani])
* Protein language model embeddings (ESM-2 [@esm2])
* Sequence charge (biobox [@biobox])
* Local energetic frustration and density (frustratometer [@frustratometer; @frustratometer_r])
* backbone 15N NMR chemical shift (LEGOLAS [@legolas])
* Curvature, writhing, torsion, arc length, phi and psi (Melodia [@melodia])
* Mean B-factor over the atoms of the residue (biobox [@biobox])
* Root mean square fluctuation across the conformations available for a structure (biobox [@biobox])

# Usage

To exemplify the usage of `resdy` we gather, curate, and featurise all the proteins of the organism **Organismus importantissimus**. While `resdy` allows processing unreviewed UNIPROT codes, its default behaviour is to only process reviewed ones.

```
import resdy as RD

UP = RD.Uniprot()
UP.get_organism_proteins(code='UP000001811')

P = RD.PDB(gap=10)
P.gather_proteins(uniprot_df=UP.df)

M = RD.Measure(df_input=P.df, residue_of_interest='LYS',
            features=['depth', 'sasa', 'propka', 'aev'])
M.measure_data()
M.save_state()

#add a line from analyser, that we can display below?
```

The code gathers … reviewed UNIPROT codes, associated with … proteins and … individually featurized aminoacids.
On a computer with … CPUs and a … GPU this terminates in ….

![here we could maybe display some violin plots of the output features? \label{fig:greatfigure2}.](greatfigure2.png)




# Acknowledgements

We thank Hao Man, Breanna Voss, and Grace Carter for testing the code, and Martin Cann for sharing his expertise in the area of protein modifications.

# References