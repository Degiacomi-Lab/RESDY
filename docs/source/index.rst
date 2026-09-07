.. coolpackage documentation master file, created by
   sphinx-quickstart on Sat Jan 16 17:07:10 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to coolpackage documentation!
=====================================

This toolkit is subdivided in a set of classes that together operate as a pipeline enabling the rapid featurisation of amino acids from collections of protein structures. Features and associated protein metadata can be explored with dedicated analysis and visualisation tools.

The pipeline is constituted by five steps, each documented in its own section below.

1. **Uniprot**, gathering the list of protein structures to work on.
2. **Protein structure processing**, downloading those structures from the PDB and from AlphaFold, and patching them.
3. **Measure**, calculating one or more features for every residue of interest.
4. **Analysis**, exploring the resulting measurements.
5. **Viewer**, displaying them interactively.

.. note::

   Steps 4 and 5 are work in progress. Their documentation will be completed in a later release.

.. toctree::
   :maxdepth: 1
   :caption: 1. Uniprot

   uniprot

.. toctree::
   :maxdepth: 1
   :caption: 2. Protein structure processing

   protein
   patcher

.. toctree::
   :maxdepth: 1
   :caption: 3. Measure

   measure
   autoapi/index

.. toctree::
   :maxdepth: 1
   :caption: 4. Analysis (work in progress)

   analysis
   preprocessing
   aggregation

.. toctree::
   :maxdepth: 1
   :caption: 5. Viewer (work in progress)

   viewer

.. toctree::
   :maxdepth: 1
   :caption: Appendix

   FAQ


Please see coolpackage's `Github page <https://github.com/Degiacomi-Lab/carbamylation>`_
for installation instructions and example scripts.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
