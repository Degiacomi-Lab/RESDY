.. RESDY documentation master file, created by
   sphinx-quickstart on Sat Jan 16 17:07:10 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RESDY documentation!
===============================

This toolkit is subdivided in a set of classes that together operate as a pipeline enabling the rapid featurisation of amino acids from collections of protein structures. Features and associated protein metadata can be explored with dedicated analysis and visualisation tools.
The pipeline is constituted by five steps, each documented in its own section below.

.. note::
   Steps 4 and 5 are work in progress. Their documentation will be completed in a later release.

1\. **Uniprot**, gathering the list of protein structures to work on.

.. toctree::
   :maxdepth: 1

   uniprot

2\. **Protein structure processing**, downloading those structures from the PDB and from AlphaFold, and patching them.

.. toctree::
   :maxdepth: 1

   protein
   patcher

3\. **Measure**, calculating one or more features for every residue of interest.

.. toctree::
   :maxdepth: 1

   measure
   autoapi/index

4\. **Analysis**, exploring the resulting measurements.

.. toctree::
   :maxdepth: 1

   analysis
   preprocessing
   aggregation

5\. **Viewer**, displaying them interactively.

.. toctree::
   :maxdepth: 1

   viewer

**Appendix**

.. toctree::
   :maxdepth: 1

   FAQ


Please see RESDY's `Github page <https://github.com/Degiacomi-Lab/RESDY>`_
for installation instructions and example scripts.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
