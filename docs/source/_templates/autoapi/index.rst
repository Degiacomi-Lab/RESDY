.. This file is the AutoAPI landing page template. It overrides the one shipped with
   sphinx-autoapi (see ``autoapi_template_dir`` in conf.py). Edit the title, the prose and
   the table below freely. Only the toctree loop at the bottom is machinery, everything
   else is plain reStructuredText. Note that Jinja delimiters are parsed even inside
   reStructuredText comments such as this one, so do not write any here.

   The table maps the keys accepted by the "features_dict" argument of the Measure class
   onto the class that computes them. It is written by hand, because the mapping lives in
   Measure._setup_measures() and cannot be recovered from the feature modules themselves.
   Keep the two in step when a feature is added or renamed: a class target that no longer
   exists does not raise a warning, it just renders as plain text without a link.

   The toctree loop iterates over the modules AutoAPI reports rather than over top-level
   objects: "resdy" and "resdy.features" are packages, so the only top-level
   object AutoAPI reports is "resdy" itself and filtering on is_top_level_object
   would yield a single row. The two package pages are dropped through exclude_patterns in
   conf.py, feature.py through autoapi_ignore, and error_reporting is filtered out of the
   loop below, so that the table links straight to each feature.

Feature measurements
====================

Each module of the ``features`` folder defines a single class, dedicated to the measurement of one
specific feature. All these classes are called by the :class:`Measure <resdy.measure.Measure>` class,
which orchestrates the featurisation of a collection of protein structures.

Features are requested through the ``features_dict`` argument of
:class:`Measure <resdy.measure.Measure>`, a dictionary with one entry per feature. The key is the
name of the feature and the value is a dictionary of keyword arguments passed to the constructor of
its class, left empty to keep the defaults, for instance
``Measure(df_input=pdb.df, features_dict={'propka': {}, 'sasa': {'probe': 1.4}, 'rmsf': {}})``.
A list of names is not accepted. The name is looked up as the upper-cased class name, so that any
class added to the ``features`` folder can be requested by its own name. The table below lists the
names accepted, and the class that computes each.

.. list-table::
   :header-rows: 1
   :widths: 18 18 64

   * - Name
     - Class
     - Measurement
   * - ``propka``
     - :py:class:`~resdy.features.PROPKA.PROPKA`
     - pKa of the residues of interest, calculated with PROPKA3.
   * - ``pkaANI``
     - :py:class:`~resdy.features.PKAANI.PKAANI`
     - pKa of the residues of interest, calculated with pKaANI.
   * - ``sasa``
     - :py:class:`~resdy.features.SASA.SASA`
     - Solvent accessible surface area of the side chain of the residue of interest, calculated on
       a substructure cut out around its anchor atom (NZ for lysine).
   * - ``depth``
     - :py:class:`~resdy.features.DEPTH.DEPTH`
     - Distance of the residue from the protein surface.
   * - ``das``
     - :py:class:`~resdy.features.DAS.DAS`
     - Dynamically accessible surface of the anchor atom of the residue (NZ for lysine), i.e. the
       number of positions it can take within the structure.
   * - ``aev``
     - :py:class:`~resdy.features.AEV.AEV`
     - ANI-2x atomic environment vector of the NZ atom, a 1008-element description of its local
       environment.
   * - ``legolas``
     - :py:class:`~resdy.features.LEGOLAS.LEGOLAS`
     - 15N NMR chemical shift, calculated with LEGOLAS.
   * - ``seqcharge``
     - :py:class:`~resdy.features.SEQCHARGE.SEQCHARGE`
     - Summed charge of the amino acids on either side of the residue of interest.
   * - ``flexibility``
     - :py:class:`~resdy.features.FLEXIBILITY.FLEXIBILITY`
     - B-factor normalised over the structure, averaged over the atoms of the residue of interest.
   * - ``frustration``
     - :py:class:`~resdy.features.FRUSTRATION.FRUSTRATION`
     - Frustration metric, calculated with Frustratometer.
   * - ``density``
     - :py:class:`~resdy.features.FRUSTRATION.FRUSTRATION`
     - Local density, calculated with Frustratometer alongside the frustration metric.
   * - ``melodia``
     - :py:class:`~resdy.features.STRUCTURE.STRUCTURE`
     - All the Melodia structural descriptors at once.
   * - ``curvature``, ``writhing``, ``torsion``, ``arc_length``, ``phi``, ``psi``
     - :py:class:`~resdy.features.STRUCTURE.STRUCTURE`
     - The individual Melodia structural descriptors, requested one by one.
   * - ``evolution``
     - :py:class:`~resdy.features.EVOLUTION.EVOLUTION`
     - ESM language model embedding of the sequence around the residue of interest.
   * - ``rmsf``
     - :py:class:`~resdy.features.RMSF.RMSF`
     - Root mean square fluctuation of the residue across all the structures curated for its
       Uniprot entry.
   * - ``secondarystructure``
     - :py:class:`~resdy.features.SECONDARYSTRUCTURE.SECONDARYSTRUCTURE`
     - DSSP secondary structure assignment of the residue of interest, as a letter or a number.

``features_dict={'all': {}}`` is a shorthand for every feature in the table except
``secondarystructure`` and the individual Melodia descriptors (which ``melodia`` covers), each with
its default settings. The list is ``ALL_FEATURES`` in :mod:`resdy.measure`.

.. toctree::
   :titlesonly:
   :hidden:

{% for page in pages|selectattr("type", "equalto", "module")|sort(attribute="name") %}
{% if not page.name.endswith("error_reporting") %}
   {{ page.include_path }}
{% endif %}
{% endfor %}
