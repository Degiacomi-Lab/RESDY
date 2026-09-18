.. This file is the AutoAPI landing page template. It overrides the one shipped with
   sphinx-autoapi (see ``autoapi_template_dir`` in conf.py). Edit the title, the prose and
   the table below freely. Only the toctree loop at the bottom is machinery, everything
   else is plain reStructuredText. Note that Jinja delimiters are parsed even inside
   reStructuredText comments such as this one, so do not write any here.

   The table maps the strings accepted by the "features" argument of the Measure class
   onto the class that computes them. It is written by hand, because the mapping lives in
   Measure._setup_measures() and cannot be recovered from the feature modules themselves.
   Keep the two in step when a feature is added.

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

A feature is requested by passing its name to the ``features`` argument of
:class:`Measure <resdy.measure.Measure>`, for instance
``Measure(df_input=pdb.df, features=['propka', 'sasa', 'rmsf'])``. The table below lists every name
accepted, and the class that computes it.

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
     - :py:class:`~resdy.features.pkaani.PKAANI`
     - pKa of the residues of interest, calculated with pKaANI.
   * - ``sasa``
     - :py:class:`~resdy.features.sasa.SASA`
     - Solvent accessible surface area of the NZ atom, calculated on a substructure of the atoms
       within 15 Angstrom of it.
   * - ``depth``
     - :py:class:`~resdy.features.depth.Depth`
     - Distance of the residue from the protein surface.
   * - ``das``
     - :py:class:`~resdy.features.das.DAS`
     - Dynamically accessible surface of the NZ atom, i.e. the number of positions it can take
       within the structure.
   * - ``aev``
     - :py:class:`~resdy.features.aev.AEV`
     - ANI-2x atomic environment vector of the NZ atom, a 1008-element description of its local
       environment.
   * - ``legolas``
     - :py:class:`~resdy.features.nmr.NMR`
     - 15N NMR chemical shift, calculated with LEGOLAS.
   * - ``aev_legolas``
     - :py:class:`~resdy.features.nmr.NMR`
     - The ANI-2x atomic environment vectors dumped by LEGOLAS while it calculates the shifts.
   * - ``seqcharge``
     - :py:class:`~resdy.features.charge.Charge`
     - Summed charge of the amino acids on either side of the residue of interest.
   * - ``flexibility``
     - :py:class:`~resdy.features.flexibility.Flexibility`
     - Average B-factor over the atoms of the residue of interest.
   * - ``frustration``
     - :py:class:`~resdy.features.frustration.Frustration`
     - Frustration metric, calculated with Frustratometer.
   * - ``density``
     - :py:class:`~resdy.features.frustration.Frustration`
     - Local density, calculated with Frustratometer alongside the frustration metric.
   * - ``melodia``
     - :py:class:`~resdy.features.structure.Structure`
     - All the Melodia structural descriptors at once.
   * - ``curvature``, ``writhing``, ``torsion``, ``arc_length``, ``phi``, ``psi``
     - :py:class:`~resdy.features.structure.Structure`
     - The individual Melodia structural descriptors, requested one by one.
   * - ``evolution``
     - :py:class:`~resdy.features.evolution.Evolution`
     - ESM language model embedding of the sequence around the residue of interest.
   * - ``rmsf``
     - :py:class:`~resdy.features.rmsf.RMSF`
     - Root mean square fluctuation of the residue across all the structures curated for its
       Uniprot entry.

``'all'`` may be given instead of a list of names, as a shorthand for the preset list defined in
:class:`Measure <resdy.measure.Measure>`.

.. toctree::
   :titlesonly:
   :hidden:

{% for page in pages|selectattr("type", "equalto", "module")|sort(attribute="name") %}
{% if not page.name.endswith("error_reporting") %}
   {{ page.include_path }}
{% endif %}
{% endfor %}
