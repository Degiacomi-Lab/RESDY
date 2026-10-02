##########################
Frequently Asked Questions
##########################


How do I add a new feature?
---------------------------

Write a class with a ``calculate`` method and drop the file into ``resdy/features/``. The
folder is scanned when the package is imported, so nothing else needs editing: no registry,
no list of names, no change to :class:`Measure <resdy.measure.Measure>`.

Two naming rules apply. The class name must be the name you will ask for in
``features_dict``, upper-cased, because that is how the feature is looked up. The file name
is free, so ``my_neighbour_count.py`` holding ``class NEIGHBOURS`` is picked up by
``features_dict={'neighbours': {}}``. Naming the file after the class is a sensible habit,
not a requirement.

``calculate`` receives the path of one curated PDB file and returns a
:class:`pandas.DataFrame` with a ``Chain`` column, a ``Resid`` column, and one column named
exactly as the feature was requested. Those three columns are the whole contract: the values
are merged into the output table on chain and residue number.

The constructor is given ``include_modified``, ``aa_properties`` and ``error_filename`` by
``Measure``, so it has to accept all three. ``aa_properties`` is what tells the feature which
residue to look at and at which atom, and is the only thing standing between a feature and a
hard-coded lysine.

.. code-block:: python

    import numpy as np
    import pandas as pd
    import biobox as bb
    from .error_reporting import report_error_to_file


    class NEIGHBOURS():
        '''
        Number of heavy atoms within a cutoff of the residue of interest.
        '''

        def __init__(self, include_modified=False,
                     aa_properties=None,
                     error_filename='measure_errors.txt',
                     cutoff=8.0):
            self.include_modified = include_modified
            self.aa_properties = aa_properties
            self.error_filename = error_filename
            self.record_errors = error_filename != 'no_record'
            self.cutoff = cutoff

        def calculate(self, path):
            try:
                M = bb.Molecule()
                M.import_pdb(path, include_hetatm=True)

                # the anchor atom of every residue of interest, one index per residue
                idx = M.atomselect('*',
                                   self.aa_properties['non_modified_codes'],
                                   self.aa_properties['atom_select_names_nonmod'],
                                   use_resname=True, get_index=True)[1]
                all_coords = M.atomselect('*', '*', '*', get_index=True)[0]

                counts = [int(np.sum(np.linalg.norm(all_coords - all_coords[i], axis=1)
                                     < self.cutoff)) - 1 for i in idx]

                return pd.DataFrame({'Chain': list(M.data['chain'][idx]),
                                     'Resid': list(M.data['resid'][idx]),
                                     'neighbours': counts})

            except Exception as e:
                if self.record_errors:
                    report_error_to_file('NEIGHBOURS', path, str(e), self.error_filename)
                print(f'NEIGHBOURS calculation failed: {e}')
                return pd.DataFrame(columns=['Chain', 'Resid', 'neighbours'])

Saved as ``resdy/features/NEIGHBOURS.py``, it is then used like any shipped feature::

    import resdy as RD

    M = RD.Measure(df_input=df, features_dict={'neighbours': {}, 'sasa': {}})
    M.measure_data()

``resdy/features/FEATURE.py`` is a fuller template, including the block that handles
``include_modified``.


How do I pass parameters to a feature?
--------------------------------------

The value of each entry of ``features_dict`` is passed to the feature's constructor as
keyword arguments, so any parameter with a default becomes user-settable::

    features_dict={'neighbours': {'cutoff': 12.0},
                   'depth': {'calculation_type': 'AtomDepth'},
                   'sasa': {}}

An empty dictionary means "use the defaults".

``include_modified``, ``aa_properties`` and ``error_filename`` are supplied by
``Measure`` and must not be listed: passing one of them raises a duplicate-argument error.
``outdir`` is supplied to ``propka`` and ``legolas``, and ``df_proteins`` to ``rmsf``, unless
you give them yourself.

.. note::
   A parameter name that does not match the constructor removes the feature from the run
   rather than stopping it. The message ``Failed to add measure feature ...`` is printed and
   the output simply has no such column, which is easy to miss in a long log. Check the names
   against the class if a column you asked for is absent.


How do I say which residues my feature can handle?
--------------------------------------------------

Three optional class attributes let a feature declare this itself, so that no central table
has to be edited when a feature is added:

``SUPPORTED_RESIDUES``
    Set of three-letter codes the feature can act on. When absent, the feature is assumed to
    work on any residue. ``Measure`` drops the feature, with a message, when the residue of
    interest is not in the set.

``RESIDUE_KWARGS``
    Constructor arguments that depend on the residue, keyed on three-letter code. They are
    applied underneath ``features_dict``, so a value the user asks for always wins.

``UNSUPPORTED_REASON``
    A string, or a dictionary keyed on three-letter code, quoted back to the user when the
    feature is dropped.

.. code-block:: python

    class NEIGHBOURS():
        SUPPORTED_RESIDUES = {'LYS', 'ARG', 'HIS'}
        RESIDUE_KWARGS = {'HIS': {'cutoff': 6.0}}
        UNSUPPORTED_REASON = 'the shell was calibrated on basic side chains'

:class:`DAS <resdy.features.DAS.DAS>` is the worked example among the shipped features: it
carries its half-sphere radii per residue and explains, residue by residue, why it cannot act
on the six it excludes.


What are the current limitations on a new feature?
--------------------------------------------------

*One column per feature.* Only the column named after the feature is merged into the output
table, and any other column of the returned dataframe is discarded without a message.
Features producing several numbers, such as melodia and frustration, are handled by special
cases inside ``Measure`` and cannot yet be reproduced by a dropped-in file.

*The folder, not an arbitrary path.* Features are discovered in ``resdy/features/`` only, so
with RESDY installed rather than run from a clone this means writing into the installed
package.

*Modified residues.* When ``include_modified=True``, the returned dataframe must also carry a
``Modified`` boolean column, or merging the result raises a :class:`KeyError`.
