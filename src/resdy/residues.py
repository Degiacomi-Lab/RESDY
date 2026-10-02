'''
Residue presets: which residue codes and which atom the pipeline works on.

Every feature computes its value at one *anchor atom* per residue of interest. The anchor
is not merely a label: SASA is summed over it, DAS grows its half sphere from it, AEV
centres its substructure on it, and depth is measured to it. It is therefore the reactive
atom of the side chain where there is one (NZ for lysine, SG for cysteine), the outermost
well defined side-chain atom otherwise, and CB where the side chain has no reactive group
at all (CA for glycine, which has nothing else).

This module holds the presets for the twenty standard amino acids and a single resolver that
validates presets and user dictionaries by the same rule::

    import resdy as RD
    M = RD.Measure(df_input=df, residue_of_interest='TYR')

:data:`AA_PRESETS` is a plain dictionary and can be read, listed or extended by the user.

It holds residue identity only, and deliberately knows nothing about any individual feature.
Whether a feature can act on a residue, and what it needs in order to do so, is declared by
the feature class itself through three optional class attributes, so that a feature added by
dropping a file into ``features/`` needs no edit here:

``SUPPORTED_RESIDUES``
    The set of three-letter codes the feature can act on, or None (the default) for all of
    them. :class:`Measure <resdy.measure.Measure>` drops the feature for any other residue.
``RESIDUE_KWARGS``
    Per-residue constructor arguments, keyed on three-letter code, merged underneath whatever
    the user passed in ``features_dict``. :class:`DAS <resdy.features.DAS.DAS>` uses this to
    carry its half-sphere radii.
``UNSUPPORTED_REASON``
    A string, or a dictionary keyed on three-letter code, quoted back to the user when the
    feature is dropped.
'''

import copy


#: The four keys an ``aa_properties`` dictionary is required to carry.
AA_PROPERTY_KEYS = ('non_modified_codes', 'modified_codes',
                    'atom_select_names_nonmod', 'atom_select_names_modified')

#: Presets for the twenty standard amino acids, keyed on the three-letter code.
#:
#: ``non_modified_codes`` lists the unmodified forms, the canonical three-letter code
#: **first** (:class:`PROPKA <resdy.features.PROPKA.PROPKA>` and
#: :class:`PKAANI <resdy.features.PKAANI.PKAANI>` parse their output with a regular
#: expression built from that first entry alone), followed by the protonation-state
#: aliases that force fields and MD engines write. ``modified_codes`` lists common
#: modified forms whose CCD parent is the residue in question, so each is a genuine
#: modified form rather than a lookalike. The lists are representative rather than
#: exhaustive: the CCD holds several hundred cysteine and lysine derivatives, and these
#: are a starting point to be extended as particular chemistry is studied.
#:
#: A modified component does not always reuse its parent's atom names, so three residues
#: name a second anchor: citrulline (CIR) is numbered from the backbone outwards and its
#: ureido carbon is C7, gamma-carboxyglutamate (CGU) carries CD1 and CD2 and has no CD,
#: and selenomethionine (MSE) has SE in place of SD.
AA_PRESETS = {
    'ALA': {'non_modified_codes': ['ALA'],
            'modified_codes': [],
            'atom_select_names_nonmod': ['CB'],
            'atom_select_names_modified': ['CB']},
    'ARG': {'non_modified_codes': ['ARG'],
            'modified_codes': ['CIR', '2MR', 'AGM'],
            'atom_select_names_nonmod': ['CZ'],
            'atom_select_names_modified': ['CZ', 'C7']},
    'ASN': {'non_modified_codes': ['ASN'],
            'modified_codes': ['MEN'],
            'atom_select_names_nonmod': ['ND2'],
            'atom_select_names_modified': ['ND2']},
    'ASP': {'non_modified_codes': ['ASP', 'ASH'],
            'modified_codes': [],
            'atom_select_names_nonmod': ['CG'],
            'atom_select_names_modified': ['CG']},
    'CYS': {'non_modified_codes': ['CYS', 'CYX'],
            'modified_codes': ['CSO', 'CSD', 'OCS', 'CSX', 'CME',
                               'CSS', 'SNC', 'SMC', 'CAS', 'CAF'],
            'atom_select_names_nonmod': ['SG'],
            'atom_select_names_modified': ['SG']},
    'GLN': {'non_modified_codes': ['GLN'],
            'modified_codes': ['MEQ'],
            'atom_select_names_nonmod': ['NE2'],
            'atom_select_names_modified': ['NE2']},
    'GLU': {'non_modified_codes': ['GLU', 'GLH'],
            'modified_codes': ['CGU'],
            'atom_select_names_nonmod': ['CD'],
            'atom_select_names_modified': ['CD', 'CD1']},
    'GLY': {'non_modified_codes': ['GLY'],
            'modified_codes': [],
            'atom_select_names_nonmod': ['CA'],
            'atom_select_names_modified': ['CA']},
    'HIS': {'non_modified_codes': ['HIS', 'HID', 'HIE', 'HIP'],
            'modified_codes': ['HIC', 'NEP'],
            'atom_select_names_nonmod': ['NE2'],
            'atom_select_names_modified': ['NE2']},
    'ILE': {'non_modified_codes': ['ILE'],
            'modified_codes': [],
            'atom_select_names_nonmod': ['CD1'],
            'atom_select_names_modified': ['CD1']},
    'LEU': {'non_modified_codes': ['LEU'],
            'modified_codes': ['NLE'],
            'atom_select_names_nonmod': ['CG'],
            'atom_select_names_modified': ['CG']},
    'LYS': {'non_modified_codes': ['LYS', 'LYSN', 'LYN'],
            'modified_codes': ['LYE', 'KCX', 'MLY', 'MLZ', 'M3L', 'ALY', 'LYR'],
            'atom_select_names_nonmod': ['NZ'],
            'atom_select_names_modified': ['NZ', 'N07']},
    'MET': {'non_modified_codes': ['MET'],
            'modified_codes': ['MSE', 'OMT', 'MHO'],
            'atom_select_names_nonmod': ['SD'],
            'atom_select_names_modified': ['SD', 'SE']},
    'PHE': {'non_modified_codes': ['PHE'],
            'modified_codes': ['MEA'],
            'atom_select_names_nonmod': ['CZ'],
            'atom_select_names_modified': ['CZ']},
    'PRO': {'non_modified_codes': ['PRO'],
            'modified_codes': ['HYP'],
            'atom_select_names_nonmod': ['CG'],
            'atom_select_names_modified': ['CG']},
    'SER': {'non_modified_codes': ['SER'],
            'modified_codes': ['SEP'],
            'atom_select_names_nonmod': ['OG'],
            'atom_select_names_modified': ['OG']},
    'THR': {'non_modified_codes': ['THR'],
            'modified_codes': ['TPO'],
            'atom_select_names_nonmod': ['OG1'],
            'atom_select_names_modified': ['OG1']},
    'TRP': {'non_modified_codes': ['TRP'],
            'modified_codes': ['TRO', 'TRQ', '0AF'],
            'atom_select_names_nonmod': ['NE1'],
            'atom_select_names_modified': ['NE1']},
    'TYR': {'non_modified_codes': ['TYR'],
            'modified_codes': ['PTR', 'TYS', 'NIY', 'DAH', 'TYY'],
            'atom_select_names_nonmod': ['OH'],
            'atom_select_names_modified': ['OH']},
    'VAL': {'non_modified_codes': ['VAL'],
            'modified_codes': [],
            'atom_select_names_nonmod': ['CB'],
            'atom_select_names_modified': ['CB']},
}

def _as_list(value, field):
    '''
    Normalise one field of an ``aa_properties`` dictionary to a list of strings.

    An empty string or None becomes an empty list, and a bare string becomes a
    one-element list, so that a user dictionary written either way behaves the same.

    :param value: The value given for the field.
    :param field: Name of the field, used in the error message.
    :type field: str
    :returns: The normalised list.
    :rtype: list
    :raises TypeError: if the value is neither a string, None, nor a sequence of strings.
    '''
    if value is None or (isinstance(value, str) and value == ''):
        return []
    if isinstance(value, str):
        return [value]
    try:
        out = list(value)
    except TypeError:
        raise TypeError(f'aa_properties["{field}"] should be a list of three-letter codes '
                        f'or atom names, got {type(value).__name__}')
    if not all(isinstance(v, str) for v in out):
        raise TypeError(f'aa_properties["{field}"] should contain only strings, got {out}')
    return out


def resolve_residue(spec):
    '''
    Turn a preset name or a user dictionary into a validated ``aa_properties`` dictionary.

    .. rubric:: Method

    - A string is matched against :data:`AA_PRESETS` case-insensitively.
    - A dictionary is checked by key *set*, so the four keys may be given in any order,
      and each field is normalised with :func:`_as_list`.
    - The result is a deep copy, so that a caller cannot alter the shipped presets.

    :param spec: Three-letter code of one of the twenty standard amino acids, or a full
        dictionary with the keys listed in :data:`AA_PROPERTY_KEYS`.
    :type spec: str, dict
    :returns: The validated properties of the residue of interest.
    :rtype: dict
    :raises KeyError: if a preset name is unknown, or a dictionary key is missing or extra.
    :raises ValueError: if the residue codes or the unmodified atom names are empty, or if
        the argument is neither a string nor a dictionary.
    '''
    if isinstance(spec, str):
        key = spec.strip().upper()
        if key not in AA_PRESETS:
            raise KeyError(f'Residue of interest given not known: {spec}. Available presets '
                           f'are {", ".join(sorted(AA_PRESETS))}; alternatively pass a '
                           f'dictionary with the keys {", ".join(AA_PROPERTY_KEYS)}.')
        return copy.deepcopy(AA_PRESETS[key])

    if isinstance(spec, dict):
        given = set(spec)
        wanted = set(AA_PROPERTY_KEYS)
        if given != wanted:
            missing = sorted(wanted - given)
            extra = sorted(given - wanted)
            raise KeyError(f'aa_properties dictionary should hold exactly the keys '
                           f'{", ".join(AA_PROPERTY_KEYS)} (in any order)'
                           + (f'; missing: {", ".join(missing)}' if missing else '')
                           + (f'; unexpected: {", ".join(extra)}' if extra else '')
                           + '. A key with nothing to declare can be set to [].')
        out = {k: _as_list(spec[k], k) for k in AA_PROPERTY_KEYS}
        if not out['non_modified_codes']:
            raise ValueError('aa_properties["non_modified_codes"] cannot be empty: it names '
                             'the residue to measure.')
        if not out['atom_select_names_nonmod']:
            raise ValueError('aa_properties["atom_select_names_nonmod"] cannot be empty: it '
                             'names the atom every feature is computed at.')
        return out

    raise ValueError(f'Unknown option given to residue_of_interest parameter: {spec}; please '
                     f'enter either a three-letter code or a dictionary. See the class '
                     f'documentation.')


def residue_key(aa_properties):
    '''
    The three-letter code a set of properties describes.

    This is the key a feature class is asked about when it declares ``SUPPORTED_RESIDUES``,
    ``RESIDUE_KWARGS`` or ``UNSUPPORTED_REASON``. The canonical code is taken as the first
    entry of ``non_modified_codes``, which is the same entry the pKa parsers build their
    regular expression from.

    :param aa_properties: A dictionary as returned by :func:`resolve_residue`.
    :type aa_properties: dict
    :returns: The upper-case three-letter code.
    :rtype: str
    '''
    return aa_properties['non_modified_codes'][0].upper()
