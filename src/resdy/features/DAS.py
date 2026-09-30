import os
import pandas as pd
import biobox as bb
from .error_reporting import report_error_to_file

class DAS():
    '''
    Dynamically accessible surface area (DAS) values for structures.

    biobox builds the half sphere as concentric shells around CA and then joins the anchor atom
    to the shell points it can reach, so the shells have to span the range of positions the
    anchor can actually occupy. biobox's own defaults, ``[6.3, 5.9, 5.4, 4.8]``, are lysine
    lengths, and a short side chain whose anchor sits 1.5 to 2.6 A from CA cannot reach them, so
    the feature returns 1, the anchor itself, whatever the environment. :data:`RESIDUE_KWARGS`
    therefore carries a set of shells per residue.

    Three of the supported residues are marginal rather than comfortable: with their shells,
    43% of ASP, 41% of LEU and 23% of CYS sites still return the anchor alone. Their values are
    usable as a ranking within the residue but are not comparable with those of a long side
    chain, which is true across residues in any case, since biobox derives its point density
    from the outermost shell.
    '''

    #: Radii of the concentric shells, in Angstrom, for each residue the feature can act on.
    #:
    #: Derived as the 98th percentile anchor-to-CA distance of the residue, measured over the 19
    #: relaxed structures of ``demo/curated``, scaled by the ratios of the lysine shells to their
    #: outermost value. LYS keeps biobox's literal list, so lysine output is unchanged.
    RESIDUE_KWARGS = {
        'ARG': {'radii': [6.3, 5.9, 5.4, 4.8]},
        'ASN': {'radii': [3.8, 3.5, 3.2, 2.9]},
        'ASP': {'radii': [2.6, 2.5, 2.3, 2.0]},
        'CYS': {'radii': [2.8, 2.7, 2.4, 2.2]},
        'GLN': {'radii': [5.0, 4.7, 4.3, 3.8]},
        'GLU': {'radii': [4.0, 3.7, 3.4, 3.0]},
        'HIS': {'radii': [4.7, 4.4, 4.1, 3.6]},
        'ILE': {'radii': [4.0, 3.7, 3.4, 3.0]},
        'LEU': {'radii': [2.7, 2.5, 2.3, 2.0]},
        'LYS': {'radii': [6.3, 5.9, 5.4, 4.8]},
        'MET': {'radii': [4.2, 4.0, 3.6, 3.2]},
        'PHE': {'radii': [5.3, 5.0, 4.6, 4.1]},
        'TRP': {'radii': [4.7, 4.4, 4.0, 3.6]},
        'TYR': {'radii': [6.6, 6.2, 5.7, 5.1]},
    }

    #: Residues the feature can act on. A residue is left out when the measurement showed the
    #: feature has nothing to say about it: DAS was computed for all twenty residues over the 19
    #: relaxed structures of ``demo/curated``, and a residue is excluded when the median count
    #: came back as 1, the anchor atom alone.
    SUPPORTED_RESIDUES = set(RESIDUE_KWARGS)

    #: Why the feature is unavailable, per residue, quoted back to the user when it is dropped.
    UNSUPPORTED_REASON = {
        'GLY': 'glycine has no side chain, and biobox raises on a CA anchor',
        'ALA': 'the CB anchor sits 1.5 A from CA, inside biobox\'s 2.0 A clash threshold, so '
               'every shell point is rejected (100% of 376 residues returned the anchor alone)',
        'VAL': 'the CB anchor sits 1.5 A from CA, inside biobox\'s 2.0 A clash threshold, so '
               'every shell point is rejected (100% of 398 residues returned the anchor alone)',
        'PRO': 'CG is locked in the pyrrolidine ring and 2.4 A from CA (63% of 232 residues '
               'returned the anchor alone, median 1)',
        'SER': 'the OG anchor is 2.4 A from CA and the shells barely clear the clash threshold '
               '(54% of 394 residues returned the anchor alone, median 1)',
        'THR': 'the OG1 anchor is 2.4 A from CA and the shells barely clear the clash threshold '
               '(63% of 444 residues returned the anchor alone, median 1)',
    }

    def __init__(self, include_modified=False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                  'modified_codes': ['LYE', 'KCX'],
                                  'atom_select_names_nonmod': ['NZ'],
                                  'atom_select_names_modified': ['NZ', 'N07']},
                 radii = None, error_filename = 'measure_errors.txt'):
        '''
        Initialise the DAS class, include any global variables that are required from measures in
        here.

        :param include_modified: Toggle to include residues which have been modified within the
            featurisation
        :type include_modified: bool
        :param radii: Radii of the concentric shells the half sphere is built from, in Angstrom,
            outermost first. Defaults to None, meaning biobox's own list, which holds lysine
            lengths. :class:`Measure <resdy.measure.Measure>` fills this in from
            :data:`RESIDUE_KWARGS` for the residue of interest, so it only needs setting by hand
            when studying a residue that has no entry there. At least two shells are required:
            biobox takes the largest gap between consecutive shells and raises when there is
            none.
        :type radii: list
        :param aa_properties: Properties of the amino acid of interest to investigate modification
            sites for. Defaults to lysine for carbamylation. Properties are the 3 letter codes for
            non modified ('non_modified_codes') and modified ('modified_codes') and the atom names
            for non modified ('atom_select_names_nonmod') and modified
            ('atom_select_names_modified')
        :type aa_properties: dict
        :param error_filename: Name of the text file passed through from overall measures to write
            any errors from calculating features out to.
        :type error_filename: str
        '''
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.radii = list(radii) if radii is not None else None
        if self.radii is not None and len(self.radii) < 2:
            raise ValueError(f'DAS needs at least two shell radii, got {self.radii}: biobox '
                             f'takes the largest gap between consecutive shells and raises '
                             f'when there is none.')
        self.error_filename = error_filename
        if self.error_filename != 'no_record':
            self.record_errors = True
        else:
            self.record_errors = False

    def calculate(self, path):
        '''
        Calculate the Dynamically Accessible Surface (DAS) of the anchor atom of the residue of
        interest (NZ for lysine). This is effectively the number of positions that the anchor
        atom can take within the structure of the protein.

        .. rubric:: Method

        - Uses biobox functionality to calculate the value
        - Create a molecule for the protein structure from the bb.Molecule class
        - Use the bb.Xlink class to setup the linking module
        - Use the hidden method .__get_half_sphere() to work out the das value
        - As the density of points in the sphere of the anchor atom is constant for a given set
          of radii, the das value is the number of points that are accessible

        :param path: The path of the pdb file that DAS is being calculated for.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number and DAS output. Outline::

                Chain   Resid   das
                x       x       x
        :rtype: pandas.DataFrame
        '''

        # 1: Load in the structure and locate all the NZ atoms within the lysines, calculate the list of chains and list of resids to go with this
        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)

            if self.include_modified:
                idx_nz = M.atomselect('*',
                                      (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']),
                                      self.aa_properties['atom_select_names_modified'],
                                      use_resname=True, get_index=True)[1]
                # due to wider selection criteria, possible to get more than 1 hit per residue of interest, remove duplicates
                key_res_chain = zip(list(M.data['resid'].values[idx_nz]), list(M.data['chain'].values[idx_nz]))
                pairs_seen, keep_pos = set(), []
                for pair, pos in zip(key_res_chain, range(len(idx_nz))):
                    if pair not in pairs_seen:
                        pairs_seen.add(pair)
                        keep_pos.append(pos)
                idx_nz = idx_nz[keep_pos]

            else:
                idx_nz = M.atomselect('*',
                                      self.aa_properties['non_modified_codes'],
                                      self.aa_properties['atom_select_names_nonmod'],
                                      use_resname=True, get_index=True)[1]

            lys_res_nums = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
            list_modified = list(a in  self.aa_properties['modified_codes'] for a in list(M.data['resname'][idx_nz]))
        except Exception as e:
            if self.record_errors: report_error_to_file('DAS 1', path, str(e), self.error_filename)
            print(f'DAS Calculation: 1 - could not load and identify the NZ '
                  f'atoms within the lysines of the structure: {str(e)}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'das'])

        # 2: Setup the Xlink module and create the half spheres
        try:
            XL = bb.Xlink(M)
            das_output = []
        except Exception as e:
            if self.record_errors: report_error_to_file('DAS 2', path, str(e), self.error_filename)
            print(f'DAS Calculation: 2 - Failed to setup the Xlink biobox class: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'das'])

        for i, lys_nz_idx in enumerate(idx_nz):
            try:
                # i is the index of the anchor atom of the residue of interest. The shells the
                # half sphere is built from are the ones of that residue, and biobox's own
                # (lysine) list is used when none was given
                if self.radii is None:
                    half_sphere_coords = XL._get_half_sphere(i=lys_nz_idx)
                else:
                    half_sphere_coords = XL._get_half_sphere(i=lys_nz_idx, radii=self.radii)
                # as the density of points created by the get half sphere is constant for any setup,
                # therefore can just count the number of coordinates that are returned for a measure for SASA Path
                das_output.append(len(half_sphere_coords))
            except Exception as e:
                if self.record_errors: report_error_to_file('DAS 2', path, str(e), self.error_filename)
                print(f'DAS Calculation: 2 - Failed to calculate the half spheres for the NZ atoms on lysine no {lys_res_nums[i]}: {e}')
                das_output.append(None)

        # 3: Create dataframe to return
        df_das = pd.DataFrame(columns=["Chain", "Resid", "das"])
        try:
            df_das['Chain'] = list_chains
            df_das['Resid'] = lys_res_nums
            df_das['das'] = das_output
            if self.include_modified: df_das['Modified'] = list_modified
        except Exception as e:
            if self.record_errors: report_error_to_file('DAS 3', path, str(e), self.error_filename)
            print(f'DAS Calculation: 3 - Failed to create datafame to append to the overall dataframe: {e}')

        return df_das


if __name__ == '__main__':
    das = DAS(include_modified=True)
    print(das.calculate(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
