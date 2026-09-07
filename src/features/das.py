import os
import pandas as pd
import biobox as bb
from src.features.error_reporting import report_error_to_file

class DAS():
    '''
    Class to house the different methods for calculating dynamically accessible surface
    area (das) values for structures
    '''

    def __init__(self, include_modified=False, aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                                                'modified_codes': ['LYE', 'KCX'],
                                                                'atom_select_names_nonmod': ['NZ'],
                                                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the DAS class, include any global variables that are required from
        measures in here.

        Parameters
        ----------
        include_modified : bool
            Toggle to include residues which have been modified within the featurisation
        aa_properties -> dict
            Properties of the amino acid of interest to investigate modification sites for.
            Defaults to lysine for carbamylation. Properties are the 3 letter codes for
            non modified ('non_modified_codes') and modified ('modified_codes') and the atom
            names for non modified ('atom_select_names_nonmod') and modified ('atom_select_names_modified')
        error_filename : str
            Name of the text file passed through from overall measures to write any errors from
            calculating features out to.
        '''
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        if self.error_filename != 'no_record': self.record_errors = True
        else: self.record_errors = False

    def calculate(self, path):
        '''
        Calculate the Dynamically Accessible Surface (DAS) of the NZ atom in the lysine structure
        This is effectively the number of positions that the NZ atom can take within the structure
        of the protein.

        Method
        ------
        Uses biobox functionality to calculate the value
        Create a molecule for the protein structure from the bb.Molecule class
        Use the bb.Xlink class to setup the linking module
        Use the hidden method .__get_half_sphere() to work out the das value
        As the density of points in the sphere of the NZ atom of the lysine is constant,
            the das value is the number of points that are accessible


        Parameters
        ----------
        path : string
            The path of the pdb file that DAS is being calculated for.

        Returns
        -------
        df_das : dataframe
            Dataframe with information on chain, residue number and DAS output. Outline:
            Chain   Resid   das
            x       x       x

        Example
        -------
        >> print(calculate_das(1ubq.pdb))
        Chain  Resid  das
        0     A      6   36
        1     A     11   47
        2     A     27   22
        3     A     29   37
        4     A     33   46
        5     A     48   34
        6     A     63   30
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
            print(f'DAS Calculation: 1 - could not load and identify the NZ atoms within the lysines of the structure: {e}')

        # 2: Setup the Xlink module and create the half spheres
        try:
            XL = bb.Xlink(M)
            das_output = []
        except Exception as e:
            if self.record_errors: report_error_to_file('DAS 2', path, str(e), self.error_filename)
            print(f'DAS Calculation: 2 - Failed to setup the Xlink biobox class: {e}')

        for i, lys_nz_idx in enumerate(idx_nz):
            try:
                # the parameteres (pts_surf, thresh, radii) for the _get_half_sphere are already set
                # for lysine residues therefore the only parameter that needs to be set is i: this
                # is the index of the atom of interest within the lysine
                half_sphere_coords = XL._get_half_sphere(i=lys_nz_idx)
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
