from multiprocessing import cpu_count
from multiprocessing import Manager
from multiprocessing.pool import Pool
from contextlib import redirect_stdout
from ast import literal_eval
import pandas as pd
import numpy as np
import biobox as bb
from features.error_reporting import report_error_to_file


# AEV packages
try:
    from ase import Atoms
    import torch
    import torchani
    aev_packages_available = True
except Exception as e:
    aev_packages_available = False
    print(f'Packages required for AEV calculation are not available, '
          f'will not be able to calculate AEVs. Error: {e}')


class AEV():
    '''
    Class to house the different methods for calculating representations for the local structure
    of lysines within the protein structure termed atomic environment vectors (AEVs). Potential AEV
    representations are currently:
    1. ANI-2x AEVs
    2. LEGOLAS ANI-2x AEVs - currently calculated through the LEGOLAS nmr package as an add-on
    3. Coarse-grain representation AEVs 
    '''

    def __init__(self, include_modified=False, aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                                                'modified_codes': ['LYE', 'KCX'],
                                                                'atom_select_names_nonmod': ['NZ'],
                                                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the AEV class, provides general global variables and information taken forward from
        the overall measures class in here.

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

        if not aev_packages_available:
            raise ImportError('>> Packages required for AEV calculations (ase/torch/torchani) are '
                              'not available, aev will be removed from features to calculate.')

        # Preparation of AEV computer
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.ANI = torchani.models.ANI2x(periodic_table_index=True).to(device=self.device)


    def calculate_aevs(self, path):
        '''
        Calculate the Atomic Environment Vectors (AEVs) of the NZ atom within the lysine structure

        Method
        ------
        Uses the ANI-2x AEV calculator to calculate the AEVs
        Option available to use cuaev accelerated AEV calculation, can also just be run with a cpu
        For each NZ atom withing the lysines of the protein, a substructure is created 
            including all atoms within a cutoff distance
        The cutoff distance is set at 6A currently as this was the minimum distance needed
            for all information and agrees with pkaANI cutoff set
        The AEV is a vector with length 1008 representing the environment for the lysine

        Parameters
        ----------
        path : string
            The path of the pdb file that SASA is being calculated for.

        Returns
        -------
        df_aevs : dataframe
            Dataframe with information on chain, residue number and AEV output. Outline:
            Chain   Resid   aev
            x       x       [x]

        Example
        -------
        >> print(calculate_aevs(1ubq.pdb))
        Chain Resid                                                aev
        0     A     6  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        1     A    11  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        2     A    27  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        3     A    29  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        4     A    33  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        5     A    48  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        6     A    63  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        '''

        df_aevs = pd.DataFrame(columns=["Chain", "Resid", "aev"])

        # 1: prepare the biobox structure, take the species and coordinates and convert to Atoms structure, find the locations of the NZ atoms within the lysines in the strucure
        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)

            if self.include_modified:
                coords_nz, idx_nz = M.atomselect('*',
                                                 (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']),
                                                 self.aa_properties['atom_select_names_modified'],
                                                 use_resname=True, get_index=True)
                # due to wider selection criteria, possible to get more than 1 hit per residue of interest, remove duplicates
                key_res_chain = zip(list(M.data['resid'].values[idx_nz]), list(M.data['chain'].values[idx_nz]))
                pairs_seen, keep_pos = set(), []
                for pair, pos in zip(key_res_chain, range(len(idx_nz))):
                    if pair not in pairs_seen:
                        pairs_seen.add(pair)
                        keep_pos.append(pos)
                idx_nz = idx_nz[keep_pos]

            else:
                coords_nz, idx_nz = M.atomselect('*',
                                                 self.aa_properties['non_modified_codes'],
                                                 self.aa_properties['atom_select_names_nonmod'],
                                                 use_resname=True, get_index=True)

            all_coords, idx = M.atomselect('*','*','*', get_index=True)
            list_resids = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.data['resname'][idx_nz]))
        except Exception as e:
            if self.record_errors: report_error_to_file('AEV 1', path, str(e), self.error_filename)
            print(f'AEV Calculations: 1 - could not create atomic structure representation: {e}')
            return

        # 2: Iterate over the protein structure to cut out substructures and calculate an AEV at each of these.
        try:
            for lys_idx, lys_coord in enumerate(coords_nz):
                # 2.1: for the NZ atom of the lysine, find all the atoms within the cutoff distance and create a substructure
                distance_cut_off = 6  # current cutoff for substructure from analysis done on different cutoffs and matching pkaANI
                coords_euc_dists = np.linalg.norm(all_coords - lys_coord, axis=1)
                list_close_points = np.where(coords_euc_dists < distance_cut_off)[0]

                S = M.get_subset(idxs=list_close_points)
                chain = list_chains[lys_idx]
                resid = list_resids[lys_idx]
                temp_atom_species = S.data['atomtype']
                temp_coords = S.coordinates[0]
                temp_structure = Atoms(temp_atom_species, temp_coords)
                aevs = None

                # 2.2: calculate the AEV for the subset of the protein and add this to the output dataframe
                try:
                    species = self.ANI.species_to_tensor(temp_structure.get_chemical_symbols()).unsqueeze(0).to(device=self.device)
                    ani_coords = torch.tensor(temp_structure.get_positions(), dtype=torch.float32).unsqueeze(0).to(device=self.device)
                    aevs = self.ANI.aev_computer((species, ani_coords)).aevs
                    # match up the position of the lysine of interest to inside the structure cutout
                    lys_nz_subloc = np.where(list_close_points == idx_nz[lys_idx])[0]
                    if len(lys_nz_subloc) != 1:
                        raise ValueError(f'>> Could not find the specified lysine within its own structure section for {idx_nz[lys_idx]}')
                    aevs = str(aevs[0,int(lys_nz_subloc[0]),:].detach().cpu().tolist())
                except Exception as e:
                    if self.record_errors: report_error_to_file('AEV 1.1', path, str(e), self.error_filename)
                    print(f'AEV Calculations: could not create AEV for resid {idx_nz[lys_idx]} of protein {path}, error: {e}')

                # 2.3: Append the new AEV to the output dataframe
                if self.include_modified: aev_to_append = {'Chain': chain, 'Resid': resid, 'aev': aevs, 'Modified': list_modified[lys_idx]}
                else: aev_to_append = {'Chain': chain, 'Resid': resid, 'aev': aevs}
                df_aevs = pd.concat([df_aevs, pd.DataFrame([aev_to_append])], ignore_index=True)

        except torch.cuda.OutOfMemoryError:
            # potential that calculating the AEVs could overload the gpu, if too much memory, catch this and skip the file
            print(f'AEV calc error: CUDA memory error with file: {path}, skipping')
            if self.record_errors: report_error_to_file('AEV 2', path, 'CUDA memory error with file', self.error_filename)
            return df_aevs
        except MemoryError:
            # potential that calculating the AEVs could overload the cpu, if too much memory, catch this and skip the file
            print(f'AEV calc error: CPU memory error with file: {path}, skipping')
            if self.record_errors: report_error_to_file('AEV 2', path, 'CPU memory error with file', self.error_filename)
            return df_aevs
        except Exception as e:
            print(f'AEV Calculations: 2 - could not create the AEVs for the protein for protein {path}, error: {e}')
            if self.record_errors: report_error_to_file('AEV 2', path, str(e), self.error_filename)
            return df_aevs

        # 4: if everything has worked, return the dataframe with the AEVs for the protein
        return df_aevs
