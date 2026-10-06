import os
import pandas as pd
import numpy as np
import biobox as bb
from .error_reporting import report_error_to_file


#: Radial cutoff of ANI-2x, in Angstrom, from the name of its parameter file
#: (``rHCNOSFCl-5.1R_16-3.5A_a8-4.params``). Used only if the loaded model does not expose it.
ANI2X_RADIAL_CUTOFF = 5.1


class AEV():
    '''
    Representations of the local structure of lysines within the protein structure, termed atomic
    environment vectors (AEVs). This class looks into calculating these using ANI-2x
    '''

    def __init__(self, include_modified=False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                  'modified_codes': ['LYE', 'KCX'],
                                  'atom_select_names_nonmod': ['NZ'],
                                  'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt', cutoff=None):
        '''
        Initialise the AEV class, provides general global variables and information taken forward
        from the overall measures class in here.

        :param include_modified: Toggle to include residues which have been modified within the
            featurisation
        :type include_modified: bool
        :param aa_properties: Properties of the amino acid of interest to investigate modification
            sites for. Defaults to lysine for carbamylation. Properties are the 3 letter codes for
            non modified ('non_modified_codes') and modified ('modified_codes') and the atom names
            for non modified ('atom_select_names_nonmod') and modified
            ('atom_select_names_modified')
        :type aa_properties: dict
        :param error_filename: Name of the text file passed through from overall measures to write
            any errors from calculating features out to.
        :type error_filename: str
        :param cutoff: Radius, in Angstrom, of the substructure cut out around the anchor atom
            before its AEV is computed. Defaults to None, meaning the radial cutoff of the
            ANI-2x model itself (5.1 A). An atom further than that from the anchor contributes
            nothing to its AEV, so a larger value gives the same result, while a smaller one
            truncates the environment the AEV describes.
        :type cutoff: float
        '''
        if cutoff is not None and not cutoff > 0:
            raise ValueError(f'AEV cutoff must be positive, got {cutoff}')
        self.cutoff = cutoff

        try:
            from ase import Atoms
            import torch
            from torchani.models import ANI2x
            from torchani.aev import AEVComputer
            from torchani.utils import ChemicalSymbolsToInts

            self.Atoms = Atoms
            self.torch = torch
            self.ANI2x = ANI2x
            self.AEVComputer = AEVComputer
            self.ChemicalSymbolsToInts = ChemicalSymbolsToInts
        except ImportError as e:
            raise ImportError(f'>> Packages required for AEV calculations (ase/torch/torchani) are '
                              f'not available, aev will be removed from features to calculate. '
                              f'Error: {str(e)}') from e

        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        if self.error_filename != 'no_record': self.record_errors = True
        else: self.record_errors = False

        self.device = None
        self.ANI_model = None
        self.species_converter = None


    def _ensure_aev_model(self):
        '''
        On the first use of the model, load the model in the desired process for the work
        '''
        if self.ANI_model is not None:
            return
        self.device = self.torch.device('cuda' if self.torch.cuda.is_available() else 'cpu')
        self.ANI_model = self.ANI2x(periodic_table_index=True).to(device=self.device)
        self.species_converter = self.ChemicalSymbolsToInts(['H', 'C', 'N', 'O', 'S', 'F', 'Cl'])

        # where the radial cutoff lives depends on the torchani version: aev_computer.radial
        # .cutoff in 2.9, aev_computer.Rcr in earlier releases
        aev_computer = self.ANI_model.aev_computer
        model_cutoff = ANI2X_RADIAL_CUTOFF
        for value in (getattr(getattr(aev_computer, 'radial', None), 'cutoff', None),
                      getattr(aev_computer, 'Rcr', None)):
            if value is not None:
                model_cutoff = float(value)
                break
        if self.cutoff is None:
            self.distance_cut_off = model_cutoff
        else:
            self.distance_cut_off = float(self.cutoff)
            if self.distance_cut_off < model_cutoff:
                print(f'>> AEV cutoff of {self.distance_cut_off} A is below the radial cutoff of '
                      f'the ANI-2x model ({model_cutoff} A), so the AEVs describe a truncated '
                      f'environment')


    def calculate(self, path):
        '''
        Calculate the Atomic Environment Vectors (AEVs) of the NZ atom within the lysine structure

        .. rubric:: Method

        - Uses the ANI-2x AEV calculator to calculate the AEVs
        - Option available to use cuaev accelerated AEV calculation, can also just be run with a cpu
        - For each NZ atom within the lysines of the protein, a substructure is created including
          all atoms within a cutoff distance
        - The cutoff distance is the radial cutoff of ANI-2x (5.1 A) unless set otherwise at
          construction: atoms further from the anchor contribute nothing to its AEV
        - The AEV is a vector with length 1008 representing the environment for the lysine

        :param path: The path of the pdb file that the AEVs are being calculated for.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number and AEV output. Outline::

                Chain   Resid   aev
                x       x       [x]
        :rtype: pandas.DataFrame
        '''

        self._ensure_aev_model()

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
                # the loop below walks coords_nz while indexing the deduplicated lists, so the
                # coordinates have to be filtered alongside the indices
                coords_nz = coords_nz[keep_pos]

            else:
                coords_nz, idx_nz = M.atomselect('*',
                                                 self.aa_properties['non_modified_codes'],
                                                 self.aa_properties['atom_select_names_nonmod'],
                                                 use_resname=True, get_index=True)

            all_coords, _ = M.atomselect('*','*','*', get_index=True)
            list_resids = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.data['resname'][idx_nz]))
        except Exception as e:
            raise Exception(f'AEV Calculations: 1 - could not create atomic structure representation for file {path}: {e}')

        # 2: Iterate over the protein structure to cut out substructures and calculate an AEV at each of these.
        try:
            for lys_idx, lys_coord in enumerate(coords_nz):
                # 2.1: for the NZ atom of the lysine, find all the atoms within the cutoff distance and create a substructure
                coords_euc_dists = np.linalg.norm(all_coords - lys_coord, axis=1)
                list_close_points = np.where(coords_euc_dists < self.distance_cut_off)[0]

                S = M.get_subset(indices=list_close_points)
                chain = list_chains[lys_idx]
                resid = list_resids[lys_idx]
                temp_atom_species = S.data['atomtype']
                temp_coords = S.coordinates[0]
                temp_structure = self.Atoms(temp_atom_species, temp_coords)
                aevs = None

                # 2.2: calculate the AEV for the subset of the protein and add this to the output dataframe
                try:
                    species = self.species_converter(temp_structure.get_chemical_symbols()).unsqueeze(0).to(device=self.device)
                    ani_coords = self.torch.tensor(temp_structure.get_positions(), dtype=self.torch.float32).unsqueeze(0).to(device=self.device)
                    aevs = self.ANI_model.aev_computer(species, ani_coords)

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

        except self.torch.cuda.OutOfMemoryError:
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


if __name__ == '__main__':
    aev = AEV(include_modified=False)
    print(aev.calculate(f'demo{os.sep}curated{os.sep}1A6M-alt1A.pdb'))
