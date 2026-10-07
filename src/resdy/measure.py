import re
import os
import io
import copy
import logging
import datetime
import glob
import time
import inspect
import json
from datetime import date
from functools import partial
import multiprocessing as mp
from multiprocessing import cpu_count
from multiprocessing import Manager
from contextlib import redirect_stdout
import pandas as pd
import numpy as np
import biobox as bb
from .helper import require_biobox
require_biobox(bb)
from .features import *
from .residues import residue_key, resolve_residue, PROTEIN_RESNAMES
from .geometry import distance_to_other_chains


def _call_with_args(func, args):
    '''
    Call ``func`` with the arguments packed in ``args``.

    ``Pool.imap`` only accepts a single-argument callable, whereas the measuring methods take
    a file and a lock. This module-level trampoline supplies the ``starmap`` behaviour while
    keeping the results ordered and streamed, which is what lets a stalled worker be detected
    (see :meth:`Measure._run_in_parallel`). It has to live at module level rather than be a
    closure or a lambda, because the non-forking start methods pickle it to reach the worker.

    :param func: The callable to run in the worker.
    :param args: Positional arguments to unpack into ``func``.
    :returns: Whatever ``func`` returns.
    '''
    return func(*args)


#: File, in the output directory, that records the settings a measurement was made with.
MEASURE_SETTINGS_FILE = 'measure_settings.json'

#: Constructor arguments of a feature that are not settings of the measurement: they are
#: filled in by Measure from its own arguments, or are data rather than parameters.
_NOT_SETTINGS = ('self', 'include_modified', 'error_filename', 'aa_properties', 'outdir',
                 'df_proteins')


def _resolved_settings(cls, kwargs):
    """
    The settings a feature class is constructed with: its defaults, overridden by the
    arguments given.

    :param cls: feature class.
    :type cls: type
    :param kwargs: keyword arguments passed to the constructor.
    :type kwargs: dict
    :returns: {argument: value}, in a form json can write.
    :rtype: dict
    """
    settings = {}
    for name, par in inspect.signature(cls.__init__).parameters.items():
        if name in _NOT_SETTINGS or par.kind in (par.VAR_POSITIONAL, par.VAR_KEYWORD):
            continue
        settings[name] = kwargs.get(name, None if par.default is par.empty else par.default)
    return json.loads(json.dumps(settings, default=str))


#: Features requested by the 'all' shorthand in ``features_dict``.
ALL_FEATURES = ['propka', 'pkaANI', 'sasa', 'depth', 'aev', 'seqcharge', 'legolas',
                'melodia', 'frustration', 'density', 'das', 'flexibility', 'evolution',
                'rmsf']


class Measure(object):
    '''
    Class to handle functions used in calling feature functions; returns a dataframe
    containing the results at the end.
    '''

    def __init__(self, df_input, outdir="result", activate_log=False, log_path='measure_log.txt',
                 features_dict={'propka': {}, 'sasa': {}, 'depth': {'calculation_type': 'ResidDepth'},
                                'aev': {}, 'das': {}, 'seqcharge': {}},
                 residue_of_interest='LYS', parallel=False, include_modified=False,
                 report_errors= True, only_relaxed=True, num_cores=0,
                 parallel_timeout=600):
        '''
        Initialisation of the Measure class. This class provides all the resources to measure
        specific quantities for the protein structures given as input

        Every residue is featurised at its anchor atom (see residue_of_interest). A residue whose
        anchor atom is absent from the structure cannot be featurised and is left out of the
        output, and the number left out is reported per file. Beside the features, every row
        carries the metadata column 'Min_Dist_Other_Chain': the distance, in A, from the anchor
        atom to the closest heavy atom of a protein residue in any other chain of the measured
        file, or NaN when the file holds a single protein chain. Ions, waters and ligands are not
        counted as a chain. The chains are those of the curated file, i.e. the deposited
        asymmetric unit. The column is computed whatever features_dict holds, and is listed in
        ``METADATA_COLUMNS`` in resdy/helper.py.

        :param df_input: The input dataframe containing information on the structures over which the
            measurements will be done. This is usually the output given from the curation steps
            (pdb.df). The format of this file depends on if running PDB_only or not. If running PDB
            only you will just need to give one column which is a list of paths to the pdb
            files. If running fully the dataframe will contain columns of 'Uniprot_Entry',
            'PDB_Code', 'Method', 'Resolution', 'Chains'
        :type df_input: pandas.DataFrame
        :param outdir: The name of the directory where the measurement output will be written to.
        :type outdir: str
        :param activate_log: By default a log is produced for the measurements, the option here
            enables a more detailed log of the measurements work for debugging.
        :type activate_log: bool
        :param log_path: The name of the output file which contains the log of the measurements.
            This file can be used to create the measurement csv file through using the
            recover_from_log() function.
        :type log_path: str
        :param features_dict: The list of measurements that you wish to use on the given structures.
            Select which of the following options to use: 'propka', 'pkaANI', 'sasa', 'depth',
            'aev', 'das', 'seqcharge', 'flexibility', 'legolas', 'melodia',
            'curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi', 'frustration',
            'density', 'evolution', 'rmsf', 'secondarystructure'. Any class defined in a script
            added to the features folder can also be requested by its class name. 'all' is a shorthand
            for the preset list defined in _setup_measures(). A dict should be given here with a
            dict per feature included eg {'depth': {}, etc...}, inside the dict per feature should
            house any optional arguments available for that specific feature class, if defaults are
            okay, leave as {}.
        :type features_dict: dict
        :param residue_of_interest: The residue of interest to calculate measurements for. The
            three-letter code of any of the twenty standard amino acids can be given as a
            string, and is matched against the presets in :data:`resdy.residues.AA_PRESETS`,
            which name the codes of the modified forms and the atom the features are computed
            at. If you are investigating a non-standard residue, or would like more control
            over a preset, please enter a dictionary of the following format (the keys may be
            given in any order):
            {'non_modified_codes': [residue codes of standard state, canonical code first],
            'modified_codes': [codes of modfified state, can be left as [] if not investigating],
            'atom_select_names_nonmod': [atom names of interest in standard state],
            'atom_select_names_modified': [atom names of interest in modified residues]}
            A feature that cannot act on the residue given is dropped from features_dict with a
            message quoting its own reason, and one that has per-residue settings is given
            them. Both are declared by the feature class, not listed here; see
            :mod:`resdy.residues`.
        :type residue_of_interest: str, dict
        :param parallel: Option to run the measurements in parallel.
        :type parallel: bool
        :param include_modified: Option to include lysines that have been seen to be modified in the
            measurements analysis. If False, only lysines of type 'LYS' will be included in the
            measurements. If True, lysines of types 'LYE' will be included in the measurements as
            well as all 'LYS' residues. In either case, a column will be included stating if the
            measured residue is a modified one.
        :type include_modified: bool
        :param report_errors: Option to record any of the protein files which are giving errors when
            measures calculations are being performed. This will write the file and the error to a
            separate text document labelled "measures_errors_{date}.txt".
        :type report_errors: bool
        :param only_relaxed: Option to only calculate measurements for structures that are relaxed
            if there is a relaxed structure available for the structure. If set to False, measures
            will be calculated to both original and relaxed form. Default is True.

            A structure with no minimised copy is measured unrelaxed either way, which happens
            whenever minimisation was not requested, was skipped, or failed. The 'Source' column
            of the output records which copy each row was taken from, 'relaxed' or 'unrelaxed',
            so that the two are not silently mixed.
        :type only_relaxed: bool
        :param num_cores: Number of cores to use when running parallel, if this is not set (or
            equal to 0) and parallel set to true, then 0.75 times the maximum number of cores
            available will be used. Otherwise it will try and use the number of cores given
            if this is possible.
        :type num_cores: int
        :param parallel_timeout: Seconds a single structure may take in a parallel worker
            before the run is treated as stalled. The timeout is on the wait for the next
            result rather than on the run as a whole, so a long measurement campaign is not
            interrupted as long as it keeps producing results. A worker that stops making
            progress raises a RuntimeError instead of blocking the run forever, which is what
            a bare ``Pool`` does when a worker deadlocks or is killed. Raise it for features
            that are slow per structure, or set it to None to wait indefinitely. Ignored when
            not running in parallel.
        :type parallel_timeout: int, float, None
        '''

        self.activate_log = False
        if activate_log:
            self.activate_log = activate_log
            self.log_path = os.path.join(outdir, log_path)

            self.logger = logging.getLogger('MeasureLog')
            self.logger.setLevel(level = logging.DEBUG)

            formatter = logging.Formatter('%(message)s') # as simple as possible
            handler = logging.FileHandler(self.log_path, encoding = 'UTF-8')
            handler.setLevel(logging.INFO)
            handler.setFormatter(formatter)

            self.logger.addHandler(handler)

        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)
        self.df_input = df_input
        self.folder = os.path.join(outdir, "curated")
        self.only_relaxed = only_relaxed

        self.include_mod = include_modified

        # document failed pdb files
        self.wrong_pdb_file = []
        self.report_errors = report_errors
        if self.report_errors:
            self.error_filename = self._setup_report_errors_file()
        else:
            self.error_filename = 'no_record'

        self.residue_of_interest = residue_of_interest
        self.features_dict = features_dict.copy()

        if 'all' in self.features_dict:
            self.features_dict = {k: {} for k in ALL_FEATURES}

        self.aa_properties = self._setup_aa_properties(residue_of_interest)
        self._setup_measures(self.features_dict.copy())
        pd.set_option("display.max_columns", None)
        pd.reset_option('display.max_rows')

        # for restarting
        self.current_index = 0
        self.progress_index = 0
        self.pdb_only_files_to_ignore = []

        # for parallel measurements
        self.parallel = parallel
        self.num_cores = num_cores
        self.parallel_timeout = parallel_timeout
        self.n_cores_to_use = 1
        if self.parallel:
            print('>> Measurements running in parallel')
            if isinstance(self.num_cores, int):
                if self.num_cores == 0:
                    self.n_cores_to_use = max(1, int(round(cpu_count() * 0.75)))
                else:
                    if self.num_cores <= os.cpu_count():
                        self.n_cores_to_use = self.num_cores
                    else:
                        print(f'>> Given number of cores for parllel running ({self.num_cores}) is not '
                                f'possible on current setup; defaulting to 0.75 times max possible')
                        self.n_cores_to_use = max(1, int(round(cpu_count() * 0.75)))
            else:
                print(f'>> Input given to number of cores is not an integer ({str(self.num_cores)}); '
                        f'defaulting to 0.75 times max possible.')
                self.n_cores_to_use = max(1, int(round(cpu_count() * 0.75)))
        else:
            self.n_cores_to_use = 1
            print('>> Measurements running in series')

        self.files_to_analyse = []
        self.parallel_items = {}

        # Check that all files in DataFrame appear at least once in folder
        files_af=[os.path.splitext(os.path.basename(c))[0] for c in glob.glob(os.path.join(self.folder, "*pdb"))]
        files_pdb=[os.path.basename(c).split("-")[0] for c in glob.glob(os.path.join(self.folder, "*pdb"))]
        for f in df_input['PDB_Code'].values:
            if f not in files_af and f not in files_pdb:
                print(f'WARNING: {f} not found in folder {self.folder}')

        self.pkaoutdir = os.path.join(outdir, "propkaoutput")
        if not os.path.exists(self.pkaoutdir):
            os.makedirs(self.pkaoutdir)

        self.PDB_only = False

        if 'Uniprot_Entry' in self.df_input.columns:
            columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid', 'Source',
                       'Min_Dist_Other_Chain']
            self.df = pd.DataFrame(columns=columns)
        else:
            self.PDB_only = True
            columns = ['PDB_Code', 'Chain', 'Resid', 'Source', 'Min_Dist_Other_Chain']
            self.df = pd.DataFrame(columns = columns)


    @staticmethod
    def _source_of(path):
        '''
        Record which copy of a structure a measurement was taken from.

        Whether a structure is measured relaxed or unrelaxed is not a property of the
        protein but of how curation went: minimisation is skipped when openmm has no
        template for something the structure retains, and fails outright on some
        structures, and ``only_relaxed`` then falls back to the unminimised file. Without
        this column a measurements table silently mixes the two.

        :param path: file the measurement was taken from.
        :type path: str
        :returns: 'relaxed' if the file is the energy-minimised copy, else 'unrelaxed'.
        :rtype: str
        '''
        stem = os.path.splitext(os.path.basename(path))[0]
        return 'relaxed' if stem.endswith('_relaxed') else 'unrelaxed'


    def _other_chain_distances(self, M):
        '''
        Distance from each residue of interest to the closest other protein chain, and which
        residues can be featurised at all.

        A residue is featurised at its anchor atom, the ``atom_select_names_nonmod`` (or,
        for a modified code, ``atom_select_names_modified``) entry of ``aa_properties``. A
        residue whose anchor is absent cannot be featurised and is missing from the returned
        dictionary, which is how the callers leave it out of the table. When a residue has
        more than one anchor atom the shortest distance is kept.

        :param M: the structure being measured.
        :type M: biobox.Molecule
        :returns: {(chain, resid): distance in A}, NaN when the structure holds no other
            protein chain.
        :rtype: dict
        '''
        d = M.data
        resname = d['resname'].astype(str).str.strip().to_numpy()
        name = d['name'].astype(str).str.strip().to_numpy()
        anchor = (np.isin(resname, self.aa_properties['non_modified_codes'])
                  & np.isin(name, self.aa_properties['atom_select_names_nonmod']))
        if self.include_mod:
            anchor |= (np.isin(resname, self.aa_properties['modified_codes'])
                       & np.isin(name, self.aa_properties['atom_select_names_modified']))

        protein = (PROTEIN_RESNAMES | set(self.aa_properties['non_modified_codes'])
                   | set(self.aa_properties['modified_codes']))
        dist = distance_to_other_chains(d['chain'].to_numpy(), resname,
                                        d['atomtype'].to_numpy(), M.points, anchor, protein)

        out = {}
        idx = np.flatnonzero(anchor)
        for c, r, v in zip(d['chain'].to_numpy()[idx], d['resid'].to_numpy()[idx], dist):
            key = (c, int(r))
            out[key] = np.fmin(out[key], v) if key in out else v
        return out


    def _setup_measures(self, features_dict):
        '''
        Convert a dict of features into a measuring protocol. If ['all'] given as input for the
        features, this will convert the features list to a list containing all current possible
        features. Dictionary is used for this such that the user can provide optional parameters
        directly to the measurements classes without any manual editing.

        :param features_dict: The list of features that are required to measure over the set of
            proteins
        :type features_dict: dict
        '''
        if 'all' in features_dict:
            features_dict = {k: {} for k in ALL_FEATURES}
            self.features_dict = features_dict
            self.features = list(features_dict)
        self.measures = []
        #: {feature: settings it was constructed with}, written to MEASURE_SETTINGS_FILE
        self.feature_settings = {}
        melodia_features = []
        melodia_added = False
        frustration_added = False
        meas_dict = {}
        self.features = list(features_dict)
        for m in features_dict:
            if m in ['frustration', 'density']:
                try:
                    if not frustration_added:
                        frustration = FRUSTRATION(include_modified=self.include_mod,
                                                  error_filename=self.error_filename,
                                                  aa_properties=self.aa_properties)
                        self.measures.append(['frustration', frustration.calculate])
                        self.feature_settings['frustration'] = _resolved_settings(FRUSTRATION, {})
                        frustration_added = True
                except Exception as e:
                    self.features.remove(m)
                    print(f'>> Failed to add frustration/density for features calculation list; error: {e}')
            elif m == 'melodia':
                try:
                    structure = STRUCTURE(melodia_features=['all'],
                                          include_modified=self.include_mod,
                                          error_filename=self.error_filename,
                                          aa_properties=self.aa_properties)
                    self.measures.append([m, structure.calculate])
                    self.feature_settings['melodia'] = _resolved_settings(
                        STRUCTURE, {'melodia_features': ['all']})
                    melodia_added = True
                    self.features += ['curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi']
                    self.features.remove('melodia')
                except Exception as e:
                    self.features.remove(m)
                    print(f'>> Failed to add melodia for features calculation list; error: {e}')
            elif m in ['curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi']:
                melodia_features += [m]
            else:
                if m.upper() in globals().keys():
                    cls_dets = globals()[m.upper()]
                    if not (inspect.isclass(cls_dets) and callable(getattr(cls_dets, 'calculate', None))):
                        self.features.remove(m)
                        print(f'>> Feauture class {m} has no callable attribute called calculate()'
                              f', so {m} cannot be measured, please ensure the feature class has a'
                              f'calculate() function. This feature has been removed from the list '
                              f'to calculate.')
                        continue
                    if inspect.isclass(cls_dets) and hasattr(cls_dets, 'calculate') and callable(getattr(cls_dets, 'calculate')):
                        try:
                            tmp_name = m.lower()
                            kwargs = features_dict.get(m, {})
                            if tmp_name in ['propka', 'legolas'] or 'outdir' in list(kwargs):
                                kwargs['outdir'] = self.outdir
                            if tmp_name == 'rmsf' and 'df_proteins' not in list(kwargs):
                                kwargs['df_proteins'] = self.df_input

                            meas_dict[tmp_name] = globals()[m.upper()](include_modified=self.include_mod,
                                                    error_filename=self.error_filename,
                                                    aa_properties=self.aa_properties,
                                                    **kwargs)

                            if tmp_name == 'evolution':
                                meas_dict[tmp_name].check_esm_model_available()

                            self.measures.append([m, meas_dict[tmp_name].calculate])
                            self.feature_settings[tmp_name] = _resolved_settings(cls_dets, kwargs)
                        except Exception as e:
                            if self.report_errors:
                                self._report_error_to_file('Setup measures: custom measure failed to be added', 'setup', f'Custom measure {m} failed to be added')
                            self.features.remove(m)
                            print(f'Failed to add measure feature {m}, please check new scripts follow the template correctly. Error: {e}')
                            continue
                else:
                    if self.report_errors:
                        self._report_error_to_file('Setup measures: measure unknown', 'setup', f'Measure {m} unknown')
                    self.features.remove(m)
                    print(f'>> Measure {m} unknown and could not load in, removed this from features to calculate'
                          f'Available features: {sorted(a for a in globals() if a.isupper())}')

        if not melodia_added and melodia_features:
            try:
                structure = STRUCTURE(melodia_features=melodia_features,
                                      include_modified=self.include_mod,
                                      error_filename=self.error_filename,
                                      aa_properties=self.aa_properties)
                self.measures.append(['melodia', structure.calculate])
                self.feature_settings['melodia'] = _resolved_settings(
                    STRUCTURE, {'melodia_features': melodia_features})
                melodia_added = True
            except Exception as e:
                for feat in melodia_features:
                    if feat in self.features:
                        self.features.remove(feat)
                print(f'>> Failed to add melodia for features calculation list; error: {e}')


    def _setup_aa_properties(self, res_details):
        '''
        Resolve the residue of interest into the properties every feature works from, and drop
        the features that cannot act on that residue.

        .. rubric:: Method

        - :func:`resolve_residue <resdy.residues.resolve_residue>` accepts a three-letter code
          of any of the twenty standard amino acids, matched case-insensitively against
          :data:`AA_PRESETS <resdy.residues.AA_PRESETS>`, or a dictionary whose four keys may
          be given in any order. Presets and dictionaries are validated by the same rule.
        - Each requested feature is then asked, through the optional class attributes described
          in :mod:`resdy.residues`, whether it can act on that residue. A feature declaring a
          ``SUPPORTED_RESIDUES`` set that excludes the residue is dropped, quoting its own
          ``UNSUPPORTED_REASON``; one declaring ``RESIDUE_KWARGS`` is given the arguments of
          that residue. Nothing here names a particular feature, so a feature added by dropping
          a file into ``features/`` participates without an edit to this method.

        :param res_details: Three-letter code of the residue to measure, or a full properties
            dictionary.
        :type res_details: str, dict
        :returns: The validated properties of the residue of interest.
        :rtype: dict
        '''
        aa_properties = resolve_residue(res_details)
        res_code = residue_key(aa_properties)

        for name in list(self.features_dict):
            cls = self._feature_class(name)
            if cls is None:
                # an unknown feature, or one whose module failed to import; _setup_measures
                # reports it and removes it
                continue

            supported = getattr(cls, 'SUPPORTED_RESIDUES', None)
            if supported is not None and res_code not in supported:
                reason = getattr(cls, 'UNSUPPORTED_REASON', None)
                if isinstance(reason, dict):
                    reason = reason.get(res_code)
                print(f'>> {name} cannot be calculated for {res_code}'
                      + (f', because {reason}' if reason else '')
                      + f'; removing {name} from the features to calculate.')
                self.features_dict.pop(name)
                continue

            defaults = getattr(cls, 'RESIDUE_KWARGS', {}).get(res_code)
            if defaults:
                # take user's own settings over the default values for the residue,
                # rebuild dictionary to ensure not writing to caller held one
                self.features_dict[name] = {**copy.deepcopy(defaults),
                                            **self.features_dict[name]}

        return aa_properties


    def _feature_class(self, name):
        '''
        The class implementing a feature, looked up the same way :meth:`_setup_measures` does.

        :param name: Name of the feature as it appears in features_dict.
        :type name: str
        :returns: The class, or None when no class of that name is available, which covers both
            an unknown feature and one whose module failed to import for want of an optional
            dependency.
        :rtype: type
        '''
        cls = globals().get(name.upper())
        if inspect.isclass(cls) and callable(getattr(cls, 'calculate', None)):
            return cls
        return None


    def _setup_report_errors_file(self):
        '''
        Function to set up the file where errors produced through running the Measure class will be
        written to such that they are easier to look over after running, rather than trawling
        through output.
        '''
        new_file_name = os.path.join(self.outdir, f'measures_errors_{date.today()}.txt')
        while os.path.exists(new_file_name):
            if '_no' in new_file_name:
                error_file_num = int(os.path.splitext(new_file_name)[0].split('_no')[-1])
                new_file_name = os.path.join(self.outdir, f'measures_errors_{date.today()}_no{(error_file_num + 1)}.txt')
            else:
                new_file_name = os.path.join(self.outdir, f'measures_errors_{date.today()}_no{1}.txt')
        with open(new_file_name, 'w') as error_f1:
            error_f1.write(f'New measures errors file created at {datetime.datetime.now()}\n')
        return new_file_name


    def measurement_settings(self):
        """
        Everything that decides what a measurement means: the settings each feature was
        constructed with, the residue measured, and the versions of RESDY and biobox.

        :rtype: dict
        """
        from . import __version__ as resdy_version
        return json.loads(json.dumps({
            'resdy_version': resdy_version,
            'biobox_version': getattr(bb, '__version__', None),
            'residue_of_interest': self.aa_properties,
            'include_modified': self.include_mod,
            'only_relaxed': self.only_relaxed,
            'features': self.feature_settings,
        }, default=str))

    def _write_settings(self):
        """Write :meth:`measurement_settings` to MEASURE_SETTINGS_FILE in the output directory."""
        with open(os.path.join(self.outdir, MEASURE_SETTINGS_FILE), 'w') as fh:
            json.dump(self.measurement_settings(), fh, indent=2)

    def _check_settings(self):
        """
        Refuse to continue a measurement with settings other than those it was started with,
        which would leave rows measured two different ways in one table.

        :raises ValueError: when MEASURE_SETTINGS_FILE exists and differs from the present
            settings.
        """
        path = os.path.join(self.outdir, MEASURE_SETTINGS_FILE)
        if not os.path.exists(path):
            return
        with open(path) as fh:
            previous = json.load(fh)
        current = self.measurement_settings()
        if previous != current:
            changed = sorted(k for k in set(previous) | set(current)
                             if previous.get(k) != current.get(k))
            raise ValueError(f'>> The measurement being restarted was made with different '
                             f'settings ({", ".join(changed)} differ from {path}). Restart it '
                             f'with the settings recorded there, or start a new measurement.')

    def measure_data(self):
        '''
        Determine the appropriate measures function to call based on the combination of running
        PDB_only and in parallel, reducing the number of individual functions that the user will
        have to call themselves.

        The settings of the measurement are written to MEASURE_SETTINGS_FILE in the output
        directory, see :meth:`measurement_settings`.
        '''
        self._write_settings()
        match (self.PDB_only, self.parallel):
            case (False, True) | (False, False):
                # not PDB only and parallel or series:
                self.measure_dataframe()
            case (True, True) | (True, False):
                # PDB only and parallel or series:
                self.measure_PDB_only()

        plddt_record_path = os.path.join(self.folder, 'AF_PLDDT_Output.csv')
        if os.path.exists(plddt_record_path):
            df_af_plddt = pd.read_csv(plddt_record_path).drop_duplicates(
                subset=['PDB_Code', 'Chain', 'Resid'], keep='last')
            self.df['_plddt_key'] = self.df['PDB_Code'].str.replace('_relaxed', '', regex=False)
            df_af_plddt = df_af_plddt.rename(columns={'PDB_Code': '_plddt_key'})
            self.df = self.df.merge(df_af_plddt, how='left', on=['_plddt_key', 'Chain', 'Resid'])
            self.df = self.df.drop(columns='_plddt_key')
            n_af = int(self.df['PDB_Code'].str.contains('AF-').sum())
            if n_af and not int(self.df['PLDDT'].notna().sum()):
                print(f'>> {n_af} AlphaFold rows were measured but none matched a PLDDT record '
                      f'in the PLDDT record path given {plddt_record_path}')
        else:
            print(f'>> No PLDDT record file available at {plddt_record_path}, '
                  f'no PLDDT column added to measurement dataframe')

        self._cleanup_calculation_files()


    def restart_measure_data(self):
        '''
        Determine the appropriate measures function to call based on the combination of running
        PDB_only and in parallel, reducing the number of individual functions that the user will have
        to call themselves. Different to measure_data() as this will restart the measurements from
        final previous point rather than starting again. Raises a ValueError if the settings
        differ from those the measurement was started with.
        '''
        self._check_settings()
        self._write_settings()
        match (self.PDB_only, self.parallel):
            case (False, False) | (False, True):
                # not PDB only and parallel or series
                self.restart_measure()
            case (True, True) | (True, False):
                # PDB only and parallel or series
                self.restart_measure_pdb_only()

        self._cleanup_calculation_files()


    def _report_error_to_file(self, measurement_stage, path, error):
        '''
        Helper function to remove redundant code writing errors in the measurements to the
        measurement error log file.

        :param measurement_stage: The stage of measurements that has caused the error with the file,
            eg propka 1
        :type measurement_stage: str
        :param path: The path of the pdb file that the measurement has been attempted on
        :type path: str
        :param error: The error that has been produced at that step of the measurement when it has
            been attempted to extract features from the pdb file
        :type error: str
        '''
        with open(self.error_filename, 'a', encoding='utf-8') as e_f:
            e_f.writelines('--------------------------------------------------------------------------\n')
            e_f.writelines(f'{measurement_stage} calc error\n')
            e_f.writelines(path + '\n')
            e_f.writelines(error + '\n')


    def save_state(self, outname="measures.csv"):
        '''
        Function saves a csv file of all of the measurements calculated through measure_dataframe()
        File automatically saved in the output directory that has been set previously when setting
        up the measures class Option to customise the name of the output file through outname
        parameter

        :param outname: the name of the csv file that the output is written to
        :type outname: str
        '''
        # sort by uniprot code to give order to output after parallel run
        if not self.PDB_only:
            self.df = self.df.sort_values(by='Uniprot_Entry')
        self.df.to_csv(os.path.join(self.outdir, outname), index_label=False, index=False)


    @staticmethod
    def _parallel_context():
        '''
        Return the multiprocessing context the worker pool is built from.

        A pool is deliberately not built on the 'fork' start method, which is the default on
        Linux. Several features import torch (via the aev, evolution and legolas back ends),
        and torch leaves its thread pools running in the parent for the rest of the session.
        Forking a multi-threaded parent gives the child a copy of those pools with none of
        their threads, so the first call that reaches the inherited OpenMP or BLAS runtime can
        block forever. 'forkserver' is preferred, as its children are forked from a clean
        single-threaded server process, and 'spawn' is the fallback on platforms without it.

        :returns: A multiprocessing context whose start method is not 'fork', where one is
            available.
        '''
        available = mp.get_all_start_methods()
        for method in ('forkserver', 'spawn'):
            if method in available:
                return mp.get_context(method)
        return mp.get_context()


    def _ensure_logger(self):
        '''
        Reattach the log file handler if the current process does not have it.

        A ``logging.Logger`` pickles as a lookup of its name, so a worker started by one of
        the non-forking methods (see :meth:`_parallel_context`) receives the logger without
        the handler that was added to it in the parent, and anything it logs would go nowhere.
        The handler is therefore recreated on first use inside the worker. Records stay
        serialised by the lock that is held around the logging calls, so the processes do not
        interleave their writes to the file.
        '''
        if not self.activate_log:
            return

        self.logger = logging.getLogger('MeasureLog')
        if not self.logger.handlers:
            self.logger.setLevel(level=logging.DEBUG)
            formatter = logging.Formatter('%(message)s')
            handler = logging.FileHandler(self.log_path, encoding='UTF-8')
            handler.setLevel(logging.INFO)
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)


    def _run_in_parallel(self, worker, items):
        '''
        Run ``worker`` over ``items`` in a pool of processes, giving up on a stalled worker.

        Results are collected through ``imap`` rather than ``starmap`` so that they arrive one
        at a time and the wait for each one can be bounded by ``parallel_timeout``. The timeout
        therefore applies to the gap between results, not to the run as a whole. We note that a
        plain ``Pool.starmap`` never returns if one of its workers deadlocks or is killed, so a
        single wedged structure otherwise consumes the whole run.

        Items are handed out one at a time rather than in chunks, for two reasons. ``imap``
        only returns an iterator whose ``next`` accepts a timeout when the chunk size is 1,
        falling back to a plain generator otherwise, and measuring one structure is coarse
        enough work that per-item dispatch costs nothing while balancing the load better.

        :param worker: Bound method to call for each item, taking the unpacked item.
        :type worker: callable
        :param items: One sequence of positional arguments per structure to measure.
        :type items: list
        :returns: The return value of ``worker`` for each item, in the order of ``items``.
        :rtype: list
        :raises RuntimeError: if no further result arrives within ``parallel_timeout`` seconds.
        '''
        ctx = self._parallel_context()
        results = []
        with ctx.Pool(self.n_cores_to_use, maxtasksperchild=20) as pool:
            iterator = pool.imap(partial(_call_with_args, worker),
                                 [tuple(item) for item in items], chunksize=1)
            while True:
                try:
                    results.append(iterator.next(timeout=self.parallel_timeout))
                except StopIteration:
                    break
                except mp.TimeoutError:
                    pool.terminate()
                    raise RuntimeError(
                        f'A parallel measuring worker produced no result for '
                        f'{self.parallel_timeout} s after {len(results)} of {len(items)} '
                        f'structures, so it is treated as stalled and the run was stopped. '
                        f'Either raise parallel_timeout if the features requested are simply '
                        f'slow on these structures, or set parallel=False to measure in '
                        f'series.') from None
        return results


    def measure_dataframe(self):
        '''
        Function to measure specified features for all the structure files curated earlier in the
        programme. Will take a list of the required proteins, finds associated curated structures
        and runs the required measurement functions. Results are saved to memory and a log file
        produced at the same time. (M.save_state() can be used to save the data to a csv). This
        either runs in series or in parallel based on the setting of parallel in the measure class
        initialisation.

        .. rubric:: Method

        Create list of files that have been curated into the self.outdir directory. Iterate over the
        list of the files, check if structure file is
        '''
        if self.PDB_only:
            return 'Call PDB_only method instead'

        files = glob.glob(os.path.join(self.folder, "*pdb"))
        files = [file for file in files if 'pkaani' not in file]
        self.files_to_analyse = files

        total_structures = len(files)
        print(f'Total number of structures to analyse: {total_structures}')

        with Manager() as manager:
            lock = manager.Lock()
            items = []
            for i, r in self.df_input.iterrows():
                pdb_code = r["PDB_Code"]
                chains = r["Chains"]
                if isinstance(chains, str):
                    chains = [c for c in chains.split('/') if c]
                else:
                    chains = []
                method = r['Method']
                res = r['Resolution']
                uniprot_code = r["Uniprot_Entry"]
                file_details = [uniprot_code, pdb_code, method, res, chains]
                items.append([file_details, lock])

            base_cols = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid']
            df_parallel = pd.DataFrame()

            gpu_feats = ['aev', 'evolution', 'legolas']
            gpu_measurements = [a for a in self.measures if a[0] in gpu_feats]
            cpu_measurements = [a for a in self.measures if a[0] not in gpu_feats]

            self.measures = cpu_measurements
            if cpu_measurements:
                if self.parallel:
                    df_parallel = pd.concat(self._run_in_parallel(self._measure_file, items),
                                            ignore_index=True)
                else:
                    df_parallel = pd.concat([self._measure_file(file_dets, lk) for file_dets, lk in items], ignore_index=True)

            self.measures = gpu_measurements
            if gpu_measurements:
                df_gpu = pd.concat([self._measure_file(d, lock) for d, lock in items], ignore_index=True)

                if df_parallel.empty:
                    df_parallel = df_gpu
                else:
                    df_gpu = df_gpu.drop(columns=[a for a in df_gpu.columns if a in df_parallel.columns and a not in base_cols])
                    df_parallel = df_parallel.merge(df_gpu, how='left', on=[a for a in base_cols if a in df_gpu.columns])

            if not cpu_measurements and not gpu_measurements:
                print('>> No measurements are registered, nothing to calculate')

            self.measures = cpu_measurements + gpu_measurements
            # skip empty frames
            frames = [f for f in (self.df, df_parallel) if not f.empty]
            if frames:
                self.df = pd.concat(frames, ignore_index=True).reset_index(drop=True)

        try:
            self.df.drop_duplicates(subset=None, keep='first', inplace=True, ignore_index=True)
            print('\n>> Removed duplicates from measurement dataframe.')
        except Exception as e:
            print(f'\n>> Failed to remove duplicates from measurement dataframe: {e}')
            if self.report_errors:
                self._report_error_to_file('Failed to remove duplicates from measurement dataframe', 'parallel measures', str(e))


    def _measure_file(self, file_details, lock):
        '''
        Take a file and calculate the required measurements for this. Return the dataframe of the
        calculations to the overall self.df

        .. rubric:: Method

        Take the given information about the file and through the given file_list, find the files to
        analyse. Through the structure, identify all lysines and create a temporary dataframe for
        the results. Iterate over the required measurements (self.measures) and insert results into
        the temporary dataframe. Append temporary dataframe to main dataframe.

        :param file_details: list of details for the file that has been selected to be calculated,
            taking the form of [uniprot_code, pdb_code, method, res, chains]
        :type file_details: list
        :param lock: lock used to stop processes writing to output files and dataframes at the same
            time
        :type lock: multiprocessing manager lock
        '''
        uniprot_code, pdb_code, method, res, chains = file_details

        # calculate features values from all PDB files associated with specific DataFrame entry
        frames_df_list = []
        for f in self.files_to_analyse:
            if pdb_code.lower() != os.path.basename(f).split("-")[0].lower():
                if 'AF-' in f and 'AF-' in pdb_code:
                    if pdb_code.split("-")[1].lower() != os.path.splitext(os.path.basename(f))[0].split("-")[1].lower():
                        continue
                else:
                    continue

            if '_cleaned.pdb' in f:
                continue

            if self.only_relaxed:
                if ('_relaxed' not in f and
                    os.path.exists(os.path.join(self.folder, f'{os.path.splitext(os.path.basename(f))[0]}_relaxed.pdb'))):
                    continue

            terminal_out_statements = []
            tstart = time.time()
            terminal_out_statements.append(f"\n> File: {f}")

            columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid', 'Source',
                       'Min_Dist_Other_Chain']
            if self.include_mod:
                columns.append('Modified')
            records = []

            try:
                M = bb.Molecule()
                M.import_pdb(f, include_hetatm=True)
            except Exception as e:
                self.wrong_pdb_file.append(f)
                if self.report_errors:
                    self._report_error_to_file('Failed to produce bb for pdb file', 'measure file parallel', str(e))
                terminal_out_statements.append(f'Failed to produce bb for pdb file with error: {e}')
                continue

            resnames_to_explore = list(self.aa_properties['non_modified_codes'])
            if self.include_mod: resnames_to_explore += list(self.aa_properties['modified_codes'])
            _, idxs = M.atomselect('*', resnames_to_explore, ['CA'], get_index=True, use_resname=True)
            other_chain = self._other_chain_distances(M)
            no_anchor = 0
            for i in idxs:

                #save only lysine entries from chain of interest
                if M.data["chain"].values[i] not in chains:
                    continue

                key = (M.data['chain'].values[i], int(M.data['resid'].values[i]))
                if key not in other_chain:
                    no_anchor += 1
                    continue

                data = ({'Uniprot_Entry': uniprot_code,
                    'PDB_Code':os.path.splitext(os.path.basename(f))[0],
                    'Method': method,
                    'Resolution': res,
                    'Chain': M.data['chain'].values[i],
                    'Resid': M.data['resid'].values[i],
                    'Source': self._source_of(f),
                    'Min_Dist_Other_Chain': other_chain[key]})

                if self.include_mod:
                    data['Modified'] = (M.data['resname'].values[i] in self.aa_properties['modified_codes'])

                records.append(data)

            df_currentfile = pd.DataFrame.from_records(records, columns=columns)
            terminal_out_statements.append(f">> {len(df_currentfile)} lysines of interest found")
            if no_anchor:
                terminal_out_statements.append(f">> {no_anchor} residue(s) of interest left out: "
                                               f"anchor atom missing, cannot be featurised")

            for meas in self.measures:
                terminal_out_statements.append(f">> evaluating {meas[0]}...")
                try:
                    out_print_trap = io.StringIO()
                    with redirect_stdout(out_print_trap):
                        result = meas[1](f)
                    terminal_out_statements.append(out_print_trap.getvalue())
                    df_currentfile = self._combine_dataframes(df_currentfile, result, meas[0]) #insert measures into temporary DataFrame

                except Exception as e:
                    if self.report_errors:
                        self._report_error_to_file('Meas feat error', f'measure file parallel; feature {meas[0]}', str(e))
                    terminal_out_statements.append(f"ERROR: {e}")
                    continue

            processing_time = round((time.time()-tstart), 2)
            terminal_out_statements.append(f">> file processed in {processing_time} seconds.")

            with lock:
                for statement in terminal_out_statements:
                    print(statement)

                if self.activate_log:
                    self._ensure_logger()
                    if df_currentfile.empty is False:
                        pd.set_option('display.max_colwidth', None,
                                      'display.width', None,
                                      'max_seq_items', None,
                                      "display.max_rows", None)
                        try:
                            self.logger.info(df_currentfile)
                            self.logger.info('--------------------------------------------------------------------------')
                        except Exception as e:
                            if self.report_errors:
                                self._report_error_to_file('logging', 'measure file parallel', str(e))
                            print(f'Error in logging: {e}')

                        # reset the pandas display options back to default for regular displaying
                        pd.reset_option('display.max_colwidth')
                        pd.reset_option('display.width')
                        pd.reset_option('max_seq_items')
                        pd.reset_option('display.max_rows')


            if not df_currentfile.empty:
                frames_df_list.append(df_currentfile)

        if not frames_df_list:
            return pd.DataFrame(columns=['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid',
                                         'Source', 'Min_Dist_Other_Chain'])

        return pd.concat(frames_df_list, ignore_index=True)


    def recover_from_log(self, log_path):
        '''
        Take the log file produced through running measure_dataframe() and convert this to a csv

        .. rubric:: Method

        - Read in the log file (measure_log.txt)
        - Work out the columns from the header
        - If the headers can't be worked out, ask for input to match up columns
        - Read in data
        - Sets self.df to be the data output recovered from the log file.

        :param log_path: the file name for the log file to convert
        :type log_path: str
        :returns: Dataframe containing all the measurements that were in the given log file
        :rtype: pandas.DataFrame
        '''
        if self.PDB_only:
            return 'Function not callable.'

        base_columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid']
        log_to_df = pd.DataFrame(columns=base_columns)
        log_path = os.path.join(self.outdir, log_path)
        print(f'Recovering measured data from file: {log_path}')
        print('WARNING: could take up to a few minutes depending on the number of measurements completed.')

        test_lines = 0
        columns_all_set = False
        potential_col_names = {'1': 'propka', '2': 'pkaANI', '3': 'sasa',
                               '4': 'depth', '5': 'aev', '6': 'das'}
        with open(log_path, "rb") as f:
            num_lines = sum(1 for _ in f)
        curr_line = 0

        with open(log_path) as inf:
            for line in inf:
                curr_line += 1
                # check if it is a header line, check if doesn't start with number or -
                if line[0].isalpha() or line[0] == ' ':
                    parts = line.split()
                    # check that columns have been written to the log file correctly
                    if columns_all_set:
                        continue
                    elif len(parts) <= 6 and not columns_all_set:
                        print('Columns were not set correctly in the log file.')
                        print(f'The first 6 columns are assumed to be: {base_columns}')
                        continue
                    elif len(parts) >= 6 and not columns_all_set:
                        # if all seems correct with the writing check that all the columns can be found in the current columns, if not, add in
                        for part in parts:
                            if part not in base_columns:
                                base_columns.append(part)
                        continue

                line_splitter_bool = all(a == '-' for a in line.strip())
                if line_splitter_bool:
                    line = ''

                # split information into parts keeping the AEV as one unit
                parts = re.split(r'([\w.,\/-]+)|(\[.+?\])', line)
                if len(parts) == 0:
                    continue

                parts = [elmnt.strip() for elmnt in parts if elmnt is not None]
                parts = [elmnt for elmnt in parts if elmnt != '']
                parts = parts[1:]

                if len(parts) == 0:
                    continue
                num_parts = len(parts)
                if num_parts != len(base_columns):
                    while num_parts != len(base_columns):
                        print('Need to set a column header')
                        print(f'Options for columns are: {potential_col_names}')
                        print(f'Please enter the number corresponding to the header required for the column which contains the following value: {parts[len(base_columns)]}')
                        new_header_val = input('Enter the number for the new column header: ')
                        while True:
                            if not new_header_val.isnumeric():
                                new_header_val = input('Enter the number for the new column header: ')
                            elif 1 <= int(new_header_val) <= len(potential_col_names):
                                break
                            else:
                                new_header_val = input('Enter the number for the new column header: ')
                        base_columns.append(potential_col_names[new_header_val])
                    columns_all_set = True
                data = dict(zip(base_columns, parts))
                log_to_df = pd.concat([log_to_df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                test_lines += 1
                last_reported = 0
                percent_prog = round((curr_line/num_lines)*100, 2)
                if percent_prog != last_reported:
                    last_reported = percent_prog
                    print(f'Progress analysing log file: {percent_prog} %\r', end='', flush=True)

        self.df = log_to_df
        print('Data recovered from log file')
        print(f'Numer of measurements read: {len(log_to_df)}')
        return log_to_df


    def restart_measure(self, log_path='measure_log.txt'):
        '''
        A function to restart the measurements calculations Useful if the initial run of the
        measurements crashes or gets stuck Works out how far along the simulation was by running an
        analysis of the measures log file

        .. rubric:: Method

        - Read over the measures log file and collate a list of files that have been analysed
        - Remove the final value from the list as this may not have been done properly
        - Remove completed files from files to do
        - Restart measure_dataframe() with the new list

        :param log_path: The name of the measures log file By default takes the name
            'measures_log.txt'
        :type log_path: str
        '''

        if self.PDB_only:
            return 'restart_measure() function not callable when using PDB_only'

        print('Restarting measurements')
        # 1. Analyse the measures log file to create a list of files that were analysed
        df_prev = self.recover_from_log(log_path=log_path)
        print('Recovered measurement data from log file to workout which protein are left...')

        # remove the last protein from list incase it wasn't completed fully
        measured_proteins = list(self.df['Uniprot_Entry'])
        if measured_proteins:
            final_protein = measured_proteins[-1]
            proteins_completed = list(set(measured_proteins))
            proteins_completed = [c for c in proteins_completed if c != final_protein]

            # 2. Update df_input to only have the files which haven't been analysed yet
            idx_to_remove = []
            for i, r in self.df_input.iterrows():
                if r['Uniprot_Entry'] in proteins_completed:
                    idx_to_remove.append(i)
            self.df_input = self.df_input.drop(idx_to_remove)

        else:
            print(f'>> No measurements of proteins were identified in the given log file: '
                  f'{log_path}, starting again from scratch')

        # 3. Restart the measure_dataframe() with the new file list
        print(f'Continuing measurements. {len(self.df_input)} proteins to measure.')
        self.measure_dataframe()


    def _combine_dataframes(self, target, to_merge, col_name):
        '''
        Function to combine the dataframe produced by a measurement function into the main dataframe
        containing all the measurements. target is a DataFrame to be filled with data, to_merge
        contains the data. Values to insert are indexed in both array by two columns: Chain and
        Resid.

        :param target: DataFrame to be filled with data.
        :type target: pandas.DataFrame
        :param to_merge: to_merge contains the new data to merge.
        :type to_merge: pandas.DataFrame
        :param col_name: Name of the column which the new data is from.
        :type col_name: str
        '''
        to_merge = to_merge.reset_index(drop=True)

        if self.include_mod and 'Modified' not in to_merge.columns:
            raise KeyError(f'>> Running include modified but Modified column not present for feature: {col_name}')

        for i, r in target.iterrows():

            chain_value = r["Chain"]
            resid_value = r["Resid"]
            if self.include_mod:
                modified_value = r['Modified']

            if self.include_mod:
                idx = np.where((to_merge["Chain"] == chain_value) &
                               (to_merge["Resid"].astype(int) == resid_value) &
                               (to_merge["Modified"].astype(bool) == modified_value))
            else:
                idx = np.where((to_merge["Chain"] == chain_value) &
                               (to_merge["Resid"].astype(int) == resid_value))

            if len(idx[0]) == 0:
                continue
            
            if len(idx[0]) > 1:
                out_print = (f'{col_name}: {len(idx[0])} rows match Chain {chain_value}, Resid '
                             f'{resid_value}; an insertion code has been not managed properly, '
                             f'no guess employed; continuing')
                if self.report_errors:
                    self._report_error_to_file('ambiguous resid insertion key in combining measurse dataframes', col_name, out_print)
                print('>> ' + out_print)
                continue


            # account for measurements that have special cases
            if col_name == 'melodia':
                melodia_features = ['curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi']
                for feature in melodia_features:
                    if feature in self.features and feature in to_merge.columns:
                        target.at[i, feature] = to_merge.loc[idx[0][0], feature]
            elif col_name == 'frustration':
                for feature in ('frustration', 'density'):
                    if feature in self.features and feature in to_merge.columns:
                        target.at[i, feature] = to_merge.loc[idx[0][0], feature]
            else:
                target.at[i, col_name] = to_merge.loc[idx[0][0], col_name]
        return target


    def measure_PDB_only(self):
        '''
        Function to measure specified features for a set of pdb files. Takes a list of pdb files,
        finds associated curated structures and runs the required measurement functions. Results are
        saved to memory and a log file produced at the same time if required. Timing is kept to
        update the predicted time remaining as it goes along. M.save_state() can be used to save
        the data to a csv.
        '''
        if not self.PDB_only:
            print('Called measure_PDB_only() when running not on PDB_only. Call '
                  'measure_dataframe() instead or change to run PDB_only.')
            return 'You have called the wrong method for measuring data, call the general measures function instead'

        files = glob.glob(os.path.join(self.folder, "*pdb"))
        files = [file for file in files if 'pkaani' not in file]
        self.files_to_analyse = files

        self.df_input = self.df_input.drop_duplicates()
        if 'completed' not in self.df_input.columns:
            self.df_input['completed'] = False

        with Manager() as manager:
            lock = manager.Lock()
            items = []
            for pdb_idx, row in self.df_input.iterrows():

                if row['completed']:
                    continue

                pdb_code = row['PDB_Code']
                items.append([pdb_code, lock])

            base_cols = ['PDB_Code', 'Chain', 'Resid']
            df_parallel = pd.DataFrame()

            gpu_feats = ['aev', 'evolution', 'legolas']
            gpu_measurements = [a for a in self.measures if a[0] in gpu_feats]
            cpu_measurements = [a for a in self.measures if a[0] not in gpu_feats]

            self.measures = cpu_measurements
            if cpu_measurements:
                if self.parallel:
                    df_parallel = pd.concat(self._run_in_parallel(self._measure_file_pdbonly, items),
                                            ignore_index=True)
                else:
                    df_parallel = pd.concat([self._measure_file_pdbonly(file_dets, lk) for file_dets, lk in items], ignore_index=True)

            self.measures = gpu_measurements
            if gpu_measurements:
                df_gpu = pd.concat([self._measure_file_pdbonly(d, lock) for d, lock in items], ignore_index=True)

                if df_parallel.empty:
                    df_parallel = df_gpu
                else:
                    df_gpu = df_gpu.drop(columns=[a for a in df_gpu.columns if a in df_parallel.columns and a not in base_cols])
                    df_parallel = df_parallel.merge(df_gpu, how='left', on=[a for a in base_cols if a in df_gpu.columns])

            if not cpu_measurements and not gpu_measurements:
                print('>> No measurements are registered, nothing to calculate')

            self.measures = cpu_measurements + gpu_measurements
            # skip empty frames
            frames = [f for f in (self.df, df_parallel) if not f.empty]
            if frames:
                self.df = pd.concat(frames, ignore_index=True).reset_index(drop=True)

        try:
            self.df.drop_duplicates(subset=None, keep='first', inplace=True, ignore_index=True)
            print('\n>> Removed duplicates from measurement dataframe.')
        except Exception as e:
            print(f'\n>> Failed to remove duplicates from measurement dataframe: {e}')
            if self.report_errors:
                self._report_error_to_file('Failed to remove duplicates from measurement dataframe', 'parallel measures', str(e))



    def _measure_file_pdbonly(self, pdb_code, lock):
        '''
        Take a file and calculate the required measurements for this. Return the dataframe of the
        calculations to the overall self.df. This is the worked for measure_PDB_only.

        :param file_details: list of details for the file that has been selected to be calculated,
            taking the form of [pdb_code]
        :type file_details: list
        :param lock: lock used to stop processes writing to output files and dataframes at the same
            time
        :type lock: multiprocessing manager lock
        '''
        file_df_list = []
        for f in self.files_to_analyse:
            if (pdb_code.lower() != os.path.basename(f).split("-")[0].lower()) and (pdb_code.lower() != os.path.splitext(os.path.basename(f))[0].lower()):
                continue

            if self.only_relaxed:
                f_stem = os.path.splitext(os.path.basename(f))[0]
                if '_relaxed' not in f and os.path.exists(os.path.join(self.folder, f'{f_stem}_relaxed.pdb')):
                    continue

            if f in self.pdb_only_files_to_ignore:
                continue

            terminal_out_statements = []
            tstart = time.time()
            terminal_out_statements.append(f"\n> File: {f}")

            columns = ['PDB_Code', 'Chain', 'Resid', 'Source', 'Min_Dist_Other_Chain']
            if self.include_mod:
                columns.append('Modified')
            records = []

            try:
                M = bb.Molecule()
                M.import_pdb(f, include_hetatm=True)
            except Exception as e:
                terminal_out_statements.append(f'Failed to create biobox molecule for file {f} with error: {e}')
                if self.report_errors:
                    self._report_error_to_file(f'Failed to create bb molecule for file: {f}', 'measure pdb only', str(e))
                self.wrong_pdb_file.append(f)
                continue

            if self.include_mod:
                _, idxs = M.atomselect("*", (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']),
                                        ['CA'], get_index=True, use_resname=True)
            else:
                _, idxs = M.atomselect("*",  self.aa_properties['non_modified_codes'],
                                        ['CA'], get_index=True, use_resname=True)

            other_chain = self._other_chain_distances(M)
            no_anchor = 0
            for i in idxs:
                key = (M.data['chain'].values[i], int(M.data['resid'].values[i]))
                if key not in other_chain:
                    no_anchor += 1
                    continue

                data = ({'PDB_Code': os.path.splitext(os.path.basename(f))[0],
                    'Chain': M.data['chain'].values[i],
                    'Resid': M.data['resid'].values[i],
                    'Source': self._source_of(f),
                    'Min_Dist_Other_Chain': other_chain[key]})

                if self.include_mod:
                    data['Modified'] = (M.data['resname'].values[i]
                                        in self.aa_properties['modified_codes'])

                records.append(data)

            df_currentfile = pd.DataFrame.from_records(records, columns=columns)
            terminal_out_statements.append(f">> {len(df_currentfile)} lysines of interest found")
            if no_anchor:
                terminal_out_statements.append(f">> {no_anchor} residue(s) of interest left out: "
                                               f"anchor atom missing, cannot be featurised")

            for meas in self.measures:
                terminal_out_statements.append(f">> evaluating {meas[0]}...")
                try:
                    out_print_trap = io.StringIO()
                    with redirect_stdout(out_print_trap):
                        result = meas[1](f)
                    terminal_out_statements.append(out_print_trap.getvalue())
                    df_currentfile = self._combine_dataframes(df_currentfile, result, meas[0]) #insert measures into temporary DataFrame

                except Exception as e:
                    if self.report_errors:
                        self._report_error_to_file('Meas feat error', f'measure file parallel; feature {meas[0]}', str(e))
                    terminal_out_statements.append(f"ERROR: {e}")
                    continue

            processing_time = round((time.time()-tstart), 2)
            terminal_out_statements.append(f">> file processed in {processing_time} sec.")

            with lock:
                for statement in terminal_out_statements:
                    print(statement)

                if self.activate_log:
                    self._ensure_logger()
                    if df_currentfile.empty is False:

                        pd.set_option('display.max_colwidth', None,
                                    'display.width', None,
                                    'max_seq_items', None,
                                    "display.max_rows", None)
                        try:
                            self.logger.info(df_currentfile)
                            self.logger.info('--------------------------------------------------------------------------')
                        except Exception as e:
                            if self.report_errors:
                                self._report_error_to_file('logging error', 'measure pdb only', str(e))
                            print(f'Error in logging measurement: {e}')

                        pd.reset_option('display.max_colwidth')
                        pd.reset_option('display.width')
                        pd.reset_option('max_seq_items')
                        pd.reset_option('display.max_rows')

            if not df_currentfile.empty:
                file_df_list.append(df_currentfile)

        if not file_df_list:
            return pd.DataFrame(columns=['PDB_Code', 'Chain', 'Resid', 'Source', 'Min_Dist_Other_Chain'])

        return pd.concat(file_df_list, ignore_index=True)


    def restart_measure_pdb_only(self, log_path='measure_log.txt'):
        '''
        A function to restart the measurements calculations for the PDB only function. Useful if the
        initial run of the measurements crashes or gets stuck. Works out how far along the
        simulation was by running an analysis of the measures log file.

        .. rubric:: Method

        - Read over the measures log file and collate a list of files that have been analysed.
        - Remove the final value from the list as this may not have been done properly.
        - Remove completed files from files to do.
        - Restart measure_dataframe() with the new list.

        :param log_path: The name of the measures log file By default takes the name
            'measures_log.txt'
        :type log_path: str
        '''

        if not self.PDB_only:
            print('>> restart_measure_pdb_only() function not callable when not running PDB only')
            return 'restart_measure_pdb_only() function not callable when not running PDB only'

        print('>> Preparing to restart measurements on PDB only')
        # 1. Analyse the measures log file to create a list of files that were analysed
        tmp_log_path = os.path.join(self.outdir, f'{os.path.splitext(log_path)[0]}_tmp.txt')
        log_path = os.path.join(self.outdir, log_path)
        print(f'>> Finding measured proteins from log file: {log_path}')
        proteins_completed = []
        with open(log_path, "rb") as f:
            num_lines = sum(1 for _ in f)

        if num_lines == 0:
            print('>> No previous measures data is found in the specified log file. Make sure the log file stated is correct or run measure_pdb_only() from start.')
            return

        curr_line = 0
        final_prot_start_line = None
        last_reported = 0
        with open(file=log_path, mode='r') as lpf:
            for line in lpf:
                curr_line += 1
                if line[0].isalpha() or line[0] in [' ', '-']:
                    continue
                parts = line.split()
                protein_code = parts[1].split('/')[-1]
                proteins_completed.append(protein_code)
                final_prot_start_line = curr_line
                percent_prog = round((curr_line/num_lines)*100, 2)
                if percent_prog != last_reported:
                    last_reported = percent_prog
                    print(f'Progress analysing log file: {percent_prog} %\r', end='', flush=True)
        print('>> Measured proteins recovered from log file')

        if final_prot_start_line is None:
            print(f'>> No complete measurement sets found in the log')
        else:
            curr_line = 0
            with open(file=log_path, mode='r') as lpf, open(file=tmp_log_path, mode='w') as npf:
                for line in lpf:
                    curr_line += 1
                    if curr_line < final_prot_start_line:
                        npf.write(line)
            os.replace(src=tmp_log_path, dst=log_path)


        # 2. remove the last protein from list incase it wasn't completed fully
        final_protein = proteins_completed[-1]
        proteins_completed = [c for c in proteins_completed if c != final_protein]
        proteins_completed = list(set(proteins_completed))
        self.progress_index = len(proteins_completed)

        # 3. Update df_input to only have the files which haven't been analysed yet
        if 'completed' not in self.df_input.columns:
            self.df_input['completed'] = False
        old_len_df_input = len(self.df_input)
        idx_to_remove = []
        files_pdbs = [os.path.basename(a) for a in glob.glob(os.path.join(self.folder, "*pdb"))]
        for i, r in self.df_input.iterrows():
            matched_pdb_files = [a.replace('.pdb', '') for a in files_pdbs if r['PDB_Code'] in a]
            for recover_file in proteins_completed:
                # case 1: exact match code and file - for measuring data from simulations mainly
                if r['PDB_Code'] == recover_file:
                    idx_to_remove.append(i)
                    self.df_input.at[i, 'completed'] = True
                    break
                # case 2: PDB files renamed by curation that are not the PDB code alone
                if recover_file in matched_pdb_files:
                    if recover_file not in self.pdb_only_files_to_ignore:
                        full_recover_file = self.folder + os.sep + recover_file + '.pdb'
                        self.pdb_only_files_to_ignore.append(full_recover_file)

            perc_prog_remove = round((i/len(self.df_input))*100, 2)
            print(f'Progress removing measured files: {perc_prog_remove} %\r', end='', flush=True)

        self.df_input = self.df_input.drop(idx_to_remove)
        new_len_df_input = len(self.df_input)
        lines_df_input_removed = old_len_df_input - new_len_df_input
        files_left_to_calc = len(self.df_input) - len(self.pdb_only_files_to_ignore)
        print(f'>> {lines_df_input_removed} exact matches in PDB codes removed from the input '
              f'list that have already been calculated.')
        print(f'>> {len(self.pdb_only_files_to_ignore)} files to ignore in measurements that '
              f'have already been calculated.')
        print(f'>> Continuing measurements. {files_left_to_calc} proteins to measure.')
        self.measure_PDB_only()


    def recover_from_log_PDB_only(self, log_path):
        '''
        Take the log file produced through running measure_PDB_only() and convert this to a csv

        .. rubric:: Method

        - Read in the log file (measure_log.txt) or other given name.
        - Work out the columns from the header.
        - If the headers can't be worked out, ask for input to match up columns.
        - Read in data.
        - Sets self.df to be the data output recovered from the log file.

        :param log_path: the file name for the log file to convert
        :type log_path: str
        :returns: Dataframe containing all the measurements that were in the given log file
        :rtype: pandas.DataFrame
        '''
        if not self.PDB_only:
            return 'Function not callable.'

        base_columns = ['PDB_Code', 'Chain', 'Resid']
        log_to_df = pd.DataFrame(columns=base_columns)
        log_path = os.path.join(self.outdir, log_path)
        print(f'Recovering measured data from file: {log_path}')
        print('WARNING: could take up to a few minutes depending on the number of measurements completed.')

        test_lines = 0
        columns_all_set = False
        potential_col_names = {'1': 'propka', '2': 'pkaANI', '3': 'sasa',
                               '4': 'depth', '5': 'aev', '6': 'das',
                               '7': 'Other'}
        with open(log_path, "rb") as f:
            num_lines = sum(1 for _ in f)
        curr_line = 0

        with open(log_path) as inf:
            set_header_line = ''
            dataframe_columns = []
            for line in inf:
                curr_line += 1
                # check if it is a header line, check if doesn't start with number or -
                if line[0].isalpha() or line[0] == ' ':
                    if line != set_header_line:
                        parts = line.split()
                        # check that columns have been written to the log file correctly
                        if len(parts) <= 3 and not columns_all_set:
                            print('Columns were not set correctly in the log file.')
                            print(f'The first 6 columns are assumed to be: {base_columns}')
                            continue
                        elif len(parts) >= 3 and not columns_all_set:
                            # if all seems correct with the writing check that all the columns can be found in the current columns, if not, add in
                            for part in parts:
                                if part not in base_columns:
                                    base_columns.append(part)
                                if part not in dataframe_columns:
                                    dataframe_columns.append(part)
                                    if len(log_to_df) != 0:
                                        log_to_df[part] = None
                            continue

                line_splitter_bool = all(a == '-' for a in line.strip())
                if line_splitter_bool:
                    line = ''

                # split information into parts keeping the AEV as one unit
                parts = re.split(r'([\w.,\/-]+)|(\[.+?\])', line)
                if len(parts) == 0:
                    continue
                parts = [elmnt.strip() for elmnt in parts if elmnt is not None]
                parts = [elmnt for elmnt in parts if elmnt != '']
                parts = parts[1:]

                if len(parts) == 0:
                    continue
                num_parts = len(parts)

                # Case 1: Setting the columns when the columns have been messed up and aren't the same as the data in the log file
                if num_parts > len(base_columns):
                    while num_parts != len(base_columns):
                        print('Need to set a column header')
                        if len(potential_col_names) != 0:
                            print(f'Options for columns are: {potential_col_names}')
                            print(f'Please enter the number corresponding to the header required for the column which contains the following value: {parts[len(base_columns)]}')
                            new_header_val = input('Enter the number for the new column header: ')
                            while True:
                                if not new_header_val.isnumeric():
                                    new_header_val = input('Enter the number for the new column header: ')
                                elif 1 <= int(new_header_val) <= len(potential_col_names):
                                    break
                                else:
                                    new_header_val = input('Enter the number for the new column header: ')
                            if new_header_val == '7':
                                new_header_name = input('Other selected, please enter a unique name for the column: ')
                                base_columns.append(new_header_name)
                            else:
                                base_columns.append(potential_col_names[new_header_val])
                                del potential_col_names[new_header_val]
                        else:
                            new_header = input(f'No more suggested columns available, please enter your column name for the column containing this value:  {parts[len(base_columns)]}')
                            base_columns.append(new_header)
                    columns_all_set = True
                data = dict(zip(base_columns, parts))
                for col in dataframe_columns:
                    if col not in data:
                        data[col] = None

                log_to_df = pd.concat([log_to_df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                test_lines += 1
                last_reported = 0
                percent_prog = round((curr_line/num_lines)*100, 2)
                if percent_prog != last_reported:
                    last_reported = percent_prog
                    print(f'Progress analysing log file: {percent_prog} %\r', end='', flush=True)

        self.df = log_to_df
        print('Data recovered from log file')
        print(f'Numer of measurements read: {len(log_to_df)}')
        return log_to_df


    def _cleanup_calculation_files(self):
        '''
        Function to remove any temporary or result files created through the calculation of the
        measurements within this class. While all are meant to have been moved at the time of
        calculation, occasionally this fails and leaves some behind. This applies particularly
        with Modeller, Legolas and propka.
        '''
        def _mv_files(files, dest):
            '''
            Generic function for moving list of files over to the destination
            '''
            if not files:
                return
            dest_path = os.path.join(self.outdir, dest)
            os.makedirs(dest_path, exist_ok=True)
            for f in files:
                try:
                    os.rename(f, os.path.join(dest_path, f))
                except Exception as e:
                    print(f'> Failed to move {f} to destination {dest_path} with error: {e}')

        print('\n>> Cleaning up leftover files from measures calculations...')
        dir_files = [f for f in os.listdir() if os.path.isfile(os.path.join(os.getcwd(),f))]
        nmr_cs_file, nmr_parquet_file, propka_pka_file, propka_error_file = [], [], [], []
        for f in dir_files:
            if f.endswith('_cs.csv'):
                nmr_cs_file.append(f)
            elif f.endswith('_cs.parquet'):
                nmr_parquet_file.append(f)
            elif f.endswith('.pka'):
                propka_pka_file.append(f)
            elif f.endswith('_propka_errors.txt'):
                propka_error_file.append(f)
        if nmr_cs_file:
            _mv_files(nmr_cs_file, 'legolas')
        if nmr_parquet_file:
            _mv_files(nmr_parquet_file, 'legolas')
        if propka_pka_file:
            _mv_files(propka_pka_file, 'propkaoutput')
        if propka_error_file:
            _mv_files(propka_error_file, 'propkaoutput')
        if self.report_errors:
            try:
                num_lines = 10
                with open(self.error_filename, 'r') as f:
                    num_lines = sum(1 for _ in f)
                if num_lines <= 2:
                    os.remove(self.error_filename)
            except Exception as e:
                print(f'>> Failed to cleanup the error file ({self.error_filename}) '
                      f'for the measurements run with error: {e}')
        print('>> Unused file cleanup complete.')


if __name__ == "__main__":

    file_one = "Demo{os.sep}curated{os.sep}1M2E-alt-1.pdb"

    from .uniprot import Uniprot
    from .protein import PDB

    UP = Uniprot()
    UP.get_protein_data("P0CG47")

    df = UP.df.iloc[7:9]

    PDB = PDB()
    PDB.gather_proteins(df)

    M = Measure(PDB.df)
    M.measure_dataframe()
