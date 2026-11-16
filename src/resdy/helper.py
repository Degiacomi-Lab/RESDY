import shutil
import sys
import os
import contextlib
import multiprocessing as mp
import threading


# Oldest biobox RESDY works with. Earlier versions lack ``Xlink.get_half_sphere``, which the
# DAS feature calls. Versions before 1.1.5 also wrote the occupancy and B-factor columns of a
# pdb file into each other's places, and computed solvent accessible surface areas that were
# too small.
MIN_BIOBOX_VERSION = (2, 0, 0)

#: Columns of the measurements table that are not features
METADATA_COLUMNS = ('Uniprot_Entry', 'PDB_Code', 'Chain', 'Resid', 'Method', 'Resolution',
                    'Source', 'Modified', 'class', 'Class', 'PLDDT', 'Largest_Gap',
                    'Min_Dist_Other_Chain', 'num_structures', 'numbering', 'Dist_To_Modelled',
                    'Dist_To_Nearest_Terminus', 'In_Modelled_Region')

# Prefixes of the provenance columns that Analysis.add_metadata_to_measures joins onto the
# measurements table, one per curation log. Any column which carries one of these prefixes
# is not considered a feature in analysis.
METADATA_PREFIXES = ('geometry_', 'minimisation_', 'assembly_')

RESIDUE_KEY = ['Uniprot_Entry', 'Resid']
RESIDUE_KEY_CHAIN = ['Uniprot_Entry', 'Chain', 'Resid']

FEATURE_LABELS = {'depth': 'Depth (Å)',
                'sasa': 'Solvent Accessible Surface Area (Å\u00b2)',
                'propka': 'pKa (PROPKA3)',
                'pkaani': 'pKa (pKaANI)',
                'flexibility': 'Flexibility (B-factor) (Å\u00b2)',
                'curvature': 'Curvature',
                'arc_length': 'Arc Length',
                'das': 'Dynamically Accessible Surface Area',
                'legolas': '15N NMR Backbone Shift (ppm) (LEGOLAS)',
                'phi': 'Phi',
                'psi': 'Psi',
                'seqcharge': 'Sequence Charge',
                'torsion': 'Torsion Angle',
                'writhing': 'Writhing',
                'frustration': 'Frustration',
                'rmsf': 'Root Mean Square Fluctuation',
                'density': 'Density',
                'aev': 'Atomic Environment Vector',
                'secondarystructure': 'Secondary Structure'}

AGG_LABELS = {'avg': 'Average',
            'med': 'Median',
            'min': 'Minimum',
            'max': 'Maximum',
            'rand': 'Random',
            'sd': 'Standard Deviation',
            'range': 'Range'}

STANDARD_RESIDS = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
                    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
                    'TYR': 'Y', 'VAL': 'V'}


@contextlib.contextmanager
def get_worker_lock(parallel, ctx=None):
    '''
    Return a lock suitable for workers that are going to run when parallel running is required.
    This stops creation of a lock when running in series, reducing the chance of hanging scripts.

    :param parallel: True when work is given to a multiprocessing pool with either protein curation
        or sidechain measurements desired to be ran in parallel
    :type parallel: bool
    :param ctx: Multiprocessing context that the pool is built off such that the Manager created in
        the function sorting the parallelisation and this share the same start method. None will
        use the default multiprocessing context
    '''
    mp.freeze_support()
    if parallel:
        manager_cls = ctx.Manager if ctx is not None else mp.Manager
        with manager_cls() as manager:
            yield manager.Lock()
    else:
        yield threading.Lock()


@contextlib.contextmanager
def mp_hide_main_from_children():
    '''
    In order for the get the multiprocessing context methods to work across all systems
    and reduce the potential for freezing children, need to prevent all children from
    re-importing the user's __main__. By doing this, gets around for placing the onus
    of using freeze_support on the user. This does have a few small caveats here but
    nothing that will affect the actual results coming out.
    '''
    _hide_lock = threading.RLock()
    _saved_main = {}
    _hide_depth = 0
    main = sys.modules["__main__"]

    with _hide_lock:
        if _hide_depth == 0:
            _saved_main["__file__"] = getattr(main, "__file__", None)
            _saved_main["__spec__"] = getattr(main, "__spec__", None)
            main.__file__ = None
            main.__spec__ = None
        _hide_depth += 1

    try:
        yield

    finally:
        with _hide_lock:
            _hide_depth -= 1
            if _hide_depth == 0:
                if _saved_main["__file__"] is not None:
                    main.__file__ = _saved_main["__file__"]
                main.__spec__ = _saved_main["__spec__"]


def non_feature_columns(columns):
    '''
    The names that are not features in a table with the given columns: every name of
    :data:`METADATA_COLUMNS`, and the columns of the table that carry one of the
    :data:`METADATA_PREFIXES`.

    :param columns: column names of a measurements table.
    :type columns: iterable
    :returns: names that must not be treated as features.
    :rtype: list
    '''
    return list(METADATA_COLUMNS) + [c for c in columns if c not in METADATA_COLUMNS
                                     and str(c).startswith(METADATA_PREFIXES)]


def require_biobox(biobox):
    '''
    Raise if the biobox given is older than :data:`MIN_BIOBOX_VERSION`. This ensures
    that use of the RESDY package will get correct output and won't see known bugs from
    previous biobox versions pass through.

    :param biobox: the imported biobox module.
    :type biobox: module
    :raises ImportError: when the version is older than required.
    '''
    version = getattr(biobox, '__version__', None)
    if not isinstance(version, str):
        return
    parts = tuple(int(''.join(c for c in p if c.isdigit()) or 0) for p in version.split('.')[:3])
    # '2.0' counts as 2.0.0
    parts = parts + (0,) * (3 - len(parts))
    if parts < MIN_BIOBOX_VERSION:
        required = '.'.join(map(str, MIN_BIOBOX_VERSION))
        raise ImportError(f'RESDY needs biobox {required} or later, found {version}. Install it '
                          f'with: conda install -c conda-forge "biobox>={required}"')


def get_download_tool():
    '''
    Get the name of the tool used for commandline download

    :returns: Name of the tool to use. Only wget and curl are supported; None is returned if
        neither is found.
    :rtype: str
    '''

    if shutil.which('wget') is not None:
        return 'wget'

    elif shutil.which('curl') is not None:
        return 'curl'

    else:
        return None


def sort_residue_key(df, include_chain=False):
    '''
    Sort out the key for aggregating the residues. Standard procedure for aggregating the residue
    data is to only aggregate on UNIPROT code and resid id, however there may be times aggregating
    by chain is well is necessary. This will return the key to aggregate on. 
    '''
    key = RESIDUE_KEY_CHAIN if include_chain else RESIDUE_KEY
    if 'Uniprot_Entry' not in df.columns:
        key = ['PDB_Code'] + [k for k in key if k != 'Uniprot_Entry']
    return [k for k in key if k in df.columns]


@staticmethod
def parallel_context():
    '''
    Return the multiprocessing context the worker pool is built from. This function is used
    to get the available context and if possible prefer 'forkserver' method first and then
    'spawn' if this is not available. If neither preferred option is available, it will fall
    back onto whatever multiprocessing context is available. 

    With several features importing torch, and torch leaving its thread pools running in the
    parent for the rest of the session. Forking a multi-threaded parent gives the child a copy
    of those pools with none of their threads, so the first call that reaches the inherited
    OpenMP or BLAS runtime can block forever.

    :returns: A multiprocessing context whose start method is not 'fork', where one is
        available.
    '''
    available = mp.get_all_start_methods()
    for method in ('forkserver', 'spawn'):
        if method in available:
            return mp.get_context(method)
    print(">> Could not get 'forkserver' or 'spawn' multiprocessing context," \
        "defaulting to 'fork'")
    return mp.get_context()



class ShutUp(object):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, *args):
        sys.stdout.close()
        sys.stdout =  self._stdout


if __name__ == '__main__':
    mp.freeze_support()
