import shutil
import sys
import os
import contextlib
import multiprocessing as mp
import threading


#: Oldest biobox RESDY works with. Earlier versions lack ``Xlink.get_half_sphere``, which the
#: DAS feature calls. Versions before 1.1.5 also wrote the occupancy and B-factor columns of a
#: pdb file into each other's places, and computed solvent accessible surface areas that were
#: too small.
MIN_BIOBOX_VERSION = (2, 0, 0)

#: Columns of the measurements table that are not features: the row key, the structure it
#: came from, and the metadata recorded beside the features.
METADATA_COLUMNS = ('Uniprot_Entry', 'PDB_Code', 'Chain', 'Resid', 'Method', 'Resolution',
                    'Source', 'Modified', 'class', 'Class', 'PLDDT', 'Largest_Gap',
                    'Min_Dist_Other_Chain')


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

class ShutUp(object):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, *args):
        sys.stdout.close()
        sys.stdout =  self._stdout
