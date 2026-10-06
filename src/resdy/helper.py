import shutil
import sys, os

#: Oldest biobox RESDY works with. Earlier versions lack ``Xlink.get_half_sphere``, which the
#: DAS feature calls. Versions before 1.1.5 also wrote the occupancy and B-factor columns of a
#: pdb file into each other's places, and computed solvent accessible surface areas that were
#: too small.
MIN_BIOBOX_VERSION = (2, 0, 0)


def require_biobox(biobox):
    '''
    Raise if the biobox given is older than :data:`MIN_BIOBOX_VERSION`.

    :param biobox: the imported biobox module.
    :type biobox: module
    :raises ImportError: when the version is older than required.
    '''
    version = getattr(biobox, '__version__', None)
    if not isinstance(version, str):
        # a mocked module, as in the documentation build
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
