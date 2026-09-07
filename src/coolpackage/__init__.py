'''
Rapid featurisation of amino acids from collections of protein structures.

The pipeline is constituted by five steps, one class each:
:class:`Uniprot <coolpackage.uniprot.Uniprot>` gathers the list of structures to work on,
:class:`PDB <coolpackage.protein.PDB>` downloads and curates them,
:class:`Measure <coolpackage.measure.Measure>` calculates the features,
:class:`Analysis <coolpackage.analysis.Analysis>` explores the results, and
:class:`Viewer <coolpackage.viewer.Viewer>` displays them interactively.

All of them are reachable from the package itself::

    import coolpackage as CPN

    UP = CPN.Uniprot()

They are resolved on first use rather than on import, so that ``import coolpackage``
succeeds on a machine where the optional dependency of one step is missing. The
ImportError of a step is raised when that step is first reached, naming the package that
is absent.
'''

import importlib

__version__ = '1.0'

# public name -> module it lives in
_EXPORTS = {
    'Uniprot': '.uniprot',
    'PDB': '.protein',
    'Measure': '.measure',
    'Analysis': '.analysis',
    'Preprocessing': '.preprocessing',
    'Aggregation': '.aggregation',
    'Viewer': '.viewer',
}

__all__ = sorted(_EXPORTS)


def __getattr__(name):
    '''
    Import the module holding ``name`` the first time the name is looked up.

    See :pep:`562`. The resolved object is cached in the module globals, so the import
    happens once.

    :param name: Attribute looked up on the package.
    :type name: str
    :returns: The class exported under that name.
    :raises AttributeError: if the name is not one of the classes the package exports.
    '''
    if name in _EXPORTS:
        obj = getattr(importlib.import_module(_EXPORTS[name], __name__), name)
        globals()[name] = obj
        return obj
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')


def __dir__():
    '''
    :returns: The names the package exports, including the ones not yet imported.
    :rtype: list
    '''
    return sorted(set(globals()) | set(_EXPORTS))
