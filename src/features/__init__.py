'''
Initialise the list of features to be used in the measurements
'''

import os
import importlib
import inspect

__all__ = []
# resolved from __file__ rather than from the working directory, so that the package can be
# imported from anywhere (Sphinx, for instance, builds from the docs folder)
feature_folder = os.path.dirname(os.path.abspath(__file__))
scripts_to_ignore = ['feature.py', '__init__.py', 'error_reporting.py']
feat_folder_scripts = [a for a in os.listdir(feature_folder) if a not in scripts_to_ignore
                       and a.endswith('.py')
                       and (os.path.isfile(os.path.join(feature_folder, a)))]

for f in feat_folder_scripts:
    mod_name = f'{__name__}.{f.split(".")[0]}'
    module = importlib.import_module(mod_name)

    for name, obj in inspect.getmembers(module, inspect.isclass):
        if obj.__module__ == mod_name:
            globals()[name] = obj
            __all__.append(name)
