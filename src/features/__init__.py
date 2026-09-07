'''
Initialise the list of features to be used in the measurements
'''

import os
import importlib
import inspect

__all__ = []
feature_folder = f'src{os.sep}features'
scripts_to_ignore = ['feature.py', '__init__.py']
feat_folder_scripts = [a for a in os.listdir(feature_folder) if a not in scripts_to_ignore
                       and (os.path.isfile(os.path.join(f'src{os.sep}features', a)))]

for f in feat_folder_scripts:
    mod_name = f'{__name__}.{f.split(".")[0]}'
    module = importlib.import_module(mod_name)

    for name, obj in inspect.getmembers(module, inspect.isclass):
        if obj.__module__ == mod_name:
            globals()[name] = obj
            print(name, globals()[name])
            __all__.append(name)
