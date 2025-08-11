'''
File for generalised reporting of errors from feature calculations.
Errors are written out into the text file 'measures_errors_{date}_{x}.txt' where
date is the date which the measures has been run and x is the run number that day.
'''

def _report_error_to_file(measurement_stage, path, error, error_filename='measure_errors.txt'):
    '''
    Helper function to remove redundant code writing errors in the measurements to the
    measurement error log file.
    
    Parameters
    ----------
    measurement_stage -> string
        The stage of measurements that has caused the error with the file, eg propka 1
    path -> string
        The path of the pdb file that the measurement has been attempted on
    error -> string
        The error that has been produced at that step of the measurement when it has been
        attempted to extract features from the pdb file
    error_filename -> string
        The name of the file to write the errors to. If the default name is used, then
        nothing will be written to the file. Default name is 'measure_errors.txt', custom
        names are of the format 'measure_errors_{date}_{x}.txt' where date is the date
        that the measures have been run and x is the run number that day.
    
    Example
    -------
    self._report_error_to_file('propka 1', path, e)
    '''
    if error_filename != 'measure_errors.txt':
        with open(error_filename, 'a') as e_f:
            e_f.writelines('-----------------------------------------------------------------------\n')
            e_f.writelines(f'{measurement_stage} calc error\n')
            e_f.writelines(path + '\n')
            e_f.writelines(error + '\n')
