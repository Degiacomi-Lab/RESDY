from os import remove
from re import X
import biobox as bb
import pandas as pd
import fileinput
import re
import os
import subprocess
import split as sp





def remove_kcx_from_fasta(pdb_code, chain):
    path = 'curate_PDB/raw/' + pdb_code + '_' + chain + '.fasta'
    with fileinput.FileInput(path, inplace = True) as f:
        for line in f:
            if("KCX" in line):
                line = line.replace('(KCX)', 'K')
                print(line, end ='')
            else:
                print(line, end ='')
    return()

if __name__ == "__main__":
    try:
        print(clean('4E3T', 'A'))
    except Exception as e:
        print("ERROR: %s"%e)


#def clean_hetatm()