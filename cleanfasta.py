from os import remove
from re import M, X
import biobox as bb
import pandas as pd
import fileinput
import re
import os
import subprocess





#This cleans the fasta file, getting rid of any non-cannonical AAs that may cause an issue.
def clean_fasta(new_name):

#Firstly it opens the fasta file

    path = 'curate_PDB/raw/' + new_name

#Next it rewrites it replacing the amino acids that may cause an issue.

    with fileinput.FileInput(path, inplace = True) as f:
        for line in f:
            if("KCX" in line):
                line = line.replace('(KCX)', 'K')
                print(line, end = '')
            if("MSE" in line):
                line = line.replace('(MSE)', 'M')
                print(line, end = '')
            else:
                print(line, end = '')

    return
