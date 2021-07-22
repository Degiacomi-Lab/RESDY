import inquirer
import urllib.request, urllib.parse, urllib.error
import re
from bs4 import BeautifulSoup
import csv
import pandas as pd
import os
import subprocess
import numpy as np
import biobox as bb
from biobox.measures.calculators import sasa

from single_pdb import get_single_pdb
from uniprot import get_pdbs
from pka_sasa import calculate_pKa_and_SASA


list_of_pdbs = list()

questions = [
  inquirer.List('Choice',
                message="Do you want to search a whole organism's proteome or a single protein?",
                choices=['Whole Proteome', 'Single Protein'],
            ),
]
answers = inquirer.prompt(questions)

if answers["Choice"] == 'Whole Proteome':

############NOT BEST WAY OF DOING THIS PROBABLY##############
    question_2 = [
  inquirer.List('Choice',
                message="Is the organism a virus?",
                choices=['Yes', 'No'],
            ),
]
    answer_2 = inquirer.prompt(question_2)
    answer_2 = answer_2['Choice']
    
    name_of_organism = input('Name of organism:')
    name_of_organism = name_of_organism.replace(' ', '+')
    code = input('Uniprot code:')
    other_uniprot_code = input('Other uniprot code:')
    number_of_chromosomes = input('How many chromosomes does the species have?')

    get_pdbs(answer_2, name_of_organism, code, other_uniprot_code, number_of_chromosomes)

elif answers ['Choice'] == 'Single Protein':
    PDBCODE = input('What is the PDB code?')
    get_single_pdb(PDBCODE)


#DOWNLOADS PDB FILES

file = open('pdb_codes.csv')
csv_reader = csv.reader(file)
next(csv_reader)

for row in csv_reader:
    pdb = row[0]
    subprocess.check_call("wget https://files.rcsb.org/download/" + pdb + ".pdb", shell=True)

calculate_pKa_and_SASA(pdb)
