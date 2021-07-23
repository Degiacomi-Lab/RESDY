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
from uniprot import get_pdbs
from pka_sasa import calculate_pKa_and_SASA
skip = 0
running = True
while running:
    try:
        list_of_pdbs = list()
        print('Hello, press Ctrl + C at anytime to return to the start')
        questions = [
        inquirer.List('Choice',
                        message="Do you want to search a whole organism's proteome or a single protein?",
                        choices=['Whole Proteome', 'Single Protein', 'Quit'],
                    ),
        ]
        answers = inquirer.prompt(questions)

        if answers["Choice"] == 'Whole Proteome':


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

            try:
                print('Ctrl + c at anytime to return to start')
                list_of_pdbs = get_pdbs(answer_2, name_of_organism, code, other_uniprot_code, number_of_chromosomes)
            except:
                print('Try again')
                skip = 1
            if len(list_of_pdbs) == 0:
                skip = 1
            

        elif answers ['Choice'] == 'Single Protein':
            PDBCODE = input('What is the PDB code?')
            list_of_pdbs.append(PDBCODE)

        elif answers ['Choice'] == 'Quit':
            skip = 1
            raise Exception

        number_of_pdbs = str(len(list_of_pdbs))

        #DOWNLOADS PDB FILES

        if skip == 0:
            for pdb in list_of_pdbs:
                subprocess.check_call("wget https://files.rcsb.org/download/" + pdb + ".pdb", shell=True)

            calculate_pKa_and_SASA(pdb)

    except KeyboardInterrupt:
            pass

    except:
        if skip == 1:
            print('Goodbye')
            running = False