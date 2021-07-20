#Gives option of whether you want to scan one/ a few proteins or whole proteome
#Maybe also gives option of input PDB codes or uniprot codes
#Opens uniprot.py
#Opens get_data.py
#Opens propkaparser.py
#Opens SurfAcc.py
#Prints results in nice format

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

    name_of_organism = input('Name of organism:')
    name_of_organism = name_of_organism.replace(' ', '+')
    code = input('Uniprot code:')
    other_uniprot_code = input('Other uniprot code:')
    number_of_chromosomes = input('How many chromosomes does the species have?')
    list_of_pdbs = list()

    starting_number = 0
    list_of_entries = list()


    list_of_pdbs = list()

    number_int = int(number_of_chromosomes)

    n = 0



    #FINDS PDB CODES FOR SINGLE CHROMOSOME ORGANISMS

    if number_of_chromosomes == '1':
        try:
            web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3achromosome&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=0&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
            html = urllib.request.urlopen(web_url)
            soup = BeautifulSoup(html, 'html.parser')
            

            for line in soup:
            
                line = str(line)

                #this regex doesn't work for single chromsome organsisms
                number_of_prot = re.findall('var resultsize = .*;', line)
                
            number_of_prot_2 = number_of_prot[0]
            number_of_prot_2 = int(number_of_prot_2[17:-1])
            print('Number of Uniprot Entries Identified: ', number_of_prot_2)


            while starting_number < number_of_prot_2:

                list_of_entries.append(starting_number)
                starting_number = starting_number + 25


            list_of_entries.append(number_of_prot_2)

                    
            for number in list_of_entries:

                strnum = str(number)
                web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3achromosome&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=' + strnum + '&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
                html = urllib.request.urlopen(web_url)
                list_UNIPROT_codes = list()
                soup = BeautifulSoup(html, 'html.parser')

                for line in soup:
                    line = str(line)
                    stuff_2 = re.findall('id="\w*"><td', line)

                for protein in stuff_2:
                    protein = protein[4:-5]
                    list_UNIPROT_codes.append(protein)

                for proteins in list_UNIPROT_codes:
                        url_2 = 'https://www.uniprot.org/uniprot/' + proteins + '.txt'
                            
                        html_2 = urllib.request.urlopen(url_2)

                        for line in html_2:
                            line = str(line)
                            stuff_2 = re.findall('PDB; ....;', line)
                                

                            for PDBCODE in stuff_2:
                                PDBCODE = PDBCODE[5:-1]
                                list_of_pdbs.append(PDBCODE)

                                                
                n = n + 1
                list_of_pdbs = list(dict.fromkeys(list_of_pdbs))
                print ("Number of PDB structures obtained: ", len(list_of_pdbs))
                print('Next Page', str(n))
                number = number + 25

        except:
            print('SEARCH FAILED')







    #FINDS PDB CODES FOR MULTI ORGANISM SPECIES

    number_int = int(number_of_chromosomes)

    if number_of_chromosomes != '1':
        web_url = 'https://www.uniprot.org/proteomes/' + code
        html = urllib.request.urlopen(web_url)
        x = 0
        soup = BeautifulSoup(html, 'html.parser')
        table = soup.find_all('table')
        for line in table:
            line = str(line)
            ID = re.findall('name="[\w ]*" type', line)

        IDS = list()

        for chromosome in ID:
            chromosome = chromosome[6:-6]
            chromosome = chromosome.replace(' ', '+')
            IDS.append(chromosome)

        try:
            for chromosome in IDS:
                print(chromosome)
                chromosome = chromosome.lower()
                web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a%22' + chromosome + '%22&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=0&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
                html = urllib.request.urlopen(web_url)
                soup = BeautifulSoup(html, 'html.parser')
                print('URL good')

                for line in soup:
                    line = str(line)
                    number_of_prot = re.findall('var resultsize = .*;', line)
                number_of_prot_2 = number_of_prot[0]
                number_of_prot_2 = int(number_of_prot_2[17:-1])
                print('Number of Uniprot Entries Identified: ', number_of_prot_2)

                while starting_number < number_of_prot_2:

                    list_of_entries.append(starting_number)
                    starting_number = starting_number + 25


                list_of_entries.append(number_of_prot_2)

                    
                for number in list_of_entries:

                        strnum = str(number)
                        web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a%22' + chromosome + '%22&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=' + strnum + '&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
                        html = urllib.request.urlopen(web_url)
                        list_UNIPROT_codes = list()
                        soup = BeautifulSoup(html, 'html.parser')

                        for line in soup:
                            line = str(line)
                            stuff_2 = re.findall('id="\w*"><td', line)

                        for protein in stuff_2:
                            protein = protein[4:-5]
                            list_UNIPROT_codes.append(protein)

                        for proteins in list_UNIPROT_codes:
                                url_2 = 'https://www.uniprot.org/uniprot/' + proteins + '.txt'
                                
                                html_2 = urllib.request.urlopen(url_2)

                                for line in html_2:
                                    line = str(line)
                                    stuff_2 = re.findall('PDB; ....;', line)
                                    

                                    for PDBCODE in stuff_2:
                                        PDBCODE = PDBCODE[5:-1]
                                        list_of_pdbs.append(PDBCODE)


                                                    

                        list_of_pdbs = list(dict.fromkeys(list_of_pdbs))
                        print ("Number of PDB structures obtained: ", len(list_of_pdbs))

                        n = n + 1
                        print('Next Page', str(n))
                        number = number + 25

        except:
            print('SEARCH FAILED')

###############################################

elif answers["Choice"] == 'Single Protein':
    list_of_pdbs.append(input('PDB code:'))
    
df = pd.DataFrame(list_of_pdbs)
df.to_csv('pdb_codes.csv', index=False)



##################################
#Need to add in get_data.py here
##################################

#Code from here onwards finds the pKa of lysine epsilon side chains then determines the surface accesiblility of lysines with pKa < 9
#Need to change the .csv file to whatever comes out of get_data.py
file = open('pdb_codes.csv')
csv_reader = csv.reader(file)
next(csv_reader)

for row in csv_reader:
    pdb = row[0]

    
    try:
        pdb_code = pdb + '.pdb'

        process = subprocess.Popen(['python', '-m', 'propka', pdb_code],
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        lysines = list()
        pkas = list()
        pkafile = pdb + '.pka'
        propres = open(pkafile)
        chain = list()

        for line in propres:
            if re.search('^   LYS' , line):
                line = line.strip()

                res = re.findall('LYS(.* [A-Z])', line)
                chain.append((res[0])[-1:])
                res = (res[0])[:-1]

                
                if res[0] == ' ':
                    res = res[1:]

                lysines.append(res)

                pka = line[13:18]
                pkas.append(pka)

        df = pd.DataFrame({'resid':lysines, 'chain':chain, 'pKa':pkas})
        df.sort_values(by=['pKa'], inplace=True)
        convert_dict = {'pKa': float}
        df = df.astype(convert_dict)
        convert_dict = {'chain': float}
        df = df.where(df['pKa'] < 9)
        df = df.dropna()
        if df.empty:
            print('No lysine residues with an epsilon amino group pKa < 9 were found in ' + pdb)
            break
        else:
            print('The following lysine residues have epsilon amino group pKa value(s) < 9')
            print(df)




        M = bb.Molecule()
        M.import_pdb(pdb_code)
        df_2 = M.data


        resid_list = df['resid'].tolist()


        for ResidueID in resid_list:
            ResidueID = str(ResidueID)
            pts, indices = M.atomselect("*",  ResidueID, ["CA"],  use_resname=False, get_index=True)


            x = sasa(M, targets=indices, probe=1.4, n_sphere_point=960, threshold=0.05)
            print(x[0])

    except:
        print(pdb_code, 'Failed testing')



    


