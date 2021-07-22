import inquirer
import urllib.request, urllib.parse, urllib.error
import re
from bs4 import BeautifulSoup
import csv
import pandas as pd
import os
import subprocess
import numpy as np

def get_pdbs(answer_2, name_of_organism, code, other_uniprot_code, number_of_chromosomes):

    list_of_entries = list()
    list_of_pdbs = list()
    starting_number = 0
    number_int = int(number_of_chromosomes)
    n = 0

    if number_of_chromosomes == '1':
        try:
            if answer_2 == 'Yes':
                chromosome = 'genome'
            else:
                chromosome = 'chromosome'
            web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a' + chromosome + '&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=0&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
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
                web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a' + chromosome + '&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=' + strnum + '&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
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
                        #Need to fix end.

        except:
            print('SEARCH FAILED')

    df = pd.DataFrame(list_of_pdbs)
    df.to_csv('pdb_codes.csv', index=False)

    return()

if __name__ == "__main__":
    try:    
        get_pdbs('No', 'Homo+sapiens+(Human)', 'UP000005640', '5b9606', '26')
    except:
        print('Error')