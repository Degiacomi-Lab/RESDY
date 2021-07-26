import inquirer
import urllib.request, urllib.parse, urllib.error
import re
from bs4 import BeautifulSoup
import csv
import pandas as pd
import os
import subprocess
import numpy as np
import sys


def get_pdbs(name_of_organism, code):
    list_of_entries = list()
    list_of_pdbs = list()
    list_clean = list()
    starting_number = 0
    n = 0
    end = 0
    try:
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

    except Exception as e:
        print("ERROR: %s"%e)
        print('Uniprot code is invalid')
        #Returns to the start

    for chromosome in IDS:
        try:
            print(chromosome)
            chromosome = chromosome.lower()
            web_url = 'https://www.uniprot.org/uniprot/?query=proteome:' + code + '+AND+proteomecomponent:%22' + chromosome + '%22&sort=score'
            html = html = urllib.request.urlopen(web_url)
            soup = BeautifulSoup(html, 'html.parser')

            for line in soup:
                line = str(line)
                code_other = re.findall('%[\w ()\d -]*%5d%22', line)

            for other_uniprot_code in code_other:
                other_uniprot_code = (code_other[0])[1:-6]

            if len(IDS) > 1:
                web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a%22' + chromosome + '%22&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=0&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
                html = urllib.request.urlopen(web_url)
                soup = BeautifulSoup(html, 'html.parser')

            if len(IDS) == 1:
                web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a' + chromosome + '&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=0&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
                html = urllib.request.urlopen(web_url)
                soup = BeautifulSoup(html, 'html.parser')

            for line in soup:
                line = str(line)
                list_of_species = re.findall('>[\w ()\d -]*</a></td><td class="number', line)
                for species in list_of_species:
                    species = species[1:-26]
                    species = species.replace(' ', '+')
                    list_clean.append(species)


            if list_clean[0] != name_of_organism:
                end = 1
                raise Exception('Likely due to incorrect information added')
            if list_clean[1] != name_of_organism:
                end = 1
                raise Exception('Likely due to incorrect information added')
            if list_clean[2] != name_of_organism:
                end = 1
                raise Exception('Likely due to incorrect information added')


            for line in soup:
                line = str(line)
                number_of_prot = re.findall('var resultsize = .*;', line)
            number_of_prot_2 = number_of_prot[0]
            number_of_prot_2 = int(number_of_prot_2[17:-1])
            print('Number of Uniprot Entries Identified: ', number_of_prot_2)
        

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed obtaining Uniprot codes for ' + chromosome)
            if end == 1:
                break
            else:
                continue

        while starting_number < number_of_prot_2:

            list_of_entries.append(starting_number)
            starting_number = starting_number + 25


        list_of_entries.append(number_of_prot_2)

                            
        for number in list_of_entries:
            try:
                strnum = str(number)
                web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a%22' + chromosome + '%22&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=' + strnum + '&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
                html = urllib.request.urlopen(web_url)
                list_UNIPROT_codes = list()
                soup = BeautifulSoup(html, 'html.parser')

                for line in soup:
                    line = str(line)
                    messy_list_of_uniprot_codes = re.findall('id="\w*"><td', line)

                for messy_protein_code in messy_list_of_uniprot_codes:
                    protein_code = messy_protein_code[4:-5]
                    list_UNIPROT_codes.append(protein_code)

            except Exception as e:
                print("ERROR: %s"%e)
                print('Failed to obtain Uniprot codes on:' + web_url)
                continue
                
            if len(list_UNIPROT_codes) == 0:
                print('No Uniprot codes were found on: ' + web_url)
                continue

            for protein_code_clean in list_UNIPROT_codes:
                try:
                    url_2 = 'https://www.uniprot.org/uniprot/' + protein_code_clean + '.txt'
                                            
                    html_2 = urllib.request.urlopen(url_2)

                    for line in html_2:
                        line = str(line)
                        messy_pdb_codes = re.findall('PDB; ....;', line)
                                                

                        for PDBCODE in messy_pdb_codes:
                            PDBCODE = PDBCODE[5:-1]
                            list_of_pdbs.append(PDBCODE)

                except Exception as e:
                    print("ERROR: %s"%e)
                    print('Failed to obtain PDB codes for' + protein_code_clean)

                                                            

            list_of_pdbs = list(dict.fromkeys(list_of_pdbs))
            print ("Number of PDB structures obtained: ", len(list_of_pdbs))

            number = number + 25
            print('Number of Uniprot Entries Searched: ' + str(number) + '\n')
            if number > number_of_prot_2:
                number = number_of_prot_2
            print((str((number/number_of_prot_2)*100))[0:3] + '%')

    return list_of_pdbs

if __name__ == "__main__":
    try:    
        print(get_pdbs('Severe+acute+respiratory+syndrome+coronavirus+2+(2019-nCoV)+(SARS-CoV-2)', 'UP000464024'))
    except Exception as e:
                print("ERROR: %s"%e)