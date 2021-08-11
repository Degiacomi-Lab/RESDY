import urllib.request, urllib.parse, urllib.error
import pandas as pd
import re
import subprocess
import biobox as bb
import os
from get_data_from_uniprot import get_chains

#This module is called if the user wants to investigate a single pdb code.
#It obtains the resolution/method obtained data from uniprot and the chain info from the PDB (via the get_chains module)

def construct_single_pdb_df(UNIPROT_code_pdb, PDBCODE_inpt, results_df):

#Firstly it opens the correct uniprot page.
    try:
        url_2 = 'https://www.uniprot.org/uniprot/' + UNIPROT_code_pdb + '.txt'
        html_2 = urllib.request.urlopen(url_2)

    except Exception as e:
        print('Error %s'%e)
        print('Failed to obtain uniprot entry for uniprot code: ' + UNIPROT_code_pdb)
        return()

#Next it parses and looks for PDB entries
    for line in html_2:
        try:
            line = str(line)
            messy_entry = re.findall('PDB; [\w -. ; \d /]*=', line)
            for entry in messy_entry:
                try:
                    words = entry.split()
                    PDBCODE = (words[1])[:-1]

#It only obtains the data from uniprot if the PDB code matches that which was input by the user.

                    if PDBCODE_inpt == PDBCODE:
                        unique_values_chain = get_chains(PDBCODE_inpt)
                        method_obtained = (words[2])[:-1]
                        resolution = (words[3])[:-1]
                        for i in range(len(unique_values_chain)):
                            try:
                                data = ({'Uniprot Entry': UNIPROT_code_pdb, 'PDB Code': PDBCODE, 'Method Structure Obtained by': method_obtained, 'Resolution': resolution, 'Chains': unique_values_chain[i]})
                                results_df = results_df.append(data, ignore_index=True)
                            except Exception as e:
                                print("Error %s"%e)
                                print('Error constructing df')
                                continue
                except:
                    continue


        except Exception as e:
            print("Error %s"%e)
            return()

    if (len(results_df) == 0):
        print('FAILURE- likely due to PDB code not beloning to uniprot entry.')

    
    return(results_df)