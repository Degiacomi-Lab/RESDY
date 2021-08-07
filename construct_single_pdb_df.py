import urllib.request, urllib.parse, urllib.error
import pandas as pd
import re
import subprocess
import biobox as bb
import os
from get_chains import get_chains

def construct_single_pdb_df(UNIPROT_code_pdb, PDBCODE_inpt, results_df):
    set_up_df = True
    try:
        url_2 = 'https://www.uniprot.org/uniprot/' + UNIPROT_code_pdb + '.txt'
        html_2 = urllib.request.urlopen(url_2)



    except Exception as e:
        print('Error %s'%e)
        print('Failed to obtain uniprot entry for uniprot code: ' + UNIPROT_code_pdb)
        set_up_df = False

    if set_up_df == True:
        for line in html_2:
            line = str(line)
            messy_entry = re.findall('PDB; [\w -. ; \d /]*=', line)
            for entry in messy_entry:
                try:
                    words = entry.split()
                    PDBCODE = (words[1])[:-1]
                    if PDBCODE_inpt == PDBCODE:
                        unique_values_chain = get_chains(PDBCODE_inpt)
                        method_obtained = (words[2])[:-1]
                        resolution = (words[3])[:-1]
                        for i in range(len(unique_values_chain)):
                            data = ({'Uniprot Entry': UNIPROT_code_pdb, 'PDB Code': PDBCODE, 'Method Structure Obtained by': method_obtained, 'Resolution': resolution, 'Chains': unique_values_chain[i]})
                            results_df = results_df.append(data, ignore_index=True)
                except:
                    continue

    if (len(results_df) == 0):
        print('FAILURE- likely due to PDB code not beloning to uniprot entry.')

    
    return(results_df)