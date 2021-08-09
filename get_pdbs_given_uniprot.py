import pandas as pd
import urllib.request, urllib.parse, urllib.error
import re

def get_pdbs_given_uniprot_code(list_UNIPROT_codes):

    columns = ['Uniprot Entry', 'PDB Code', 'Method Structure Obtained by', 'Resolution', 'Chains']
    df = pd.DataFrame(columns=columns)

    for protein_code_clean in list_UNIPROT_codes:
        try:
            url_2 = 'https://www.uniprot.org/uniprot/' + protein_code_clean + '.txt'
            html_2 = urllib.request.urlopen(url_2)
        except Exception as e:
            print('Error %s'%e)
            print('Failed to obtain data for Uniprot entry: ' + protein_code_clean)
            continue

    
        for line in html_2:
            try:
                line = str(line)
                messy_entry = re.findall('PDB; [\w -. ; \d /]*=', line)


                for entry in messy_entry:
                    words = entry.split()
                    chain_ent_num = (len(words) - 1)
                    chain_info = words[chain_ent_num]
                    chain_info = chain_info[:-1]
                    chain_info = chain_info.split('/')
                    PDBCODE = (words[1])[:-1]
                    method_obtained = (words[2])[:-1]
                    resolution = (words[3])[:-1]
                    for i in range(len(chain_info)):
                        if len(chain_info[i]) != 1:
                            chain = (chain_info[i])[:1]
                        else:
                            chain = chain_info[i]
                        data = ({'Uniprot Entry': protein_code_clean, 'PDB Code': PDBCODE, 'Method Structure Obtained by': method_obtained, 'Resolution': resolution, 'Chains': chain})
                        df = df.append(data, ignore_index=True)

                        if len(messy_entry) == 0:
                            AF_code = 'AF-' + protein_code_clean + 'F1-model_v1'
                            data = ({'Uniprot Entry': protein_code_clean, 'PDB Code': AF_code, 'Method Structure Obtained by': 'Predicted', 'Resolution': 'N/A', 'Chains': 'N/A'})
                            df = df.append(data, ignore_index=True)

            except Exception as e:
                print('Error %s'%e)
                print('Failed to obtain data for Uniprot entry ' + protein_code_clean)
    print(df)

    return(df)

if __name__ == "__main__":
    try:    
        list_UNIPROT_codes = ['P50897']
        print(get_pdbs_given_uniprot_code(list_UNIPROT_codes))
    except Exception as e:
        print("ERROR: %s"%e)