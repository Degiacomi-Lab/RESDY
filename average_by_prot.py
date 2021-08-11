import pandas as pd
import os

#This code calculates the mean and standard deviation for the sasa and pKa values from all the structures of a given protein.

def average_prot(all_pka_sasa_res):
    print('Averaging Uniprot data...')

#First creates a dataframe and puts all uniprot entries in a list.

    try:
        columns = ['resid', 'chain', 'plddt', 'Uniprot Entry', 'pKa', 'pKa stdev', 'PDB Code', 'sasa', 'sasa stdev']
        avgd_pka_sasa = pd.DataFrame(columns=columns)
        list_of_uniprot_codes = all_pka_sasa_res['Uniprot Entry']
        list_of_uniprot_codes_no_dup = list()
    except Exception as e:
        print('Failed to construct dataframe')
        print("ERROR: %s"%e)
        return()

#Next removes all duplicates from the list of uniprot entries.

    try:
        for entry in list_of_uniprot_codes:
            if entry not in list_of_uniprot_codes_no_dup:
                list_of_uniprot_codes_no_dup.append(entry)
            else:
                continue

    except Exception as e:
        print('Failed to obtain Uniprot Entries from dataframe')
        print("ERROR: %s"%e)
        return()
#Next constructs a df for all the lysines of all the structures of a given uniprot entry.

    for uniprot_code in list_of_uniprot_codes_no_dup:

        try:
            df_one_uniprot_code = all_pka_sasa_res.where(all_pka_sasa_res['Uniprot Entry'] == uniprot_code)
            df_one_uniprot_code = df_one_uniprot_code[df_one_uniprot_code['Uniprot Entry'].notna()]
            list_chains = df_one_uniprot_code['chain']
            list_chains_no_dup = list()

#Next makes a list of all the chains in the structure

            for entry in list_chains:
                if entry not in list_chains_no_dup:
                    list_chains_no_dup.append(entry)
                else:
                    continue

        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return()


#Next makes a df including all the lysines in a specific chain in all uniprot entries.

        for chain in list_chains_no_dup:

            try:
                df_one_chain = df_one_uniprot_code.where(df_one_uniprot_code['chain'] == chain)
                df_one_chain= df_one_chain[df_one_chain['chain'].notna()]

#Next makes a list of the lysine resids in that chain and gets rid of all duplicates.

                list_of_resids = df_one_chain['resid']
                list_of_resids_no_dup = list()
                for entry in list_of_resids:
                    if entry not in list_of_resids_no_dup:
                        list_of_resids_no_dup.append(entry)
                    else:
                        continue
                    
            except Exception as e:
                print('Failed to construct dataframe for chain ' + chain)
                print("ERROR: %s"%e)
                continue
#Next calculates the mean and standard deviation for each individual lysine's pKa/sasa in all the available structures and puts them in a df.

            for residue in list_of_resids_no_dup:
                try:
                    df_one_resid = df_one_chain.where(df_one_chain['resid'] == residue)
                    df_one_resid = df_one_resid[df_one_resid['resid'].notna()]
                except Exception as e:
                    print('Failed to construct dataframe for residue ' + residue)
                    print("ERROR: %s"%e)
                    continue

                try:
                    df_one_resid["sasa"] = pd.to_numeric(df_one_resid["sasa"], downcast="float")

                    avg_pka = df_one_resid['pKa'].mean()
                    avg_sasa = df_one_resid['sasa'].mean()
                    stddev_pka = df_one_resid['pKa'].std()
                    stdev_sasa = df_one_resid['sasa'].std()

                except Exception as e:
                    print('Failed to average data for ' + residue)
                    print("ERROR: %s"%e)
                    continue

                try:
                    avg_plddt = df_one_resid['plddt'].mean()
                except:
                    avg_plddt = 'N/A'

                try:
                    list_of_pdbs = df_one_resid['PDB Code']
                    list_of_pdbs_no_dup = list()

                    for entry in list_of_pdbs:
                        if entry not in list_of_pdbs_no_dup:
                            list_of_pdbs_no_dup.append(entry)
                        else:
                            continue

                    PDB_codes_avgd = ''
                    for i in range(len(list_of_pdbs_no_dup)):
                        PDB_codes_avgd = PDB_codes_avgd + '/' + list_of_pdbs_no_dup[i]
    
                    d = {'resid': residue, 'chain': chain, 'plddt': avg_plddt, 'Uniprot Entry': uniprot_code, 'pKa': avg_pka, 'pKa stdev':stddev_pka, 'PDB Code': PDB_codes_avgd, 'sasa': avg_sasa, 'sasa stdev': stdev_sasa}
                    avgd_pka_sasa = avgd_pka_sasa.append(d, ignore_index=True)
    
            
                except Exception as e:
                    print('Failed to construct final dataframe for residue ' + residue)
                    print("ERROR: %s"%e)
                    continue
                
    #Lastly appends into df and saves as Output/results.csv

        try:
            avgd_pka_sasa['sasa stdev'] = avgd_pka_sasa['sasa stdev'].fillna(0)
            avgd_pka_sasa['pKa stdev'] = avgd_pka_sasa['pKa stdev'].fillna(0)

            if not os.path.exists("Output"):
                os.mkdir("Output")

            avgd_pka_sasa.to_csv('Output/results.csv')
            
            print('Done')

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed to average data')
            return()


    return(avgd_pka_sasa)






if __name__ == "__main__":
    try:
        columns = ['Uniprot Entry', 'resid', 'chain', 'pKa', 'sasa', 'PDB Code']
        all_pka_sasa_res = pd.DataFrame(columns=columns)
        data = {'Uniprot Entry': 'ABC', 'resid': 4, 'chain': 'A', 'pKa': 7.3, 'sasa': 4, 'PDB Code': '4jjf'}
        all_pka_sasa_res = all_pka_sasa_res.append(data, ignore_index=True)
        data_2 = {'Uniprot Entry': 'ABC', 'resid': 4, 'chain': 'A', 'pKa': 7.4, 'sasa': 5, 'PDB Code': '6jig'}
        data_3 = {'Uniprot Entry': 'CBA', 'resid': 2, 'chain': 'B', 'pKa': 7.2, 'sasa': 1, 'PDB Code': '5jof'}
        data_4 = {'Uniprot Entry': 'CBA', 'resid': 2, 'chain': 'B', 'pKa': 9.3, 'sasa': 2, 'PDB Code': '9ifr'}
        all_pka_sasa_res = all_pka_sasa_res.append(data_2, ignore_index = True)
        all_pka_sasa_res = all_pka_sasa_res.append(data_3, ignore_index = True)
        all_pka_sasa_res = all_pka_sasa_res.append(data_4, ignore_index = True)
        print(all_pka_sasa_res)
        #Problem must be that df coming in is wrong datatype...
        print(average_prot(all_pka_sasa_res))
    except Exception as e:
        print("ERROR: %s"%e)
