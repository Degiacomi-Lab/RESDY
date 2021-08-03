from os import name
import numpy
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import inquirer

def analyse_data(all_pka_sasa_res, carbam_pdb_list, carbam_resid_list):
    print(carbam_pdb_list)
    print(carbam_resid_list)
    try:
        list_of_ids = list()
        all_pka_sasa_res = all_pka_sasa_res.astype({"resid": int})

    except Exception as e:
        print("ERROR: %s"%e)

    
    for i in range(len(carbam_resid_list)):

        try:
            print(type(carbam_pdb_list[i]))
            print(type(carbam_resid_list[i]))
            carbam_res = all_pka_sasa_res.where(all_pka_sasa_res['resid'] == int(carbam_resid_list[i]))
            carbam_res = carbam_res.where(carbam_res['PDB Code'] == carbam_pdb_list[i])

            carbam_res_df = carbam_res[carbam_res['resid'].notna()]
            list_ids = carbam_res_df.index.tolist()

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed to mark carbamate on resid = ' + str(i))

    try:
        list_of_ids = list_of_ids + list_ids
    except Exception as e:
            print("ERROR: %s"%e)

    try:
        all_pka_sasa_res['carbamylated'] = all_pka_sasa_res.index.isin(list_of_ids)
        print(all_pka_sasa_res)
    except Exception as e:
            print("ERROR: %s"%e)

    try:
        fig, ax = plt.subplots()
        colors = {True:'red', False:'black'}
        ax.scatter(all_pka_sasa_res['sasa'], all_pka_sasa_res['pKa'], c=all_pka_sasa_res['carbamylated'].map(colors))
        plt.title('pKa vs sasa')
        plt.xlabel('sasa')
        plt.ylabel('pKa')
        #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
        #ax.set_xticks(ax.get_xticks()[::4])
        plt.show()

    except Exception as e:
            print("ERROR: %s"%e)
            print('Error plotting data.')

    return(all_pka_sasa_res)





if __name__ == "__main__":
    try:
        columns = ['Uniprot Code', 'resid', 'chain', 'pKa', 'sasa', 'PDB Code']
        all_pka_sasa_res = pd.DataFrame(columns=columns)
        data = {'Uniprot Code': 'ABC', 'resid': '4', 'chain': 'A', 'pKa': 7.3, 'sasa': 4, 'PDB Code': '4jjf'}
        all_pka_sasa_res = all_pka_sasa_res.append(data, ignore_index=True)
        data_2 = {'Uniprot Code': 'ABC', 'resid': '4', 'chain': 'A', 'pKa': 7.4, 'sasa': 5, 'PDB Code': '6jig'}
        all_pka_sasa_res = all_pka_sasa_res.append(data_2, ignore_index = True)
        #Problem must be that df coming in is wrong datatype...
        carbam_pdb_list = ['4jjf']
        carbam_resid_list = [4]
        print(analyse_data(all_pka_sasa_res, carbam_pdb_list, carbam_resid_list))
        
    except Exception as e:
        print("ERROR: %s"%e)