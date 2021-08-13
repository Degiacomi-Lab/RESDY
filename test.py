
import get_initial_data as gd
import pandas as pd

#This code parses a .csv file to find uniprot and pdb codes to pass into the pipeline.
#If you want to just input uniprot codes put them in the first column and leave the second empty
#If you want to input pdb codes put the uniprot code in the first column and pdb code in the second.

def from_csv_file(csv_file, results_df):
#First read .csv file.
    try:
        csv_df = pd.read_csv(csv_file)
        print('.csv file successfully opened')
    except Exception as e:
        print('Error: %s'%e)
        print('Failed to find .csv file.')
#Next abstract column names
    try:
        column_names = list(csv_df.columns)
        print(column_names)
        csv_df[column_names[1]] = csv_df[column_names[1]].fillna(0)

    except Exception as e:
        print('Error: %s'%e)
        print('Failed to get data from .csv file')
#Lastly feed them into the functions which get the data about the protein and append it to the results_df dataframe.

    for i in range(len(csv_df)):
        try:
            uniprot_code = csv_df.at[i, column_names[0]]
            pdb_code = csv_df.at[i, column_names[1]]
            if pdb_code == 0:
                results_df = get_pdbs_uniprot(uniprot_code, results_df)
            else:
                results_df = construct_single_pdb_df(uniprot_code, pdb_code, results_df)
        except Exception as e:
            print('Error %s'%e)
            continue

    return(results_df)

if __name__ == "__main__":

    try:
        columns = ['Uniprot Entry', 'PDB Code', 'Method Structure Obtained by', 'Resolution', 'Chains']
        results_df = pd.DataFrame(columns=columns)
        print(from_csv_file('test.csv', results_df))
        
    except Exception as e:
                print("ERROR: %s"%e)
