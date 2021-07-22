import pandas as pd

def get_single_pdb(PDBCODE):
        list_of_pdbs = list()
        list_of_pdbs.append(PDBCODE)
        df = pd.DataFrame(list_of_pdbs)
        df.to_csv('pdb_codes.csv', index=False)

        return()

#FOR TESTING
if __name__ == "__main__":
    try:
        get_single_pdb('ABC')
    except:
        print('ERROR')